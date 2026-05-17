from __future__ import annotations

import json
import math
import re
from pathlib import Path

import cv2
import numpy as np
from PIL import Image
from rapidocr_onnxruntime import RapidOCR


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "_ocr_editable_layout.json"
ASSET_DIR = ROOT / "_editable_assets"
FULL_BACKGROUND_PAGES = {1}
FULL_BACKGROUND_EDITABLE_TEXT = {
    1: {"AI 辅助科研实践", "以LBM热对流算法与并行优化为例", "课堂学术汇报"},
}
MANUAL_CLEANUPS = {
    21: [
        {"left": 882.0, "top": 452.0, "width": 22.0, "height": 26.0, "bg": "#FFFFFF"},
    ],
}
MANUAL_CROPS = {
    2: [
        {"name": "card1_number", "rect": (356, 243, 440, 327)},
        {"name": "card1_icon", "rect": (466, 241, 566, 339)},
        {"name": "card2_number", "rect": (996, 243, 1081, 327)},
        {"name": "card2_icon", "rect": (1103, 242, 1203, 342)},
        {"name": "card3_number", "rect": (356, 533, 440, 618)},
        {"name": "card3_icon", "rect": (466, 530, 564, 632)},
        {"name": "card4_number", "rect": (996, 533, 1082, 618)},
        {"name": "card4_icon", "rect": (1109, 530, 1201, 634)},
    ],
}
MANUAL_FRAME_CROPS = {
    2: [
        {"name": "card1_frame", "rect": (326, 211, 951, 488), "erase": [(356, 243, 440, 327), (466, 241, 566, 339)]},
        {"name": "card2_frame", "rect": (967, 211, 1630, 488), "erase": [(996, 243, 1081, 327), (1103, 242, 1203, 342)]},
        {"name": "card3_frame", "rect": (326, 502, 951, 780), "erase": [(356, 533, 440, 618), (466, 530, 564, 632)]},
        {"name": "card4_frame", "rect": (967, 502, 1630, 780), "erase": [(996, 533, 1082, 618), (1109, 530, 1201, 634)]},
    ],
}
NATIVE_SHAPES = {}
SLIDE_W_PT = 960.0
SLIDE_H_PT = 540.0
IMG_W = 1672.0
IMG_H = 941.0
SX = SLIDE_W_PT / IMG_W
SY = SLIDE_H_PT / IMG_H


NOISE_WORDS = (
    "BELIKANG",
    "BEIHANG",
    "NORTHWESTERN",
    "POLYTECHNICAL",
    "UNIVERSITY",
    "1938",
)


def rgb_to_hex(rgb: tuple[int, int, int]) -> str:
    return "#{:02X}{:02X}{:02X}".format(*rgb)


def luminance(rgb: np.ndarray | tuple[int, int, int]) -> float:
    r, g, b = [float(x) for x in rgb]
    return 0.299 * r + 0.587 * g + 0.114 * b


def color_distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a.astype(float) - b.astype(float)))


def normalize_text(text: str) -> str:
    text = text.strip()
    text = text.replace("／", "/")
    text = re.sub(r"\s+", " ", text)
    if len(text) >= 5:
        text = text.strip("一—-＿_ ")
    text = text.replace("AI提效", "AI 提效")
    text = text.replace("AI不是", "AI 不是")
    text = text.replace("AI辅助科研", "AI 辅助科研")
    text = text.replace("利用GPU", "利用 GPU")
    text = text.replace("多节点/多", "多节点 / 多")
    text = text.replace("CPU/GPU/MPI", "CPU / GPU / MPI")
    replacements = {
        "以 LBM": "以LBM",
        "GPU加速": "GPU 加速",
        "CPU 多核": "CPU 多核",
        "AI 不是让科研变得不需要思考，而是让我们把更多时间留给真正需要思考的问题": "AI 不是让科研变得不需要思考，而是让我们把更多时间留给真正需要思考的问题。",
        "OpenACC：低重构成本移植 GPU": "OpenACC：低重构成本移植 GPU",
        "MPI：多节点/多 GPU 扩展": "MPI：多节点 / 多 GPU 扩展",
    }
    return replacements.get(text, text)


def should_skip(page: int, text: str, score: float, box: list[list[float]]) -> bool:
    if not text:
        return True

    upper = text.upper()
    if any(word in upper for word in NOISE_WORDS):
        return True

    xs = [p[0] for p in box]
    ys = [p[1] for p in box]
    x1, x2 = min(xs), max(xs)
    y1, y2 = min(ys), max(ys)
    w = x2 - x1
    h = y2 - y1

    if score < 0.55:
        return True

    if w < 8 or h < 8:
        return True

    # Page numbers and seal/logo text are already preserved by the raster layer.
    if re.fullmatch(r"\d{1,4}", text):
        if y1 > IMG_H * 0.84 or x1 > IMG_W * 0.88 or y2 < IMG_H * 0.25:
            return True

    if x1 > IMG_W * 0.86 and y2 < IMG_H * 0.18:
        return True

    if x2 < IMG_W * 0.23 and y2 < IMG_H * 0.25 and not text.startswith("Part"):
        return True

    # The top-left Part ribbon is already clean in the raster layer; OCR overlay
    # often creates visible mask steps because the label sits tightly in boxes.
    if y2 < IMG_H * 0.09 and x1 < IMG_W * 0.28:
        return True

    if re.fullmatch(r"[|I一\-—_ ]+", text):
        return True

    return False


def sample_colors(img_rgb: np.ndarray, box: list[list[float]]) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    xs = [int(round(p[0])) for p in box]
    ys = [int(round(p[1])) for p in box]
    x1 = max(0, min(xs) - 2)
    x2 = min(img_rgb.shape[1], max(xs) + 3)
    y1 = max(0, min(ys) - 2)
    y2 = min(img_rgb.shape[0], max(ys) + 3)
    crop = img_rgb[y1:y2, x1:x2]
    if crop.size == 0:
        return (255, 255, 255), (0, 70, 168)

    pixels = crop.reshape(-1, 3).astype(np.float32)
    if len(pixels) > 6000:
        pixels = pixels[:: max(1, len(pixels) // 6000)]

    k = min(4, max(1, len(pixels) // 8))
    if k <= 1:
        avg = tuple(int(x) for x in np.mean(pixels, axis=0))
        return avg, (0, 70, 168)

    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 20, 1.0)
    _, labels, centers = cv2.kmeans(
        pixels,
        k,
        None,
        criteria,
        3,
        cv2.KMEANS_PP_CENTERS,
    )
    counts = np.bincount(labels.flatten(), minlength=k)
    bg_idx = int(np.argmax(counts))
    bg = centers[bg_idx]

    fg_idx = bg_idx
    best = -1.0
    for i, center in enumerate(centers):
        if i == bg_idx or counts[i] < max(8, len(pixels) * 0.01):
            continue
        dist = color_distance(center, bg)
        score = dist * math.sqrt(float(counts[i]) / max(1.0, float(len(pixels))))
        if score > best:
            best = score
            fg_idx = i

    fg = centers[fg_idx] if fg_idx != bg_idx else np.array([0, 70, 168], dtype=np.float32)

    bg_rgb = tuple(int(round(max(0, min(255, c)))) for c in bg)
    fg_rgb_raw = tuple(int(round(max(0, min(255, c)))) for c in fg)

    # Snap common slide colors so the overlay stays visually clean.
    bg_lum = luminance(bg_rgb)
    fg_lum = luminance(fg_rgb_raw)
    r, g, b = fg_rgb_raw
    if fg_lum > 205 and bg_lum < 150:
        fg_rgb = (255, 255, 255)
    elif b > r + 35 and b >= g:
        fg_rgb = (0, 70, 168)
    elif fg_lum < 95:
        fg_rgb = (20, 20, 20)
    elif r > g + 35 and r > b + 35:
        fg_rgb = (190, 40, 40)
    else:
        fg_rgb = fg_rgb_raw

    if bg_lum > 232:
        bg_rgb = (255, 255, 255)

    return bg_rgb, fg_rgb


def make_item(page: int, img_rgb: np.ndarray, raw_item: list) -> dict | None:
    box, raw_text, raw_score = raw_item
    text = normalize_text(str(raw_text))
    score = float(raw_score)
    if should_skip(page, text, score, box):
        return None

    xs = [float(p[0]) for p in box]
    ys = [float(p[1]) for p in box]
    x1, x2 = min(xs), max(xs)
    y1, y2 = min(ys), max(ys)
    w = x2 - x1
    h = y2 - y1
    bg_rgb, fg_rgb = sample_colors(img_rgb, box)

    h_pt = h * SY
    font_size = max(7.0, min(44.0, h_pt * 0.84))
    if font_size > 24:
        font_size *= 0.94
    if len(text) <= 4 and font_size > 16:
        font_size *= 0.92

    fg_lum = luminance(fg_rgb)
    is_blue = fg_rgb[2] > fg_rgb[0] + 35 and fg_rgb[2] >= fg_rgb[1]
    is_white = fg_lum > 220
    bold = bool(font_size >= 15.0 and (is_blue or is_white or len(text) <= 12))

    return {
        "page": page,
        "text": text,
        "score": round(score, 4),
        "left": round(x1 * SX, 3),
        "top": round(y1 * SY, 3),
        "width": round(w * SX, 3),
        "height": round(h * SY, 3),
        "fontSize": round(font_size, 2),
        "bold": bold,
        "bg": rgb_to_hex(bg_rgb),
        "fg": rgb_to_hex(fg_rgb),
    }


def rect_intersection(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> float:
    ax1, ay1, ax2, ay2 = a
    bx1, by1, bx2, by2 = b
    ix1, iy1 = max(ax1, bx1), max(ay1, by1)
    ix2, iy2 = min(ax2, bx2), min(ay2, by2)
    if ix2 <= ix1 or iy2 <= iy1:
        return 0.0
    return (ix2 - ix1) * (iy2 - iy1)


def make_transparent_asset(crop: np.ndarray, bg_rgb: tuple[int, int, int]) -> Image.Image:
    arr = crop.astype(np.float32)
    bg = np.array((255, 255, 255), dtype=np.float32)
    diff = np.sqrt(np.sum((arr - bg) ** 2, axis=2))
    alpha = np.clip((diff - 14.0) * 255.0 / 34.0, 0, 255).astype(np.uint8)
    rgba = np.dstack([crop.astype(np.uint8), alpha])
    return Image.fromarray(rgba, "RGBA")


def save_asset(page: int, name: str, crop: np.ndarray, left_px: int, top_px: int) -> dict:
    asset_img = make_transparent_asset(crop, (255, 255, 255))
    asset_path = ASSET_DIR / f"slide{page:02d}_{name}.png"
    asset_img.save(asset_path)
    return {
        "page": page,
        "file": str(asset_path),
        "left": round(left_px * SX, 3),
        "top": round(top_px * SY, 3),
        "width": round(crop.shape[1] * SX, 3),
        "height": round(crop.shape[0] * SY, 3),
        "bg": "#FFFFFF",
        "area": round(float(crop.shape[0] * crop.shape[1]), 1),
        "density": 0.0,
        "manual": True,
    }


def save_masked_asset(
    page: int,
    name: str,
    source: np.ndarray,
    mask: np.ndarray,
    left_px: int,
    top_px: int,
    kind: str,
) -> dict | None:
    if int(np.count_nonzero(mask)) < 12:
        return None
    crop = np.full_like(source, 255)
    crop[mask > 0] = source[mask > 0]
    asset_img = make_transparent_asset(crop, (255, 255, 255))
    asset_path = ASSET_DIR / f"slide{page:02d}_{name}.png"
    asset_img.save(asset_path)
    area = float(crop.shape[0] * crop.shape[1])
    return {
        "page": page,
        "file": str(asset_path),
        "left": round(left_px * SX, 3),
        "top": round(top_px * SY, 3),
        "width": round(crop.shape[1] * SX, 3),
        "height": round(crop.shape[0] * SY, 3),
        "bg": "#FFFFFF",
        "area": round(area, 1),
        "density": round(float(np.count_nonzero(mask)) / max(1.0, area), 4),
        "kind": kind,
    }


def split_frame_candidate(
    page: int,
    idx: int,
    cand: dict,
    visual_img: np.ndarray,
    non_white: np.ndarray,
) -> list[dict]:
    x1, y1, x2, y2 = cand["rect"]
    cw, ch = x2 - x1, y2 - y1
    if cw < 120 or ch < 58:
        return []
    if cw > IMG_W * 0.95 or ch > IMG_H * 0.78:
        return []
    if max(cw / max(1, ch), ch / max(1, cw)) > 9.0:
        return []
    if float(cand["density"]) > 0.24:
        return []

    mask = non_white[y1:y2, x1:x2]
    raw_count = int(np.count_nonzero(mask))
    if raw_count < 80:
        return []

    band = max(12, min(28, int(min(cw, ch) * 0.12)))
    top = int(np.count_nonzero(mask[:band, :]))
    bottom = int(np.count_nonzero(mask[max(0, ch - band) :, :]))
    left = int(np.count_nonzero(mask[:, :band]))
    right = int(np.count_nonzero(mask[:, max(0, cw - band) :]))
    min_h = max(35, int(cw * 0.10))
    min_v = max(25, int(ch * 0.10))
    has_frame_band = (top > min_h and bottom > min_h) or (left > min_v and right > min_v)
    if not has_frame_band:
        return []

    num, labels, stats, _ = cv2.connectedComponentsWithStats(mask, 8)
    frame_mask = np.zeros_like(mask)
    detail_mask = np.zeros_like(mask)
    touch_pad = max(12, min(24, int(min(cw, ch) * 0.10)))
    frame_area = 0
    detail_area = 0

    for comp in range(1, num):
        cx, cy, bw, bh, area = [int(v) for v in stats[comp]]
        if area < 12:
            continue
        comp_mask = labels == comp
        near_edge = (
            cx <= touch_pad
            or cy <= touch_pad
            or cx + bw >= cw - touch_pad
            or cy + bh >= ch - touch_pad
        )
        long_h = bw > cw * 0.34 and bh <= max(14, int(ch * 0.16))
        long_v = bh > ch * 0.34 and bw <= max(14, int(cw * 0.13))
        if near_edge or long_h or long_v:
            frame_mask[comp_mask] = 255
            frame_area += area
        else:
            detail_mask[comp_mask] = 255
            detail_area += area

    if frame_area < 40 or detail_area < 25:
        return []
    if frame_area < raw_count * 0.10:
        return []

    out: list[dict] = []
    source = visual_img[y1:y2, x1:x2]
    frame_asset = save_masked_asset(page, f"asset{idx:02d}_frame", source, frame_mask, x1, y1, "frame")
    if frame_asset is not None:
        out.append(frame_asset)

    grouped = cv2.dilate(detail_mask, cv2.getStructuringElement(cv2.MORPH_RECT, (7, 7)), iterations=1)
    gnum, glabels, gstats, _ = cv2.connectedComponentsWithStats(grouped, 8)
    detail_index = 1
    for group in range(1, gnum):
        gx, gy, gw, gh, garea = [int(v) for v in gstats[group]]
        if garea < 20:
            continue
        pad = 6
        dx1 = max(0, gx - pad)
        dy1 = max(0, gy - pad)
        dx2 = min(cw, gx + gw + pad)
        dy2 = min(ch, gy + gh + pad)
        detail_crop_mask = detail_mask[dy1:dy2, dx1:dx2]
        detail_source = source[dy1:dy2, dx1:dx2]
        detail_asset = save_masked_asset(
            page,
            f"asset{idx:02d}_detail{detail_index:02d}",
            detail_source,
            detail_crop_mask,
            x1 + dx1,
            y1 + dy1,
            "detail",
        )
        if detail_asset is not None:
            out.append(detail_asset)
            detail_index += 1

    return out if len(out) > 1 else []


def save_frame_asset(page: int, frame_def: dict, img_rgb: np.ndarray, text_items: list[dict]) -> dict:
    x1, y1, x2, y2 = frame_def["rect"]
    crop = img_rgb[y1:y2, x1:x2].copy()
    frame_rect = (float(x1), float(y1), float(x2), float(y2))

    def erase_rect(rect: tuple[float, float, float, float], pad: int) -> None:
        ex1, ey1, ex2, ey2 = rect
        if rect_intersection(frame_rect, rect) <= 0:
            return
        rx1 = max(0, int(round(ex1 - x1)) - pad)
        ry1 = max(0, int(round(ey1 - y1)) - pad)
        rx2 = min(crop.shape[1], int(round(ex2 - x1)) + pad)
        ry2 = min(crop.shape[0], int(round(ey2 - y1)) + pad)
        if rx2 > rx1 and ry2 > ry1:
            crop[ry1:ry2, rx1:rx2] = 255

    for item in text_items:
        item_rect = (
            float(item["left"]) / SX,
            float(item["top"]) / SY,
            (float(item["left"]) + float(item["width"])) / SX,
            (float(item["top"]) + float(item["height"])) / SY,
        )
        erase_rect(item_rect, 8)

    for rect in frame_def.get("erase", []):
        erase_rect(tuple(float(v) for v in rect), 10)

    return save_asset(page, frame_def["name"], crop, x1, y1)


def extract_assets(page: int, img_rgb: np.ndarray, text_items: list[dict]) -> list[dict]:
    h, w = img_rgb.shape[:2]
    visual_img = img_rgb.copy()

    for item in text_items:
        x1 = max(0, int((float(item["left"]) / SX) - 4))
        y1 = max(0, int((float(item["top"]) / SY) - 3))
        x2 = min(w, int(((float(item["left"]) + float(item["width"])) / SX) + 5))
        y2 = min(h, int(((float(item["top"]) + float(item["height"])) / SY) + 4))
        bg_hex = str(item.get("bg", "#FFFFFF")).lstrip("#")
        bg_rgb = tuple(int(bg_hex[i : i + 2], 16) for i in (0, 2, 4))
        visual_img[y1:y2, x1:x2] = np.array(bg_rgb, dtype=np.uint8)

    diff = np.sqrt(np.sum((visual_img.astype(np.float32) - 255.0) ** 2, axis=2))
    r = visual_img[:, :, 0].astype(np.int16)
    g = visual_img[:, :, 1].astype(np.int16)
    b = visual_img[:, :, 2].astype(np.int16)
    blue_hint = (b > 90) & (b > r + 18) & (b > g - 8)
    non_white = ((diff > 20.0) | blue_hint).astype(np.uint8) * 255

    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (9, 9))
    grouped = cv2.dilate(non_white, kernel, iterations=1)
    num, labels, stats, _ = cv2.connectedComponentsWithStats(grouped, 8)

    candidates: list[dict] = []
    for i in range(1, num):
        x, y, bw, bh, area = [int(v) for v in stats[i]]
        if bw < 10 or bh < 10 or area < 45:
            continue
        if bw > w * 0.98 and bh > h * 0.98:
            continue

        x1 = max(0, x - 8)
        y1 = max(0, y - 8)
        x2 = min(w, x + bw + 8)
        y2 = min(h, y + bh + 8)

        asset_rect = (x1, y1, x2, y2)
        text_overlap = 0.0
        asset_area = float((x2 - x1) * (y2 - y1))
        for item in text_items:
            tr = (
                float(item["left"]) / SX,
                float(item["top"]) / SY,
                (float(item["left"]) + float(item["width"])) / SX,
                (float(item["top"]) + float(item["height"])) / SY,
            )
            text_overlap += rect_intersection(asset_rect, tr)
        if text_overlap / max(1.0, asset_area) > 0.80:
            continue

        crop = img_rgb[y1:y2, x1:x2]
        edge_pixels = np.concatenate(
            [
                crop[: min(4, crop.shape[0]), :, :].reshape(-1, 3),
                crop[max(0, crop.shape[0] - 4) :, :, :].reshape(-1, 3),
                crop[:, : min(4, crop.shape[1]), :].reshape(-1, 3),
                crop[:, max(0, crop.shape[1] - 4) :, :].reshape(-1, 3),
            ],
            axis=0,
        )
        bg = tuple(int(round(v)) for v in np.median(edge_pixels, axis=0))

        raw_pixels = int(np.count_nonzero(non_white[y1:y2, x1:x2]))
        density = raw_pixels / max(1.0, asset_area)
        if density < 0.003:
            continue

        candidates.append(
            {
                "rect": (x1, y1, x2, y2),
                "area": asset_area,
                "bg": bg,
                "density": density,
                "crop": crop,
            }
        )

    candidates.sort(key=lambda c: c["area"], reverse=True)
    kept: list[dict] = []
    for cand in candidates:
        rect = cand["rect"]
        duplicate = False
        for prev in kept:
            inter = rect_intersection(rect, prev["rect"])
            if inter / min(cand["area"], prev["area"]) > 0.72:
                duplicate = True
                break
        if not duplicate:
            kept.append(cand)

    ASSET_DIR.mkdir(exist_ok=True)
    for old in ASSET_DIR.glob(f"slide{page:02d}_asset*.png"):
        old.unlink()
    for old in ASSET_DIR.glob(f"slide{page:02d}_card*.png"):
        old.unlink()

    assets: list[dict] = []
    for idx, cand in enumerate(kept, 1):
        x1, y1, x2, y2 = cand["rect"]
        split_assets = split_frame_candidate(page, idx, cand, visual_img, non_white)
        if split_assets:
            assets.extend(split_assets)
            continue
        crop = visual_img[y1:y2, x1:x2]
        asset_img = make_transparent_asset(crop, cand["bg"])
        asset_path = ASSET_DIR / f"slide{page:02d}_asset{idx:02d}.png"
        asset_img.save(asset_path)
        assets.append(
            {
                "page": page,
                "file": str(asset_path),
                "left": round(x1 * SX, 3),
                "top": round(y1 * SY, 3),
                "width": round((x2 - x1) * SX, 3),
                "height": round((y2 - y1) * SY, 3),
                "bg": rgb_to_hex(cand["bg"]),
                "area": round(cand["area"], 1),
                "density": round(cand["density"], 4),
            }
        )
    if page == 2:
        assets = [
            asset
            for asset in assets
            if not (
                float(asset["left"]) > 170.0
                and float(asset["width"]) > 300.0
                and float(asset["height"]) > 140.0
                and float(asset["top"]) < 460.0
            )
        ]
    for frame_def in MANUAL_FRAME_CROPS.get(page, []):
        assets.append(save_frame_asset(page, frame_def, img_rgb, text_items))
    for crop_def in MANUAL_CROPS.get(page, []):
        x1, y1, x2, y2 = crop_def["rect"]
        crop = img_rgb[y1:y2, x1:x2]
        assets.append(save_asset(page, crop_def["name"], crop, x1, y1))
    return assets


def main() -> None:
    ocr = RapidOCR()
    ASSET_DIR.mkdir(exist_ok=True)
    for old in ASSET_DIR.glob("slide*_asset*.png"):
        old.unlink()
    slides: list[dict] = []
    for page in range(1, 22):
        img_path = ROOT / f"{page}.png"
        img = Image.open(img_path).convert("RGB")
        img_rgb = np.array(img)
        result, elapsed = ocr(str(img_path))
        items = []
        for raw_item in result or []:
            item = make_item(page, img_rgb, raw_item)
            if item is not None:
                items.append(item)
        items.sort(key=lambda x: (x["top"], x["left"]))
        mode = "full_background" if page in FULL_BACKGROUND_PAGES else "white_asset"
        assets = [] if mode == "full_background" else extract_assets(page, img_rgb, items)
        slides.append(
            {
                "page": page,
                "mode": mode,
                "image": str(img_path),
                "elapsed": elapsed,
                "itemCount": len(items),
                "assetCount": len(assets),
                "items": items,
                "assets": assets,
                "cleanups": MANUAL_CLEANUPS.get(page, []),
                "nativeShapes": NATIVE_SHAPES.get(page, []),
                "editableTexts": sorted(FULL_BACKGROUND_EDITABLE_TEXT.get(page, [])),
            }
        )
        print(
            f"slide {page:02d}: {mode}, {len(items)} editable text items, "
            f"{len(assets)} movable assets"
        )

    payload = {
        "slideWidth": SLIDE_W_PT,
        "slideHeight": SLIDE_H_PT,
        "sourceImageWidth": IMG_W,
        "sourceImageHeight": IMG_H,
        "slides": slides,
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
