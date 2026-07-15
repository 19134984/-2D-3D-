from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageStat


def page_number(path: Path) -> int:
    match = re.search(r"(\d+)$", path.stem)
    if not match:
        raise ValueError(f"No page number in {path.name}")
    return int(match.group(1))


def page_metrics(path: Path) -> dict[str, object]:
    with Image.open(path) as source:
        image = source.convert("L")
        width, height = image.size
        white = Image.new("L", image.size, 255)
        difference = ImageChops.difference(image, white)
        bbox = difference.point(lambda value: 255 if value > 10 else 0).getbbox()
        histogram = image.histogram()
        dark_pixels = sum(histogram[:245])
        ink_ratio = dark_pixels / (width * height)
        mean = ImageStat.Stat(image).mean[0]
        if bbox is None:
            margins = (width, height, width, height)
        else:
            left, top, right, bottom = bbox
            margins = (left, top, width - right, height - bottom)
        border = 3
        border_crops = (
            image.crop((0, 0, width, border)),
            image.crop((0, height - border, width, height)),
            image.crop((0, 0, border, height)),
            image.crop((width - border, 0, width, height)),
        )
        border_dark = 0
        border_total = 0
        for crop in border_crops:
            hist = crop.histogram()
            border_dark += sum(hist[:245])
            border_total += crop.width * crop.height
        return {
            "page": page_number(path),
            "file": path.name,
            "width": width,
            "height": height,
            "mean_gray": round(mean, 3),
            "ink_ratio": round(ink_ratio, 6),
            "margin_left": margins[0],
            "margin_top": margins[1],
            "margin_right": margins[2],
            "margin_bottom": margins[3],
            "border_dark_ratio": round(border_dark / border_total, 8),
        }


def make_contact_sheet(
    page_paths: list[Path],
    output_path: Path,
    columns: int = 4,
    thumb_width: int = 330,
) -> None:
    rows = math.ceil(len(page_paths) / columns)
    label_height = 30
    gap = 16
    with Image.open(page_paths[0]) as sample:
        aspect = sample.height / sample.width
    thumb_height = round(thumb_width * aspect)
    sheet_width = gap + columns * (thumb_width + gap)
    sheet_height = gap + rows * (label_height + thumb_height + gap)
    sheet = Image.new("RGB", (sheet_width, sheet_height), "#d8d8d8")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default()

    for index, path in enumerate(page_paths):
        row, column = divmod(index, columns)
        x = gap + column * (thumb_width + gap)
        y = gap + row * (label_height + thumb_height + gap)
        page = page_number(path)
        draw.rectangle((x, y, x + thumb_width, y + label_height), fill="white")
        draw.text((x + 8, y + 8), f"Page {page}", fill="black", font=font)
        with Image.open(path) as source:
            thumb = source.convert("RGB")
            thumb.thumbnail((thumb_width, thumb_height), Image.Resampling.LANCZOS)
            canvas = Image.new("RGB", (thumb_width, thumb_height), "white")
            offset = ((thumb_width - thumb.width) // 2, (thumb_height - thumb.height) // 2)
            canvas.paste(thumb, offset)
            sheet.paste(canvas, (x, y + label_height))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(output_path, quality=92)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pages_dir", type=Path)
    parser.add_argument("contacts_dir", type=Path)
    parser.add_argument("metrics_csv", type=Path)
    parser.add_argument("--group-size", type=int, default=12)
    args = parser.parse_args()

    pages = sorted(args.pages_dir.glob("page-*.png"), key=page_number)
    if not pages:
        raise SystemExit(f"No rendered pages in {args.pages_dir}")

    metrics = [page_metrics(path) for path in pages]
    args.metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.metrics_csv.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)

    for start in range(0, len(pages), args.group_size):
        group = pages[start : start + args.group_size]
        first = page_number(group[0])
        last = page_number(group[-1])
        output = args.contacts_dir / f"contact-{first:03d}-{last:03d}.png"
        make_contact_sheet(group, output)

    sizes = {(row["width"], row["height"]) for row in metrics}
    blanks = [row["page"] for row in metrics if float(row["ink_ratio"]) < 0.001]
    borders = [row["page"] for row in metrics if float(row["border_dark_ratio"]) > 0.0001]
    tight = [
        row["page"]
        for row in metrics
        if min(
            int(row["margin_left"]),
            int(row["margin_top"]),
            int(row["margin_right"]),
            int(row["margin_bottom"]),
        )
        < 8
    ]
    print(f"pages={len(pages)}")
    print(f"sizes={sorted(sizes)}")
    print(f"blank_candidates={blanks}")
    print(f"border_candidates={borders}")
    print(f"tight_margin_candidates={tight}")


if __name__ == "__main__":
    main()
