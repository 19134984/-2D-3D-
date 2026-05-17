---
name: image-to-ppt-lookalike
description: Rebuild one or more flat images, screenshots, scanned slides, paper figures, diagrams, or UI captures into PowerPoint slides that visually match the source image while making the practical parts editable. Use when the user asks for image-to-PPT conversion, editable PPT reconstruction, screenshot-to-PowerPoint, picture-to-slide, or a VBA/Python-generated PPT where "looks like the original" matters more than perfect object recovery.
---

# Image To Ppt Lookalike

## Overview

Reconstruct flat visual sources as editable PowerPoint decks using a hybrid strategy: editable native objects for text and simple geometry, cropped raster images for complex visual regions. Optimize first for visual resemblance to the source, then for reasonable editability.

## Core Principle

Treat the source image as the visual contract. Do not promise full layer recovery from a flat bitmap. Preserve complex regions as image crops when rebuilding them as shapes would reduce fidelity, take too long, or introduce errors.

Use this editability hierarchy:

1. Recreate text as PowerPoint text boxes when OCR is reliable.
2. Recreate simple rectangles, lines, arrows, circles, tables, and frames as native shapes.
3. Recreate simple charts only when data or clean structure is available.
4. Insert cropped image regions for photos, heatmaps, screenshots, dense plots, formulas, logos, icons, and highly detailed artwork.
5. Keep a hidden or off-slide copy of the original image only when useful for later manual comparison.

## Hybrid Background Strategy

Choose a per-slide reconstruction mode before rebuilding:

- **Full-background mode**: For slides with complex gradients, curved bands, dense textures, simulation fields, or large atmospheric visuals, keep the full-slide image as the visual base. Add only high-confidence editable overlays, or leave the slide as a visual background when overlays would visibly damage the design.
- **White-base asset mode**: For slides that are mostly white or flat solid background, replace the full-slide raster with a white/solid base. Extract logos, icons, decorative dots, frame lines, arrows, page-number blocks, and other non-text visuals as separate cropped transparent images or native shapes, then place editable OCR text above them.
- When moving an extracted asset would otherwise reveal the original underneath, first cover the original region with the slide background color, then reinsert the crop as its own object.
- Keep extraction pragmatic: if a crop has messy edges, loses shadows, or breaks a complex composition, revert that region or slide to full-background mode. Mixed decks are acceptable and often best.
- For white-base asset mode, avoid broad text masks after assets are placed; use no mask or very tight cleanup patches so card borders and icon lines are not accidentally covered.

## Workflow

1. Inspect the source image dimensions, aspect ratio, and visual layout. Set the PPT slide size to match the source aspect ratio before placing elements.
2. Place the source image as a full-slide temporary background layer for alignment.
3. Identify visual layers: background, major panels, text blocks, borders, arrows, simple icons, charts, formulas, and complex raster regions.
4. Choose full-background mode or white-base asset mode for each slide. Use a mixed approach across a deck when some pages are complex and others are mostly white.
5. OCR text blocks when possible. Use editable text boxes for reliable OCR, but keep text as image crops when OCR is poor, formula-heavy, or visually fragile.
6. Rebuild simple geometry with native PowerPoint shapes. Match fill, line color, line width, transparency, z-order, and approximate corner radius.
7. Crop and insert complex image regions and movable visual assets. Preserve them at source resolution whenever possible; avoid resampling small text-heavy crops.
8. Align all reconstructed elements over the temporary background, then hide/remove the background only on slides where the rebuilt base and assets fully replace it.
9. Render or screenshot the resulting slide and compare it against the source. Iterate on positions, crop boundaries, font sizes, and z-order until the slide looks close.

## Implementation Choices

Prefer the simplest generation path that works in the current environment:

- Use PowerPoint automation or VBA when the user explicitly wants VBA code, Office-native generation, or manual handoff inside PowerPoint.
- Use Python libraries such as `python-pptx`, image processing tools, and OCR when batch generation, repeatability, or non-interactive execution matters.
- Use a layout JSON intermediate when the reconstruction has many elements: store slide size, element type, coordinates, text, crop path, style, and z-order, then generate PPT from that structured layout.

Coordinates should be normalized or consistently scaled from image pixels to slide units. Keep the original image size and slide size recorded so later edits can preserve alignment.

## Reconstruction Rules

- Visual match beats theoretical editability. If a region is visually complex, insert a crop.
- Keep text editable only when the result still looks like the source. Bad OCR is worse than a clean crop.
- Rebuild frames and layout scaffolding as native shapes when simple geometry matches well. When the source has subtle corner radius, anti-aliased borders, clipped icons, or the user wants independent composition, crop the frame/border/long divider lines as a transparent bottom image layer and place separated icons, numbers, and editable OCR text above it.
- On mostly white slides, prefer separating small logos, icons, borders, arrows, and page markers into movable image objects or native shapes instead of trapping them in a full-slide background.
- Use real PowerPoint tables only for clean, grid-like tables. For dense or scanned tables, prefer a crop plus optional editable overlay labels.
- Preserve formulas as crops unless the user asks for editable equations and accepts manual correction.
- Do not use decorative reinterpretation. Match the source colors, spacing, and proportions.
- Do not create a marketing-style redesign unless the user asks for redesign. This skill is for faithful reconstruction.

## Quality Check

Before final delivery, verify:

- Slide aspect ratio matches the source.
- Main blocks align within a few pixels or a visibly acceptable tolerance.
- Text does not overflow boxes and does not shift neighboring content.
- Cropped image regions are sharp enough at the intended slide size.
- The final slide no longer shows the temporary full-slide background unless intentionally included.
- The user can edit the meaningful scaffolding: titles, ordinary labels, frames, arrows, and major boxes.

For detailed acceptance criteria, read `references/visual-matching-checklist.md`.

## User-Facing Caveat

When explaining results, describe the deck as "visually reconstructed and partially editable" unless every important element has truly been rebuilt as a native PowerPoint object.
