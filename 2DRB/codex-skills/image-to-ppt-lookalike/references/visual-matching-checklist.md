# Visual Matching Checklist

Use this checklist when reconstructing an image as a PowerPoint slide.

## Fidelity Priorities

1. Overall composition and aspect ratio.
2. Major panel positions and sizes.
3. Text placement, approximate font weight, and hierarchy.
4. Border, arrow, and line alignment.
5. Color similarity.
6. Image crop sharpness.
7. Editability of the parts the user is likely to change.

## What To Rebuild

- Rebuild as editable text: clear headings, labels, short paragraphs, legends, and callouts.
- Rebuild as editable shapes: boxes, frames, simple arrows, separators, simple table grids, flowchart nodes.
- Keep as image crops: photos, rendered simulations, heatmaps, dense plots, formulas, screenshots, logos, detailed icons, handwritten or scanned text.

## Comparison Method

If possible, export or screenshot the generated slide at the same pixel size as the source image. Compare the output visually against the source, focusing on visible drift rather than exact internal object structure.

Check for:

- Background removal mistakes.
- Crops that include unwanted neighboring pixels.
- Text boxes that wrap differently from the source.
- Misaligned borders after deleting the temporary background.
- Font substitution that changes line breaks.
- Objects placed in the wrong z-order.

## Delivery Language

Say clearly which parts are editable native PowerPoint objects and which parts are preserved as image crops. Avoid claiming that the original editable layers were recovered from a flat image.
