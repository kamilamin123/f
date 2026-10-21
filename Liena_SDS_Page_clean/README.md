# Liena SDS-PAGE reports: gel images without lane/band annotations

The Image Lab reports from 14 and 15 October 2026, with the page-1 gel image
replaced by a clean version: no lane frames, lane numbers, band lines,
MW labels or red saturation highlighting, in grayscale like Figure S6 of
Kim et al. (Biosens. Bioelectron., supplementary information). All other
pages, tables and lane-profile plots are unchanged.

- `*.pdf`: the cleaned reports (same file names as the originals)
- `gel_images/*.png`: the cleaned gel images at full resolution (1392 x 1000)

## How the images were cleaned

- **14 Oct (stain-free gel)**: the original, unannotated scan was read from
  `ChemEng UGlab. User 2026-10-14 14hr 54min.scn` and rendered with the same
  display range as the report (white bands on black). No pixels are
  reconstructed.
- **15 Oct (Coomassie gels)**: no `.scn` files were available, so the
  overlays were removed from the report image itself. The gel image is
  grayscale, so coloured pixels (blue lanes/numbers, magenta band lines) and
  the black MW labels beside the lane frame were masked and filled by
  biharmonic inpainting from the surrounding gel.

## Reproducing

```
pip install numpy scipy scikit-image pillow pypdf
python scripts/scn_image.py "<report>.scn" out.png     # 14 Oct
python scripts/clean_gel.py page1_image.jpg out.png    # 15 Oct (pdfimages -f 1 -l 1 -j)
python scripts/replace_gel.py "<report>.pdf" out.png "<cleaned>.pdf"
```
