# Liena SDS-PAGE reports: gel images without lane/band annotations

All 14 Image Lab reports (22–23 September and 14–15 October 2026), with the
page-1 gel image replaced by a clean version: no lane frames, lane numbers,
band lines, MW labels or red saturation highlighting, in grayscale like
Figure S6 of Kim et al. (Biosens. Bioelectron., supplementary information).
All other pages, tables and lane-profile plots are unchanged.

- `*.pdf`: the cleaned reports (same file names as the originals)
- `gel_images/*.png`: the cleaned gel images at full resolution (1392 x 1000)

## How the images were cleaned

- **23 Sep (Coomassie) and 14 Oct (stain-free)**: the original, unannotated
  scan was read from the matching `.scn` file and rendered with the same
  display range and polarity as the report (dark bands on white for 23 Sep,
  white bands on black for 14 Oct). No pixels are reconstructed.
- **22 Sep and 15 Oct**: no `.scn` files were available, so the overlays
  were removed from the report image itself and the masked pixels filled by
  biharmonic inpainting from the surrounding gel.
  - 15 Oct (Coomassie): the gel image is grayscale, so coloured pixels (blue
    lanes/numbers, magenta band lines) and the black MW labels beside the
    lane frame are masked.
  - 22 Sep (stain-free): the report uses Image Lab's cyan "Stain Free"
    colormap, which was converted back to grayscale (luma is proportional
    to intensity; pixels off the colormap are overlay). Saturated pixels
    (red with a white edge) are shown as full white, as in the raw scan.
    Checked against the 14 Oct raw scan, this reconstruction is within
    1.6 gray levels on average (3.9 on the inpainted pixels).

## Reproducing

```
pip install numpy scipy scikit-image pillow pypdf
python scripts/scn_image.py "<report>.scn" out.png     # 23 Sep, 14 Oct
python scripts/clean_gel.py page1_image.jpg out.png    # 22 Sep, 15 Oct (pdfimages -f 1 -l 1 -j)
python scripts/replace_gel.py "<report>.pdf" out.png "<cleaned>.pdf"
```
