# Liena SDS-PAGE reports: gel images without lane/band annotations

All 14 Image Lab reports (22–23 September and 14–15 October 2026), with the
page-1 gel image replaced by a clean version: no lane frames, lane numbers,
band lines, MW labels or red saturation highlighting, in grayscale like
Figure S6 of Kim et al. (Biosens. Bioelectron., supplementary information).
All other pages, tables and lane-profile plots are unchanged.

- `*.pdf`: the cleaned reports (same file names as the originals)
- `gel_images/*.png`: the cleaned gel images, one per report (1392 x 1000)
- `publication/`: one set per scan (the 15:02 and 14:45 reports are re-analyses
  of the same scan), LZW-compressed grayscale TIFF:
  - `*_lanes_600dpi.tif`: cropped to the analysed lanes (Image Lab lane frame
    plus 15 camera px), for a figure panel
  - `*_full_600dpi.tif`: the whole uncropped image (7.0 x 5.0 in at 600 dpi),
    e.g. for supplementary/source-data figures
  - `*_full_16bit.tif`: 16-bit master at camera resolution, for the two scans
    with raw data (23 Sep, 14 Oct), for further contrast adjustment

The camera resolution is 1392 x 1000 px (107.8 µm/px). The 600 dpi files are
3x bicubic enlargements of that, so they meet journal resolution requirements
but carry no detail beyond it.

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

To keep stained bands intact, only the pixels the band lines actually drew
over are reconstructed: JPEG smears the lines' colour a few pixels further
than their brightness, and those fringe pixels keep the real gel values.
This was tuned against the two raw scans, by cleaning their report images
and comparing with the raw data (gray levels out of 255):

| stained-band pixels next to removed lines | mean error | 95% within |
|---|---|---|
| grayscale (23 Sep) | 1.2 | 3.9 |
| stain-free (14 Oct) | 4.4 | 12 |

For stain-free the remaining error is the colormap-to-grayscale conversion
itself: bands away from any line show the same error.

For journal submission, images re-exported from the original scans in
Image Lab (annotations hidden) are still the safest source, since they
involve no reconstruction at all.

## Reproducing

```
pip install numpy scipy scikit-image pillow pypdf
python scripts/build.py <folder with the Image Lab .pdf/.scn files> <output folder>
```
