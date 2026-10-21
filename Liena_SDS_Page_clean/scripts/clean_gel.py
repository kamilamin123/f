"""Remove Image Lab lane/band annotations from a report gel image.

Overlays drawn by Image Lab on the report image:
  * blue lane frames, lane boundaries, lane numbers and "Std" label
  * magenta band lines
  * black molecular-weight labels and ticks beside the lane frame
  * red (white-edged) highlight on saturated pixels

Two report colour schemes are handled:
  * grayscale (Coomassie / white light): gray pixels are gel, coloured
    pixels are overlay.
  * "Stain Free" colormap (cyan ramp): luma is proportional to intensity
    and chroma follows the ramp, so pixels off the ramp are overlay.
    Saturated (red/white) pixels are set to full white, as in the raw scan.
Masked pixels are filled by biharmonic inpainting from the surrounding gel.
The output is grayscale: as displayed for grayscale reports, white bands on
black for stain-free reports.
"""
import sys

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.restoration import inpaint_biharmonic

# Stain Free colormap: Cb ~ CB_RAMP*Y, Cr ~ CR_RAMP*Y, and Y = STAINFREE_MAX_Y
# at full scale (fitted against the raw .scn scan of the 14 Oct gel)
CB_RAMP, CR_RAMP, STAINFREE_MAX_Y = 0.533, -0.519, 146.0
# Chroma above which a pixel's luma is taken as overlay. Tuned against the raw
# scans of 23 Sep (grayscale) and 14 Oct (stain-free): lower values mask gel
# luma under the JPEG chroma fringe (stained bands lose contrast), higher
# values leave faint line edges.
MASK_THR = 30


def load_ycbcr(path):
    im = Image.open(path)
    im.draft("YCbCr", im.size)  # decoded YCbCr, no RGB round trip
    ycc = np.asarray(im.convert("YCbCr")).astype(float)
    return ycc[..., 0], ycc[..., 1] - 128, ycc[..., 2] - 128


def is_stain_free(y, cb, cr):
    # the cyan ramp has strongly negative Cr on anything above black
    lit = y > 10
    return np.median(cr[lit] / y[lit]) < CR_RAMP / 2


def frame_lines(blue):
    """Outer lane frame: column spans (incl. JPEG bleed) of its left and right
    vertical lines, and its top and bottom rows."""
    lab, n = ndi.label(ndi.binary_dilation(blue, iterations=2))
    sizes = ndi.sum(np.ones_like(lab), lab, range(1, n + 1))
    frame = (lab == (np.argmax(sizes) + 1)) & blue
    rows = np.where(frame.any(1))[0]
    y0, y1 = rows.min(), rows.max()
    xs = np.where(frame.sum(0) > 0.5 * (y1 - y0))[0]
    breaks = np.where(np.diff(xs) > 1)[0]
    left = (xs[0], xs[breaks[0]] if len(breaks) else xs[-1])
    right = (xs[breaks[-1] + 1] if len(breaks) else xs[0], xs[-1])
    return left, right, y0, y1


def label_mask(y, left, right, y0, y1, stain_free):
    """Black MW labels/ticks in the columns left and right of the frame."""
    bg = ndi.median_filter(y, size=31)
    if stain_free:  # black text on a near-black background
        dark = (y < 2.5) & (y < bg - 2.5)
    else:
        dark = (y < bg - 45) & (y < 0.6 * bg)
    # ticks touch the frame line, so reach a few px past it
    side = np.zeros_like(dark)
    side[max(y0 - 15, 0):y1 + 15, max(left[0] - 60, 0):left[1] + 3] = True
    side[max(y0 - 15, 0):y1 + 15, right[0] - 2:right[1] + 60] = True
    # the blue frame line is dark too; drop it so ticks stay separate blobs
    side[:, left[0] - 1:left[1] + 2] = side[:, right[0] - 1:right[1] + 2] = False
    text = dark & side
    lab, n = ndi.label(text, structure=np.ones((3, 3)))
    mask = np.zeros_like(text)
    pad = 3  # covers anti-aliasing and JPEG ringing around the glyphs
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        big_enough = (lab[sl] == i).sum() >= 4  # not JPEG noise
        if h <= 30 and w <= 60 and big_enough:  # overlapping labels merge
            xa, xb = max(sl[1].start - pad, 0), sl[1].stop + pad
            # ticks run onto the frame line: carry them across it
            if xb >= left[0] - 1 and xa < left[0]:
                xb = left[1] + 3
            if xa <= right[1] + 1 and xb > right[1]:
                xa = right[0] - 3
            mask[max(sl[0].start - pad, 0):sl[0].stop + pad, xa:xb] = True
    return mask


def analyse(path, mask_thr=MASK_THR):
    """Gel value in [0, 1], overlay mask and lane frame of a report image.

    Band lines (magenta) cross stained bands, so only pixels whose chroma
    exceeds mask_thr are masked: JPEG spreads the overlay's colour (4:2:0
    chroma) further than its luma, and the luma of those fringe pixels is
    still the gel's. On the stain-free ramp, thin blue lane/frame lines keep
    too little excess chroma for that, so there blue is masked down to the
    detection threshold.
    """
    y, cb, cr = load_ycbcr(path)
    stain_free = is_stain_free(y, cb, cr)
    if stain_free:
        dcb, dcr = cb - CB_RAMP * y, cr - CR_RAMP * y
        value = np.clip(y / STAINFREE_MAX_Y, 0, 1)
        saturated = ((cr > 50) & (cb < 0)) | (y > 165)  # red fill, white edge
        value[saturated] = 1.0
        detect_thr = 20
    else:
        dcb, dcr = cb, cr
        value = y / 255.0
        saturated = np.zeros_like(y, bool)
        detect_thr = 14
    chroma = np.hypot(dcb, dcr)

    # blue overlay: excess Cb, little excess Cr (on the stain-free ramp a blue
    # line's own luma shows as up to ~+8); magenta/red have excess Cr > 50
    blue = (chroma > detect_thr) & (dcb > 12) & (dcr < 12)
    left, right, y0, y1 = frame_lines(blue)

    mask = chroma > mask_thr
    if stain_free:
        mask |= blue
    mask |= label_mask(y, left, right, y0, y1, stain_free)
    mask &= ~saturated  # saturated pixels are known: full intensity
    return value, mask, (left[0], y0, right[1], y1)


def clean(path, **kw):
    value, mask, frame = analyse(path, **kw)
    out = inpaint_biharmonic(value, mask)
    return (np.clip(out, 0, 1) * 255).round().astype(np.uint8), mask, frame


def clean_merged(paths, **kw):
    """Clean several reports of the same scan (different lane/band analyses):
    a pixel is only inpainted where every report's overlay covers it."""
    values, masks, frames = zip(*(analyse(p, **kw) for p in paths))
    known = ~np.array(masks)
    count = known.sum(0)
    value = (np.array(values) * known).sum(0) / np.maximum(count, 1)
    mask = count == 0
    out = inpaint_biharmonic(value, mask)
    return (np.clip(out, 0, 1) * 255).round().astype(np.uint8), mask, frames[0]


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    img, mask, bbox = clean(src)
    Image.fromarray(img).save(dst)
    if len(sys.argv) > 3:
        Image.fromarray((mask * 255).astype(np.uint8)).save(sys.argv[3])
    print(dst, "masked px:", int(mask.sum()), "frame:", bbox)
