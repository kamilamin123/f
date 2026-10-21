"""Remove Image Lab lane/band annotations from a grayscale gel image.

Overlays drawn by Image Lab on the report image:
  * blue lane frames, lane boundaries, lane numbers and "Std" label
  * magenta band lines
  * black molecular-weight labels and ticks beside the lane frame

The underlying gel image is grayscale, so coloured pixels are overlay.
Masked pixels are filled by biharmonic inpainting on the luma channel.
"""
import sys

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.restoration import inpaint_biharmonic


def load_ycbcr(path):
    im = Image.open(path)
    im.draft("YCbCr", im.size)  # decoded YCbCr, no RGB round trip
    ycc = np.asarray(im.convert("YCbCr")).astype(float)
    return ycc[..., 0], ycc[..., 1] - 128, ycc[..., 2] - 128


def overlay_masks(y, cb, cr, chroma_thr=14):
    chroma = np.hypot(cb, cr)
    colour = chroma > chroma_thr
    # blue overlay: Cb strongly positive, Cr negative; magenta: both positive
    blue = colour & (cb > 12) & (cr < 4)

    # lane frame = biggest blue connected component
    lab, n = ndi.label(ndi.binary_dilation(blue, iterations=2))
    sizes = ndi.sum(np.ones_like(lab), lab, range(1, n + 1))
    frame = (lab == (np.argmax(sizes) + 1)) & blue
    rows, cols = np.where(frame)
    y0, y1, x0, x1 = rows.min(), rows.max(), cols.min(), cols.max()

    # black MW labels/ticks in the columns left and right of the frame;
    # ticks touch the frame line, so reach a few px past it
    bg = ndi.median_filter(y, size=31)
    dark = (y < bg - 45) & (y < 0.6 * bg)
    side = np.zeros_like(dark)
    side[max(y0 - 15, 0):y1 + 15, max(x0 - 60, 0):x0 + 3] = True
    side[max(y0 - 15, 0):y1 + 15, x1 - 2:x1 + 60] = True
    colour = ndi.binary_dilation(colour, iterations=1)
    # the blue frame line is dark too; drop it so ticks stay separate blobs.
    # Only the strong blue core: JPEG puts weak chroma noise on black glyphs.
    blue_core = ndi.binary_dilation((cb > 35) & (cr < 5), iterations=1)
    text = dark & side & ~blue_core
    lab, n = ndi.label(text, structure=np.ones((3, 3)))
    text = np.zeros_like(text)
    pad = 3  # covers anti-aliasing and JPEG ringing around the glyphs
    for sl in ndi.find_objects(lab):
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if h <= 30 and w <= 60:  # overlapping labels merge into one blob
            text[max(sl[0].start - pad, 0):sl[0].stop + pad,
                 max(sl[1].start - pad, 0):sl[1].stop + pad] = True
    return colour | text, (x0, y0, x1, y1)


def clean(path):
    y, cb, cr = load_ycbcr(path)
    mask, bbox = overlay_masks(y, cb, cr)
    out = inpaint_biharmonic(y / 255.0, mask) * 255.0
    return np.clip(out, 0, 255).round().astype(np.uint8), mask, bbox


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    img, mask, bbox = clean(src)
    Image.fromarray(img).save(dst)
    if len(sys.argv) > 3:
        Image.fromarray((mask * 255).astype(np.uint8)).save(sys.argv[3])
    print(dst, "masked px:", int(mask.sum()), "frame:", bbox)
