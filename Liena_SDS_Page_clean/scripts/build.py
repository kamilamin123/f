"""Build the cleaned reports, gel images and publication TIFFs.

usage: python build.py <folder with the Image Lab .pdf/.scn files> <output folder>

Per report:      <out>/<report>.pdf, <out>/gel_images/<report>.png
Per scan:        <out>/publication/<date>_<time>_{full,lanes}_600dpi.tif
                 (+ _full_16bit.tif when the raw .scn scan is available)

Reports re-analysing the same scan ("...02min", "...02min2", ...) share one
publication image, inpainted only where every report's overlay covers it.
"""
import io
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image
from pypdf import PdfReader
from skimage.transform import rescale

import clean_gel
import replace_gel
import scn_image

UPSCALE = 3      # bicubic; 1392 x 1000 camera px -> 4176 x 3000 (7.0 x 5.0 in at 600 dpi)
DPI = 600
MARGIN = 15      # camera px kept around the Image Lab lane frame in the lanes crop


def page1_jpeg(pdf):
    """The report's gel image as stored (raw DCT stream, not re-encoded)."""
    xobjects = PdfReader(pdf).pages[0]["/Resources"]["/XObject"]
    gel, = (x.get_object() for x in xobjects.values()
            if x.get_object()["/Subtype"] == "/Image")
    assert gel["/Filter"] == "/DCTDecode"
    return io.BytesIO(gel._data)


def scan_name(stem):
    """Report name -> scan name: drop the copy number after "min"."""
    return re.sub(r"min\d+$", "min", stem)


def short_name(scan):
    """'ChemEng UGlab. User 2026-10-15 15hr 02min' -> '2026-10-15_15h02'."""
    d, h, m = re.search(r"(\d{4}-\d\d-\d\d) (\d+)hr (\d+)min", scan).groups()
    return f"{d}_{h}h{m}"


def save_tiff(v, path, bits=8, dpi=DPI):
    v = np.clip(v, 0, 1)
    if bits == 16:
        im = Image.fromarray((v * 65535).round().astype(np.uint16))
    else:
        im = Image.fromarray((v * 255).round().astype(np.uint8))
    im.save(path, compression="tiff_lzw", dpi=(dpi, dpi))


def build_report(pdf, scn, out):
    """Cleaned gel image and PDF for one report; returns the lane frame."""
    if scn.exists():
        _, _, frame = clean_gel.analyse(page1_jpeg(pdf))
        img = scn_image.render(scn)
    else:
        img, _, frame = clean_gel.clean(page1_jpeg(pdf))
    png = out / "gel_images" / f"{pdf.stem}.png"
    Image.fromarray(img).save(png)
    replace_gel.replace(pdf, png, out / pdf.name)
    return frame


def build_publication(scan, pdfs, scn, frame, out):
    if scn.exists():
        v = scn_image.display(scn)
        save_tiff(v, out / f"{short_name(scan)}_full_16bit.tif", bits=16, dpi=DPI // UPSCALE)
    else:
        img, _, _ = clean_gel.clean_merged([page1_jpeg(p) for p in pdfs])
        v = img / 255.0
    x0, y0, x1, y1 = frame
    lanes = v[max(y0 - MARGIN, 0):y1 + MARGIN + 1, max(x0 - MARGIN, 0):x1 + MARGIN + 1]
    for part, a in (("full", v), ("lanes", lanes)):
        big = rescale(a, UPSCALE, order=3, mode="edge", anti_aliasing=False)
        save_tiff(big, out / f"{short_name(scan)}_{part}_600dpi.tif")


def main(src, out):
    src, out = Path(src), Path(out)
    (out / "gel_images").mkdir(parents=True, exist_ok=True)
    (out / "publication").mkdir(exist_ok=True)
    scans = defaultdict(list)
    for pdf in sorted(src.glob("*.pdf")):
        scans[scan_name(pdf.stem)].append(pdf)

    with ProcessPoolExecutor() as pool:
        frames = {p: pool.submit(build_report, p, src / f"{s}.scn", out)
                  for s, ps in scans.items() for p in ps}
        frames = {p: f.result() for p, f in frames.items()}
        jobs = [pool.submit(build_publication, s, ps, src / f"{s}.scn",
                            frames[ps[0]], out / "publication")
                for s, ps in scans.items()]
        for job in jobs:
            job.result()


if __name__ == "__main__":
    main(*sys.argv[1:3])
