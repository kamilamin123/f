"""Render the raw (annotation-free) image stored in a Bio-Rad Image Lab .scn file.

Uses the same display range and polarity as the report (transform low/high
fractions of the scanner max value), in grayscale with no saturation
highlight: white bands on black for the Stain Free colormap, dark bands on
white for grayscale (e.g. Coomassie) scans.
"""
import re
import sys

import numpy as np
from PIL import Image


def part(blob, description):
    start = blob.index(b"Content-Description: " + description)
    header_end = blob.index(b"\r\n\r\n", start) + 4
    length = int(re.search(rb"Content-Length: (\d+)",
                           blob[blob.rfind(b"Content-Type", 0, start):start]).group(1))
    return blob[header_end:header_end + length]


def render(path):
    blob = open(path, "rb").read()
    header = part(blob, b"ImageHeader").decode("utf8")
    width = int(re.search(r'<size_pix width="(\d+)"', header).group(1))
    height = int(re.search(r'<size_pix width="\d+" height="(\d+)"', header).group(1))
    max_value = float(re.search(r'<scanner max_value="(\d+)"', header).group(1))
    raw = np.frombuffer(part(blob, b"ImageData"), dtype="<u2").reshape(height, width)

    item = part(blob, b"ItemHeaderTag").decode("utf8")
    t = re.search(r'<transform gamma="([\d.]+)" high_frac="([\d.]+)"[^>]*low_frac="([\d.]+)"', item)
    gamma, high, low = float(t.group(1)), float(t.group(2)) * max_value, float(t.group(3)) * max_value
    # without a colormap, zero_is="white" shows high values dark
    colormap = re.search(r"<colormap>([^<]*)</colormap>", item).group(1)
    invert = ('zero_is="white"' in header and not colormap) != ('invert="true"' in item)

    v = np.clip((raw.astype(float) - low) / (high - low), 0, 1) ** gamma
    if invert:
        v = 1 - v
    return (v * 255).round().astype(np.uint8)


if __name__ == "__main__":
    Image.fromarray(render(sys.argv[1])).save(sys.argv[2])
