"""Swap the page-1 gel image of an Image Lab report PDF for a cleaned one."""
import io
import sys

from PIL import Image
from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, NumberObject


def replace(src_pdf, clean_png, dst_pdf):
    img = Image.open(clean_png).convert("L")
    buf = io.BytesIO()
    img.save(buf, "JPEG", quality=95, subsampling=0)

    writer = PdfWriter(clone_from=PdfReader(src_pdf))
    xobjects = writer.pages[0]["/Resources"]["/XObject"]
    gels = [x.get_object() for x in xobjects.values()
            if x.get_object()["/Subtype"] == "/Image"]
    assert len(gels) == 1, f"expected one image on page 1, found {len(gels)}"
    gel = gels[0]
    assert (gel["/Width"], gel["/Height"]) == img.size

    gel._data = buf.getvalue()
    gel[NameObject("/Filter")] = NameObject("/DCTDecode")
    gel[NameObject("/ColorSpace")] = NameObject("/DeviceGray")
    gel[NameObject("/BitsPerComponent")] = NumberObject(8)
    gel[NameObject("/Length")] = NumberObject(len(gel._data))
    gel.pop("/DecodeParms", None)

    with open(dst_pdf, "wb") as f:
        writer.write(f)


if __name__ == "__main__":
    replace(*sys.argv[1:4])
