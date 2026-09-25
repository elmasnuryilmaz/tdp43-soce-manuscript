"""Export the live DOCX to Markdown (pandoc gfm) and point images at figures/main.

Embedded pictures are matched to figure files by pixel content, so the Markdown always
names the figure that the DOCX actually contains.
"""
import io, re, subprocess, zipfile
from pathlib import Path

from PIL import Image, ImageChops

ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx"
MD = ROOT / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.md"

text = subprocess.run(["pandoc", "-f", "docx", "-t", "gfm", str(DOCX)], check=True,
                      capture_output=True, text=True).stdout
z = zipfile.ZipFile(DOCX)
figures = {f: Image.open(f).convert("RGB") for f in sorted((ROOT / "figures/main").glob("Figure*.png"))}
for media in re.findall(r'src="(media/image\d+\.png)"', text):
    img = Image.open(io.BytesIO(z.read("word/" + media))).convert("RGB")
    match = [f for f, ref in figures.items() if ref.size == img.size and ImageChops.difference(ref, img).getbbox() is None]
    assert len(match) == 1, (media, match)
    text = text.replace(f'src="{media}"', f'src="../figures/main/{match[0].name}"')
MD.write_text(text, encoding="utf-8")
print(MD)
