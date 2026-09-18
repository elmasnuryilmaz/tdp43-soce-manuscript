#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Relabel the representative Fura-2 trace export in English and with the correct
CPA concentration (10 uM, confirmed by the authors). Only the text annotations are
replaced; the traces and axes are untouched."""
from PIL import Image, ImageDraw, ImageFont

SRC = "/Users/elmas/Desktop/TEZ/outputs/019ff1ef-483b-7cd0-b9df-7c3cb67ea683/Elmas_130626_Ca2_Trase_Grafikleri.png"
DST = "/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI/source_data/representative_Fura2_traces_relabelled.png"
BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
REG = "/System/Library/Fonts/Supplemental/Arial.ttf"

im = Image.open(SRC).convert("RGB")
d = ImageDraw.Draw(im)

# (box to clear, new text, font path, size, bold?)  box = (x0, y0, x1, y1)
JOBS = [
    ((225, 80, 300, 165), None, None, None),          # panel letter A
    ((1635, 15, 1715, 95), None, None, None),          # panel letter B
    ((680, 170, 940, 255), "Control", BOLD, 58),
    ((2315, 105, 2435, 175), "TDP-43 KD", BOLD, 58),
    ((790, 360, 1025, 430), "CPA (10 µM)", REG, 42),
    ((2390, 320, 2625, 390), "CPA (10 µM)", REG, 42),
    ((710, 1220, 905, 1285), "Time (s)", BOLD, 46),
    ((2285, 1290, 2465, 1355), "Time (s)", BOLD, 46),
]
for box, text, fp, size in JOBS:
    d.rectangle(box, fill="white")
    if text:
        f = ImageFont.truetype(fp, size)
        x0, y0, x1, y1 = box
        tb = d.textbbox((0, 0), text, font=f)
        tx = x0 + (x1 - x0 - (tb[2] - tb[0])) / 2 - tb[0]
        ty = y0 + (y1 - y0 - (tb[3] - tb[1])) / 2 - tb[1]
        d.text((tx, ty), text, font=f, fill="black")

im.save(DST)
print("written:", DST, im.size)
