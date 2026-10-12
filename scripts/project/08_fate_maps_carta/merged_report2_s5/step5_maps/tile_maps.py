#!/usr/bin/env python
# Side-by-side sheet of the six cohort maps: rows Npp53 / Npp53Rb1 / both; columns c10 included / excluded.
import sys
from pathlib import Path
from PIL import Image
H = Path(__file__).resolve().parent / 'maps'; suf = sys.argv[1] if len(sys.argv) > 1 else ''
rows = ['Npp53', 'Npp53Rb1', 'Npp53_Npp53Rb1']
ims = [[Image.open(H / f'{r}_{v}{suf}.png').convert('RGB') for v in ('inc', 'exc')] for r in rows]
w = max(i.width for row in ims for i in row); h = max(i.height for row in ims for i in row)
sheet = Image.new('RGB', (2 * w, 3 * h), 'white')
for ri, row in enumerate(ims):
    for ci, im in enumerate(row): sheet.paste(im, (ci * w, ri * h))
sheet.thumbnail((2 * w // 2, 3 * h // 2))
sheet.save(H / f'six_maps_sheet{suf}.png'); print(sheet.size)
