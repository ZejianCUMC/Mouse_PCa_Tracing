from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
R = Path(__file__).resolve().parent
A = Image.open(R.parents[1] / "figures/barcode_heatmaps_v4_20261005/Np53Rb1_Intact_JZ272_barcode_heatmap.png").convert("RGB"); B = Image.open(R / "heatmap/Np53Rb1_Intact_JZ272_barcode_heatmap.png").convert("RGB")
h = 1400; A = A.resize((round(A.width * h / A.height), h)); B = B.resize((round(B.width * h / B.height), h))
c = Image.new("RGB", (A.width + B.width + 40, h + 80), "white"); d = ImageDraw.Draw(c); f = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 40)
c.paste(A, (0, 80)); c.paste(B, (A.width + 40, 80)); d.text((10, 15), "a  Current (558 cells)", fill="black", font=f); d.text((A.width + 50, 15), "b  States >= 20 cells (528 cells)", fill="black", font=f)
c.save(R / "JZ272_current_vs_min20.png"); print(R / "JZ272_current_vs_min20.png")
