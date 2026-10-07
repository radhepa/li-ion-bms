"""
Download the LG 18650HG2 test files used by this project (25 C, ~9 MB) from
Mendeley Data: P. Kollmeyer et al., McMaster University, "LG 18650HG2 Li-ion
Battery Data", doi:10.17632/cp3473x7xv.3, licensed CC BY 4.0.
"""
import time
import urllib.request
from pathlib import Path

BASE = "https://data.mendeley.com/public-files/datasets/cp3473x7xv/files/{}/file_downloaded"
FILES = {
    "HPPC_549.mat": "f7542bd4-73ee-4fff-a3a8-e6f2ff0489ed",
    "C20_549.mat": "cf4a2ca6-b724-4e55-ad28-bc0fea30456c",
    "Cap1C_551.mat": "45a10363-90f3-433a-b54b-cda5989ee70f",
    "UDDS_551.csv": "58f73e68-6748-4232-a26a-176b3dab0c1f",
    "LA92_551.mat": "dec0b398-41c0-451b-ba0f-614b2059452b",
    "US06_551.mat": "5c8e6b67-974e-41f9-a898-1c8da4b85193",
    "Mixed1_551.mat": "105b96e7-39f9-4bf2-8eac-40b6f8b13784",
}
RAW = Path(__file__).resolve().parent / "raw"
RAW.mkdir(exist_ok=True)

for name, fid in FILES.items():
    dest = RAW / name
    if dest.exists() and dest.stat().st_size > 10_000:
        print(f"have  {name}")
        continue
    for attempt in range(5):           # the server sometimes answers with a transient error
        data = urllib.request.urlopen(BASE.format(fid), timeout=120).read()
        if len(data) > 10_000:
            dest.write_bytes(data)
            print(f"got   {name} ({len(data) / 1e6:.1f} MB)")
            break
        time.sleep(3)
    else:
        print(f"FAILED {name}: try again later")
