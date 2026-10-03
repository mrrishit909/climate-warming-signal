"""Step 1: two public climate series.

    ./venv/bin/python download.py   -> data/gistemp.csv, data/co2_annmean_mlo.csv, results/pull.json

gistemp.csv: NASA GISS Surface Temperature Analysis (GISTEMP v4), global land-ocean temperature index: monthly and annual deviations in degrees C from the 1951-1980 average, 1880 onwards.
  Cite: GISTEMP Team, 2026: GISS Surface Temperature Analysis (GISTEMP), version 4. NASA Goddard Institute for Space Studies. Dataset accessed 2026-10-03 at https://data.giss.nasa.gov/gistemp/ ;
  and Lenssen, N., G.A. Schmidt, M. Hendrickson, P. Jacobs, M. Menne, and R. Ruedy, 2024: A GISTEMPv4 observational uncertainty ensemble. J. Geophys. Res. Atmos., 129(17), e2023JD040179.
co2_annmean_mlo.csv: NOAA Global Monitoring Laboratory, annual mean carbon dioxide at Mauna Loa, 1959 onwards (ppm), with the Scripps Institution of Oceanography record (C. D. Keeling) before 1974.
  Freely available; NOAA asks for a citation (Lan, X., Tans, P. and K.W. Thoning, Trends in globally-averaged CO2 determined from NOAA Global Monitoring Laboratory measurements).
Both files are tiny and are updated by their providers, so the exact copies used are committed, with their checksums.
"""
import hashlib
import json
import ssl
import urllib.request
from datetime import date
from pathlib import Path

import certifi

HERE = Path(__file__).parent
FILES = {"gistemp.csv": "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv",
         "co2_annmean_mlo.csv": "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_mlo.csv"}
CTX = ssl.create_default_context(cafile=certifi.where())


def main():
    (HERE / "data").mkdir(exist_ok=True)
    (HERE / "results").mkdir(exist_ok=True)
    info = {"pulled": date.today().isoformat(), "files": {}}
    for name, url in FILES.items():
        blob = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "research-portfolio (Rishit Raj Mathur)"}), context=CTX, timeout=120).read()
        (HERE / "data" / name).write_bytes(blob)
        info["files"][name] = {"url": url, "bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest()}
    (HERE / "results" / "pull.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
