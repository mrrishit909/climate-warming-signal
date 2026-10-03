"""Step 6: charts.   ./venv/bin/python charts.py -> charts/00_cover.png, 01_temperature.png, 02_trends.png, 03_acceleration.png, 04_co2.png, 05_co2_temperature.png"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
R = HERE / "results"
INK, DIM, GRID, BG, BLUE, ORANGE, GRAY, GREEN = "#f2f2f0", "#8a8a87", "#1d1d1d", "#0b0b0b", "#3987e5", "#d95926", "#9a9a96", "#3aa876"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "text.color": INK,
                     "axes.edgecolor": "#3a3a3a", "axes.labelcolor": DIM, "xtick.color": DIM, "ytick.color": DIM,
                     "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 1,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"], "font.size": 11, "axes.titlesize": 13,
                     "axes.titlelocation": "left", "axes.titlepad": 12})
S = pd.read_csv(R / "series.csv"); T = pd.read_csv(R / "trends.csv"); ACC = pd.read_csv(R / "acceleration.csv"); SM = pd.read_csv(R / "summary.csv").set_index("key")["value"]


def save(fig, name):
    fig.tight_layout(); fig.savefig(HERE / "charts" / name, dpi=150); plt.close(fig)


def temperature(ax):
    dec = pd.read_csv(R / "decades.csv")
    ax.bar(S["year"], S["anomaly"], color=[ORANGE if v > 0 else BLUE for v in S["anomaly"]], width=0.85)
    for i, r in dec.iterrows():
        ax.hlines(r["mean"], r["decade"] - 0.4, r["decade"] + r["years"] - 0.6, color=INK, linewidth=2.2, label="decade average" if i == 0 else None)
    ax.axhline(0, color=DIM, linewidth=0.8); y, v = int(SM["record_year"]), SM["record_anomaly"]
    ax.annotate(f"{y}: +{v:.2f} °C", (y, v), xytext=(1978, 1.1), color=INK, fontsize=10, arrowprops={"arrowstyle": "-", "color": INK, "linewidth": 0.8})
    ax.set_xlim(1878, 2027); ax.set_ylabel("°C"); ax.grid(axis="x", visible=False); ax.legend(frameon=False, labelcolor=INK, fontsize=9, loc="upper left")
    ax.set_title("Global temperature against the 1951 to 1980 average, by year")


def trends(ax):
    d = S[S["year"] >= 1970]; ax.scatter(d["year"], d["anomaly"], s=22, color=GRAY, alpha=0.9, linewidths=0, zorder=2)
    for (lo, hi), col in (((1975, 1999), BLUE), ((2000, 2025), ORANGE)):
        w = S[(S["year"] >= lo) & (S["year"] <= hi)]; k, b = np.polyfit(w["year"], w["anomaly"], 1); r = T[(T["first"] == lo) & (T["last"] == hi)].iloc[0]
        ax.plot([lo, hi], [k * lo + b, k * hi + b], color=col, linewidth=2.6, label=f"{lo}-{hi}: {r['per_decade']:.2f} °C per decade", zorder=3)
    w = S[S["year"] >= 1970]; c = np.polyfit(w["year"] - 1997.5, w["anomaly"], 2); xs = np.linspace(1970, 2025, 100); ax.plot(xs, np.polyval(c, xs - 1997.5), color=INK, linestyle="--", linewidth=1.5, label="quadratic fit, 1970-2025", zorder=3)
    ax.set_xlim(1969, 2027); ax.set_ylabel("°C against 1951-1980"); ax.legend(frameon=False, labelcolor=INK, fontsize=9, loc="upper left"); ax.set_title("Two straight-line trends and a quadratic")


def acceleration(ax1, ax2):
    a = ACC.copy(); a["label"] = a["first"].astype(str) + "-" + a["last"].astype(str); a = a.iloc[::-1].reset_index(drop=True); y = np.arange(len(a))
    cols = [ORANGE if lo < 0 else BLUE for lo in a["ci_lo"]]
    ax1.hlines(y, a["ci_lo"], a["ci_hi"], color=cols, linewidth=2.2); ax1.scatter(a["c_per_decade_squared"], y, s=60, color=cols, zorder=4); ax1.axvline(0, color=DIM, linewidth=1)
    ax1.set_yticks(y); ax1.set_yticklabels(a["label"]); ax1.set_ylim(-1.5, len(a) - 0.4); ax1.grid(axis="y", visible=False); ax1.set_xlabel("acceleration term, °C per decade per decade (95% interval)")
    ax1.set_title("Is there a bend in the trend?", fontsize=12)
    ax2.hlines(y, a["rate_start"], a["rate_end"], color=GRAY, linewidth=1.5); ax2.scatter(a["rate_start"], y, s=55, color=BLUE, zorder=4, label="at the start of the window"); ax2.scatter(a["rate_end"], y, s=55, color=ORANGE, zorder=4, label="at the end")
    ax2.set_yticks(y); ax2.set_yticklabels([]); ax2.set_ylim(-1.5, len(a) - 0.4); ax2.grid(axis="y", visible=False); ax2.set_xlabel("warming rate implied by the quadratic (°C per decade)"); ax2.legend(frameon=False, labelcolor=INK, fontsize=8, loc="lower right")
    ax2.set_title("How fast is it warming?", fontsize=12)


def co2(ax1, ax2):
    d = S.dropna(subset=["co2_ppm"]); ax1.plot(d["year"], d["co2_ppm"], color=BLUE, linewidth=2.4); ax1.set_ylabel("ppm"); ax1.set_xlim(1959, 2025)
    for yr in (1959, 2025):
        v = float(d.loc[d["year"] == yr, "co2_ppm"].iloc[0]); ax1.scatter([yr], [v], s=40, color=ORANGE, zorder=5); ax1.annotate(f"{v:.0f} ppm in {yr}", (yr, v), xytext=(1962, 372) if yr < 2000 else (1984, 418), color=ORANGE, fontsize=9, arrowprops={"arrowstyle": "-", "color": ORANGE, "linewidth": 0.8})
    ax1.set_title("Carbon dioxide at Mauna Loa", fontsize=12)
    g = pd.read_csv(R / "co2_growth.csv"); lab = [f"{int(d_)}s" if y == 10 else f"{int(d_)}-25" for d_, y in zip(g["decade"], g["years"])]
    ax2.bar(range(len(g)), g["mean_increase_ppm"], color=ORANGE, width=0.7)
    for i, v in enumerate(g["mean_increase_ppm"]):
        ax2.text(i, v + 0.05, f"{v:.2f}", ha="center", color=INK, fontsize=9)
    ax2.set_xticks(range(len(g))); ax2.set_xticklabels(lab, fontsize=9); ax2.set_ylabel("ppm per year"); ax2.grid(axis="x", visible=False); ax2.set_ylim(0, 3.1); ax2.set_title("Average yearly rise", fontsize=12)


def link(ax):
    d = S.dropna(subset=["co2_ppm"]); x = np.log2(d["co2_ppm"]); k, b = np.polyfit(x, d["anomaly"], 1); r = pd.read_csv(R / "co2_link.csv").iloc[0]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("t", [BLUE, ORANGE]); ax.scatter(d["co2_ppm"], d["anomaly"], c=d["year"], cmap=cmap, s=30, linewidths=0, zorder=3)
    xs = np.linspace(x.min(), x.max(), 50); ax.plot(2 ** xs, k * xs + b, color=INK, linewidth=1.4, linestyle="--", label=f"{r['per_doubling']:.2f} °C per doubling (95% interval {r['ci_lo']:.2f} to {r['ci_hi']:.2f})")
    ax.set_xscale("log", base=2); ticks = [320, 360, 400, 440]; ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks]); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for yr in (1959, 2025):
        t = d[d["year"] == yr].iloc[0]; ax.text(t["co2_ppm"] * (1 if yr == 1959 else 1.004), t["anomaly"] + (0.12 if yr == 1959 else -0.02), str(yr), color=INK, fontsize=9, ha="center" if yr == 1959 else "left")
    ax.set_xlabel("CO₂ at Mauna Loa (ppm, log scale)"); ax.set_ylabel("°C against 1951-1980"); ax.legend(frameon=False, labelcolor=INK, fontsize=9, loc="upper left"); ax.set_title("Temperature against carbon dioxide, annual means 1959 to 2025")


def main():
    (HERE / "charts").mkdir(exist_ok=True)
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 8), gridspec_kw={"width_ratios": [1.3, 1]})
    temperature(a); acceleration_ax = b
    a.set_title("Global temperature against the 1951-1980 average", fontsize=15)
    acc = ACC.copy(); acc["label"] = acc["first"].astype(str) + "-" + acc["last"].astype(str); acc = acc.iloc[::-1].reset_index(drop=True); y = np.arange(len(acc)); cols = [ORANGE if lo < 0 else BLUE for lo in acc["ci_lo"]]
    b.hlines(y, acc["ci_lo"], acc["ci_hi"], color=cols, linewidth=3); b.scatter(acc["c_per_decade_squared"], y, s=90, color=cols, zorder=4); b.axvline(0, color=DIM, linewidth=1)
    b.set_yticks(y); b.set_yticklabels(acc["label"], fontsize=12); b.set_ylim(-0.6, len(acc) - 0.4); b.grid(axis="y", visible=False); b.set_xlabel("acceleration term, °C per decade per decade (95% interval)"); b.set_title("Is there a bend in the trend?", fontsize=15)
    fig.suptitle("Is the world's warming speeding up? Likely, but the evidence thins out without 2023 to 2025", x=0.01, ha="left", fontsize=20, color=INK, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(HERE / "charts" / "00_cover.png", dpi=100); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4.0)); temperature(ax); save(fig, "01_temperature.png")
    fig, ax = plt.subplots(figsize=(10, 4.0)); trends(ax); save(fig, "02_trends.png")
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8)); acceleration(a, b); save(fig, "03_acceleration.png")
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8)); co2(a, b); save(fig, "04_co2.png")
    fig, ax = plt.subplots(figsize=(10, 4.0)); link(ax); save(fig, "05_co2_temperature.png")


if __name__ == "__main__":
    main()
