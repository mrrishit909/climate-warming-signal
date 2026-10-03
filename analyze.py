"""Steps 2-4: audit the two series, measure the warming, ask whether it is speeding up, and see how it lines up with carbon dioxide.

    ./venv/bin/python analyze.py
        -> results/audit.csv  decades.csv  warmest.csv  trends.csv  bandwidth.csv  acceleration.csv  slope_change.csv  co2_growth.csv  co2_link.csv  series.csv  summary.csv

Temperature: NASA GISTEMP v4 annual global land-ocean anomaly (J-D column, degrees C against the 1951-1980 mean), complete years 1880 to 2025; 2026 has 8 months and is left out.
Carbon dioxide: NOAA annual mean at Mauna Loa, 1959 to 2025 (ppm).
Every regression is ordinary least squares with a Newey-West (Bartlett kernel, 5 lags, small-sample correction n/(n-k)) covariance, because yearly temperatures are not independent; intervals use the
t distribution with n-k degrees of freedom. The five lags are fixed in advance; bandwidth.csv shows 3 and 8 as well. The time axis is in decades.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

HERE = Path(__file__).parent
R = HERE / "results"
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
LAGS = 5
WINDOWS = [(1880, 1969), (1970, 2025), (1975, 1999), (2000, 2025), (2011, 2025), (1880, 2025)]
ACCEL = [(1960, 2025), (1970, 2025), (1970, 2022), (1980, 2025)]


def ols_hac(y, X, lags=LAGS):
    """coefficients, Newey-West covariance and the t quantile for df = n - k"""
    r = sm.OLS(y, X).fit(cov_type="HAC", cov_kwds={"maxlags": lags, "use_correction": True})
    return r.params, r.cov_params(), stats.t.ppf(0.975, len(y) - X.shape[1]), sm.OLS(y, X).fit()


def main():
    R.mkdir(exist_ok=True)
    g = pd.read_csv(HERE / "data" / "gistemp.csv", skiprows=1, na_values="***").set_index("Year")
    complete = g[MONTHS].notna().all(axis=1); a = g.loc[complete, "J-D"]
    c2 = pd.read_csv(HERE / "data" / "co2_annmean_mlo.csv", comment="#").set_index("year")["mean"]
    both = a.index.intersection(c2.index)

    # Step 2: audit
    audit = [("first year of temperature", a.index.min()), ("last complete year of temperature", a.index.max()), ("annual temperature values", len(a)),
             ("months of 2026 available (left out)", int(g.loc[2026, MONTHS].notna().sum())),
             ("largest gap between the annual value and the mean of its twelve months (degrees C)", (a - g.loc[complete, MONTHS].mean(axis=1)).abs().max()),
             ("mean anomaly 1951 to 1980 (degrees C)", a.loc[1951:1980].mean()),
             ("first year of carbon dioxide", c2.index.min()), ("last year of carbon dioxide", c2.index.max()), ("years with both series", len(both))]
    pd.DataFrame(audit, columns=["check", "value"]).to_csv(R / "audit.csv", index=False)
    pd.DataFrame({"year": a.index, "anomaly": a.to_numpy(), "co2_ppm": c2.reindex(a.index).to_numpy()}).to_csv(R / "series.csv", index=False)

    # Step 3: how much warming, and how fast
    dec = a.groupby((a.index // 10) * 10).agg(years="count", mean="mean").reset_index().rename(columns={"Year": "decade"}); dec.round(6).to_csv(R / "decades.csv", index=False)
    w = a.sort_values(ascending=False, kind="stable").head(12); pd.DataFrame({"rank": range(1, 13), "year": w.index, "anomaly": w.to_numpy()}).to_csv(R / "warmest.csv", index=False)
    rows, bw = [], []
    for lo, hi in WINDOWS:
        y = a.loc[lo:hi]; x = (y.index.to_numpy() - lo) / 10.0; X = sm.add_constant(x); p, V, t, o = ols_hac(y.to_numpy(), X); se = np.sqrt(V[1, 1]); n = len(y)
        rows.append({"first": lo, "last": hi, "years": n, "per_decade": p[1], "se_naive": o.bse[1], "se_newey_west": se, "ci_lo": p[1] - t * se, "ci_hi": p[1] + t * se, "p": 2 * stats.t.sf(abs(p[1] / se), n - 2)})
        if (lo, hi) in ((1970, 2025), (1975, 1999), (2000, 2025)):
            for L in (3, 5, 8):
                p, V, t, _ = ols_hac(y.to_numpy(), X, L); s = np.sqrt(V[1, 1]); bw.append({"first": lo, "last": hi, "lags": L, "per_decade": p[1], "se_newey_west": s, "ci_lo": p[1] - t * s, "ci_hi": p[1] + t * s})
    tr = pd.DataFrame(rows); tr.round(6).to_csv(R / "trends.csv", index=False); pd.DataFrame(bw).round(6).to_csv(R / "bandwidth.csv", index=False)
    e1, e2 = tr[(tr["first"] == 1975) & (tr["last"] == 1999)].iloc[0], tr[(tr["first"] == 2000) & (tr["last"] == 2025)].iloc[0]       # two windows that do not overlap
    d = e2["per_decade"] - e1["per_decade"]; sd = np.hypot(e1["se_newey_west"], e2["se_newey_west"]); df = int(e1["years"] + e2["years"] - 4); t = stats.t.ppf(0.975, df)
    pd.DataFrame([{"later": "2000-2025", "earlier": "1975-1999", "difference_per_decade": d, "se": sd, "ci_lo": d - t * sd, "ci_hi": d + t * sd, "p": 2 * stats.t.sf(abs(d / sd), df)}]).round(6).to_csv(R / "slope_change.csv", index=False)

    # Step 4: is it speeding up? a quadratic in time; the warming rate at the start and end of the window follows from it
    rows = []
    for lo, hi in ACCEL:
        y = a.loc[lo:hi]; mid = (lo + hi) / 2; x = (y.index.to_numpy() - mid) / 10.0; X = np.column_stack([np.ones(len(x)), x, x ** 2]); p, V, t, _ = ols_hac(y.to_numpy(), X); n = len(y)
        se_c = np.sqrt(V[2, 2]); row = {"first": lo, "last": hi, "years": n, "c_per_decade_squared": p[2], "se": se_c, "ci_lo": p[2] - t * se_c, "ci_hi": p[2] + t * se_c, "p": 2 * stats.t.sf(abs(p[2] / se_c), n - 3)}
        for label, yr in (("start", lo), ("end", hi)):
            xe = (yr - mid) / 10.0; v = np.array([0.0, 1.0, 2 * xe]); s = float(v @ p); se = float(np.sqrt(v @ V @ v)); row.update({f"rate_{label}": s, f"rate_{label}_lo": s - t * se, f"rate_{label}_hi": s + t * se})
        rows.append(row)
    pd.DataFrame(rows).round(6).to_csv(R / "acceleration.csv", index=False)

    # carbon dioxide: how fast it rises, and how temperature lines up with it
    inc = c2.diff().dropna(); cg = inc.groupby((inc.index // 10) * 10).agg(years="count", mean_increase_ppm="mean").reset_index().rename(columns={"year": "decade"}); cg.round(6).to_csv(R / "co2_growth.csv", index=False)
    link = []
    for lo in (1959, 1979, 2000):
        y = a.loc[lo:2025]; x = np.log2(c2.loc[lo:2025].to_numpy()); X = sm.add_constant(x); p, V, t, o = ols_hac(y.to_numpy(), X); se = np.sqrt(V[1, 1]); n = len(y)
        link.append({"first": lo, "last": 2025, "years": n, "per_doubling": p[1], "se_newey_west": se, "ci_lo": p[1] - t * se, "ci_hi": p[1] + t * se, "r2": o.rsquared})
    pd.DataFrame(link).round(6).to_csv(R / "co2_link.csv", index=False)

    b = a.loc[1880:1920].mean(); last10 = a.loc[2016:2025].mean(); top = a.sort_values(ascending=False, kind="stable")
    summ = {"record_year": int(top.index[0]), "record_anomaly": top.iloc[0], "mean_1880_1920": b, "mean_2016_2025": last10, "warming_since_1880_1920": last10 - b,
            "eleven_warmest_are_2015_to_2025": int(sorted(top.index[:11]) == list(range(2015, 2026))), "co2_1959": c2.loc[1959], "co2_2025": c2.loc[2025], "co2_rise_pct": 100 * (c2.loc[2025] / c2.loc[1959] - 1),
            "co2_increase_2020_2025": inc.loc[2020:2025].mean()}
    pd.DataFrame(list(summ.items()), columns=["key", "value"]).round(6).to_csv(R / "summary.csv", index=False)
    for f in ("audit", "decades", "warmest", "trends", "bandwidth", "slope_change", "acceleration", "co2_growth", "co2_link", "summary"):
        print("==", f); print(pd.read_csv(R / f"{f}.csv").to_string(index=False))


if __name__ == "__main__":
    main()
