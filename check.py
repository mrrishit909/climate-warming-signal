"""Step 5: the project's check. Recompute with different code, in plain Python (lists and loops: no pandas, numpy or statsmodels): the audit from the raw text of both files, the decade means, the
ranking of the warmest years, every regression by a hand-written least squares with a hand-written Newey-West covariance (Bartlett weights, small-sample correction), the intervals and p-values from
the t distribution, the acceleration terms and the warming rate at the start and end of each window, the carbon dioxide growth by decade and the link between temperature and CO2.

    ./venv/bin/python check.py
"""
import csv
import math
from pathlib import Path

from scipy import stats

HERE = Path(__file__).parent
R = HERE / "results"
LAGS = 5
WINDOWS = [(1880, 1969), (1970, 2025), (1975, 1999), (2000, 2025), (2011, 2025), (1880, 2025)]
ACCEL = [(1960, 2025), (1970, 2025), (1970, 2022), (1980, 2025)]
n_ok = 0


def same(label, got, want, abs_tol=2e-6):
    global n_ok
    n_ok += 1
    assert abs(float(got) - float(want)) <= abs_tol, f"{label}: recount {got} != results {want}"


def table(name):
    return list(csv.DictReader(open(R / f"{name}.csv")))


def inverse(M):
    k = len(M); A = [row[:] + [1.0 if i == j else 0.0 for j in range(k)] for i, row in enumerate(M)]
    for c in range(k):
        p = max(range(c, k), key=lambda r: abs(A[r][c])); A[c], A[p] = A[p], A[c]
        d = A[c][c]; A[c] = [v / d for v in A[c]]
        for r in range(k):
            if r != c:
                f = A[r][c]; A[r] = [v - f * w for v, w in zip(A[r], A[c])]
    return [row[k:] for row in A]


def mm(A, B):
    return [[sum(A[i][t] * B[t][j] for t in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def ols_hac(y, X, L=LAGS):
    """coefficients, Newey-West covariance, naive standard errors, R squared"""
    n, k = len(y), len(X[0]); XtX = [[sum(X[t][i] * X[t][j] for t in range(n)) for j in range(k)] for i in range(k)]; I = inverse(XtX)
    beta = [sum(I[i][j] * sum(X[t][j] * y[t] for t in range(n)) for j in range(k)) for i in range(k)]
    e = [y[t] - sum(X[t][j] * beta[j] for j in range(k)) for t in range(n)]; xe = [[X[t][j] * e[t] for j in range(k)] for t in range(n)]
    S = [[0.0] * k for _ in range(k)]
    for t in range(n):
        for i in range(k):
            for j in range(k):
                S[i][j] += xe[t][i] * xe[t][j]
    for lag in range(1, L + 1):
        w = 1 - lag / (L + 1)
        for t in range(lag, n):
            for i in range(k):
                for j in range(k):
                    S[i][j] += w * (xe[t][i] * xe[t - lag][j] + xe[t - lag][i] * xe[t][j])
    V = [[v * n / (n - k) for v in row] for row in mm(mm(I, S), I)]
    sse = sum(v * v for v in e); my = sum(y) / n; sst = sum((v - my) ** 2 for v in y)
    return beta, V, [math.sqrt(sse / (n - k) * I[j][j]) for j in range(k)], 1 - sse / sst


def main():
    temp, months_present = {}, {}
    lines = open(HERE / "data" / "gistemp.csv").read().splitlines(); header = lines[1].split(","); assert header[0] == "Year" and header[13] == "J-D"
    for line in lines[2:]:
        f = line.split(",")
        if len(f) < 14 or not f[0].isdigit():
            continue
        m = [None if v == "***" else float(v) for v in f[1:13]]; months_present[int(f[0])] = sum(v is not None for v in m)
        if all(v is not None for v in m):
            assert f[13] != "***"; temp[int(f[0])] = (float(f[13]), sum(m) / 12)
    co2 = {int(r[0]): float(r[1]) for r in csv.reader(l for l in open(HERE / "data" / "co2_annmean_mlo.csv") if not l.startswith("#") and l.strip()) if r[0].isdigit()}
    a = {y: v[0] for y, v in temp.items()}; years = sorted(a)

    # audit and series
    au = {r["check"]: float(r["value"]) for r in table("audit")}
    same("first year", years[0], au["first year of temperature"]); same("last complete year", years[-1], au["last complete year of temperature"]); same("annual values", len(years), au["annual temperature values"])
    same("months of 2026", months_present[2026], au["months of 2026 available (left out)"]); assert 2026 not in a and years == list(range(1880, 2026))
    same("gap", max(abs(v[0] - v[1]) for v in temp.values()), au["largest gap between the annual value and the mean of its twelve months (degrees C)"])
    same("baseline mean", sum(a[y] for y in range(1951, 1981)) / 30, au["mean anomaly 1951 to 1980 (degrees C)"])
    same("co2 first", min(co2), au["first year of carbon dioxide"]); same("co2 last", max(co2), au["last year of carbon dioxide"]); same("years with both", len(set(co2) & set(a)), au["years with both series"])
    for r in table("series"):
        y = int(r["year"]); same(f"anomaly {y}", a[y], r["anomaly"])
        assert (r["co2_ppm"] == "") == (y not in co2)
        if y in co2:
            same(f"co2 {y}", co2[y], r["co2_ppm"])

    # decades and the warmest years
    for r in table("decades"):
        d = int(r["decade"]); ys = [y for y in years if y // 10 * 10 == d]; same(f"decade {d} years", len(ys), r["years"]); same(f"decade {d} mean", sum(a[y] for y in ys) / len(ys), r["mean"])
    top = sorted(years, key=lambda y: (-a[y], y))[:12]
    for r, y in zip(table("warmest"), top):
        same(f"rank {r['rank']} year", y, r["year"]); same("anomaly", a[y], r["anomaly"])

    # trends, with the Newey-West covariance
    def line_fit(lo, hi, L):
        ys = [a[y] for y in range(lo, hi + 1)]; X = [[1.0, (y - lo) / 10.0] for y in range(lo, hi + 1)]; b, V, se0, _ = ols_hac(ys, X, L); n = len(ys); t = stats.t.ppf(0.975, n - 2); s = math.sqrt(V[1][1])
        return {"n": n, "slope": b[1], "naive": se0[1], "se": s, "lo": b[1] - t * s, "hi": b[1] + t * s, "p": 2 * stats.t.sf(abs(b[1] / s), n - 2)}
    tr = {}
    for r in table("trends"):
        lo, hi = int(r["first"]), int(r["last"]); f = line_fit(lo, hi, LAGS); tr[(lo, hi)] = f
        same(f"years {lo}-{hi}", f["n"], r["years"]); same("slope", f["slope"], r["per_decade"]); same("naive se", f["naive"], r["se_naive"]); same("nw se", f["se"], r["se_newey_west"])
        same("ci_lo", f["lo"], r["ci_lo"]); same("ci_hi", f["hi"], r["ci_hi"]); same("p", f["p"], r["p"])
    assert set(tr) == set(WINDOWS)
    for r in table("bandwidth"):
        f = line_fit(int(r["first"]), int(r["last"]), int(r["lags"])); same(f"bandwidth {r['first']}-{r['last']} L{r['lags']} slope", f["slope"], r["per_decade"])
        same("se", f["se"], r["se_newey_west"]); same("ci_lo", f["lo"], r["ci_lo"]); same("ci_hi", f["hi"], r["ci_hi"])
    s = table("slope_change")[0]; e1, e2 = tr[(1975, 1999)], tr[(2000, 2025)]; d = e2["slope"] - e1["slope"]; sd = math.hypot(e1["se"], e2["se"]); df = e1["n"] + e2["n"] - 4; t = stats.t.ppf(0.975, df)
    same("slope change", d, s["difference_per_decade"]); same("se", sd, s["se"]); same("ci_lo", d - t * sd, s["ci_lo"]); same("ci_hi", d + t * sd, s["ci_hi"]); same("p", 2 * stats.t.sf(abs(d / sd), df), s["p"])

    # acceleration: a quadratic in time
    rows = table("acceleration"); assert [(int(r["first"]), int(r["last"])) for r in rows] == ACCEL
    for r in rows:
        lo, hi = int(r["first"]), int(r["last"]); mid = (lo + hi) / 2; ys = [a[y] for y in range(lo, hi + 1)]; xs = [(y - mid) / 10.0 for y in range(lo, hi + 1)]
        b, V, _, _ = ols_hac(ys, [[1.0, x, x * x] for x in xs]); n = len(ys); t = stats.t.ppf(0.975, n - 3); sc = math.sqrt(V[2][2])
        same(f"accel {lo}-{hi} years", n, r["years"]); same("c", b[2], r["c_per_decade_squared"]); same("se", sc, r["se"]); same("ci_lo", b[2] - t * sc, r["ci_lo"]); same("ci_hi", b[2] + t * sc, r["ci_hi"])
        same("p", 2 * stats.t.sf(abs(b[2] / sc), n - 3), r["p"])
        for label, yr in (("start", lo), ("end", hi)):
            xe = (yr - mid) / 10.0; v = [0.0, 1.0, 2 * xe]; est = sum(v[i] * b[i] for i in range(3)); se = math.sqrt(sum(v[i] * V[i][j] * v[j] for i in range(3) for j in range(3)))
            same(f"rate_{label}", est, r[f"rate_{label}"]); same(f"rate_{label}_lo", est - t * se, r[f"rate_{label}_lo"]); same(f"rate_{label}_hi", est + t * se, r[f"rate_{label}_hi"])

    # carbon dioxide
    inc = {y: co2[y] - co2[y - 1] for y in sorted(co2) if y - 1 in co2}
    for r in table("co2_growth"):
        d = int(r["decade"]); ys = [y for y in inc if y // 10 * 10 == d]; same(f"co2 decade {d} years", len(ys), r["years"]); same("mean increase", sum(inc[y] for y in ys) / len(ys), r["mean_increase_ppm"])
    for r in table("co2_link"):
        lo = int(r["first"]); ys = [a[y] for y in range(lo, 2026)]; X = [[1.0, math.log2(co2[y])] for y in range(lo, 2026)]; b, V, _, r2 = ols_hac(ys, X); n = len(ys); t = stats.t.ppf(0.975, n - 2); s = math.sqrt(V[1][1])
        same(f"co2 link {lo} years", n, r["years"]); same("per doubling", b[1], r["per_doubling"]); same("se", s, r["se_newey_west"]); same("ci_lo", b[1] - t * s, r["ci_lo"]); same("ci_hi", b[1] + t * s, r["ci_hi"]); same("r2", r2, r["r2"])

    sm = {r["key"]: float(r["value"]) for r in table("summary")}
    b0 = sum(a[y] for y in range(1880, 1921)) / 41; l10 = sum(a[y] for y in range(2016, 2026)) / 10
    exp = {"record_year": top[0], "record_anomaly": a[top[0]], "mean_1880_1920": b0, "mean_2016_2025": l10, "warming_since_1880_1920": l10 - b0, "eleven_warmest_are_2015_to_2025": int(sorted(top[:11]) == list(range(2015, 2026))),
           "co2_1959": co2[1959], "co2_2025": co2[2025], "co2_rise_pct": 100 * (co2[2025] / co2[1959] - 1), "co2_increase_2020_2025": sum(inc[y] for y in range(2020, 2026)) / 6}
    assert sorted(sm) == sorted(exp)
    for k, v in exp.items():
        same(f"summary {k}", v, sm[k])
    print(f"check.py: {n_ok} recomputed numbers agree with results/ (hand-written least squares and Newey-West covariance)")


if __name__ == "__main__":
    main()
