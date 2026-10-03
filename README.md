# Is the world's warming speeding up?

Global temperature since 1880 and carbon dioxide since 1959, from two public sources. Four questions:
- how much has the world warmed, and how recent are the warmest years;
- how fast is it warming now, against 1975 to 1999 and against the years before 1970, with uncertainty that allows for one year resembling the last;
- is the rate itself rising, and does that depend on which years are included;
- how fast is carbon dioxide rising, and how closely does temperature follow it.

The questions follow GeeksforGeeks' "Climate Change Data Analysis". The data is real and current (to 2025), not a course file.

Build log: https://mrrishit909.github.io/projects/climate-warming-signal/

## Data

- **Temperature:** NASA GISS Surface Temperature Analysis (GISTEMP v4), the global land-ocean index: monthly and annual deviations in degrees C from the 1951-1980 average, 1880 onwards. Complete years 1880 to 2025 (146);
  2026 has 8 months and is left out. Cite: GISTEMP Team, 2026: GISS Surface Temperature Analysis (GISTEMP), version 4. NASA Goddard Institute for Space Studies. Dataset accessed 2026-10-03 at https://data.giss.nasa.gov/gistemp/;
  and Lenssen, N., G.A. Schmidt, M. Hendrickson, P. Jacobs, M. Menne and R. Ruedy, 2024: A GISTEMPv4 observational uncertainty ensemble. J. Geophys. Res. Atmos., 129(17), e2023JD040179.
- **Carbon dioxide:** NOAA Global Monitoring Laboratory, annual mean at Mauna Loa, 1959 to 2025 (ppm), with the Scripps Institution of Oceanography record (C. D. Keeling) before 1974. Freely available; NOAA asks for a citation.
  Measurements at Mauna Loa were suspended from November 2022 to July 2023 after an eruption, and the record for December 2022 to July 2023 is from a site on Maunakea.

Both files are tiny and their providers update them, so the exact copies used are committed, with checksums in `results/pull.json`; `download.py` fetches the current versions, which can differ after the next update. The script re-fetched the NOAA file byte for byte; the GISS server was refusing connections when I tried it on the temperature file, so that request has not been run end to end by the script (the file was fetched by hand from the same address).

```
./venv/bin/pip install pandas numpy scipy statsmodels matplotlib certifi      # in a venv
./venv/bin/python download.py    # data/*.csv
./venv/bin/python analyze.py     # results/*
./venv/bin/python check.py
./venv/bin/python charts.py
```

## Steps

1. **Download** the two series.
2. **Audit:** complete years, whether each annual value matches the mean of its twelve months, and whether the 1951-1980 anomalies average zero, as the baseline says.
3. **How much, how fast:** decade averages, the warmest years, and straight-line trends in six windows by least squares with a Newey-West covariance (Bartlett kernel, 5 lags, small-sample correction), because yearly
   temperatures are not independent. Intervals use the t distribution; the lags are fixed in advance, and 3 and 8 are shown as well. The 2000-2025 slope is compared with the 1975-1999 slope.
4. **Is it speeding up?** A quadratic in time over four windows; its curvature term, and the warming rate it implies at the start and end of each window.
5. **Carbon dioxide:** the average yearly rise by decade, and temperature regressed on log2 of the concentration, which gives degrees C per doubling.
6. **Check.** `check.py` recomputes 449 numbers in plain Python (lists and loops, without pandas, numpy or statsmodels): the audit from the raw text, the decade means and rankings, every regression with a hand-written
   least squares and a hand-written Newey-West covariance, the intervals and p-values, the acceleration terms and rates, and the carbon dioxide figures. The hand-written covariance matches statsmodels' to the
   last digit saved, and changing its weights or its small-sample correction makes the check fail.
7. **Charts.**

## Results

- **How much.** The last ten years (2016 to 2025) averaged 1.02 degrees C above the 1951-1980 average and 1.30 above the 1880-1920 average (-0.28), a common stand-in for pre-industrial that is not an official baseline. The eleven warmest
  years in the record are the eleven years 2015 to 2025: 2024 is the warmest at +1.29, then 2025 (+1.19) and 2023 (+1.17). Decade averages: -0.22 in the 1880s, -0.34 in the 1910s (the coldest), 0.04 in the 1970s, 0.25, 0.39, 0.59 and 0.81
  in the 1980s to 2010s, and 1.07 for 2020 to 2025.
- **How fast.** 0.037 degrees C per decade from 1880 to 1969 (95% interval 0.022 to 0.052); 0.202 from 1970 to 2025 (0.179 to 0.225); 0.170 from 1975 to 1999 (0.119 to 0.221) and 0.267 from 2000 to 2025 (0.207 to 0.327); 0.386 in
  the last 15 years (0.291 to 0.482). The 2000-2025 rate is 0.097 above the 1975-1999 rate (0.020 to 0.173, p = 0.014). The ordinary standard error is too small: for 1970 to 2025 the Newey-West one is 0.0115 against 0.0088, 31% larger.
  With 3 or 8 lags instead of 5 the intervals barely move (2000 to 2025: 0.203 to 0.331, 0.207 to 0.327 and 0.213 to 0.321).
- **Speeding up, but it depends on three years.** Fitting a quadratic to 1970 to 2025, the curvature is 0.0166 degrees C per decade per decade (0.0059 to 0.0272, p = 0.003): the warming rate rises from 0.11 (0.05 to 0.17) in 1970
  to 0.29 (0.23 to 0.36) in 2025. Starting in 1960 the curvature is 0.0185 (0.0111 to 0.0259) and starting in 1980 it is 0.0328 (0.0164 to 0.0492). Stopping in 2022, before the three record years, it is 0.0077, positive but with
  an interval that includes zero (-0.0019 to 0.0174, p = 0.11), and the rate at the end is 0.23 (0.17 to 0.29) instead of 0.29.
- **Carbon dioxide.** 316.0 ppm in 1959 and 427.4 in 2025, 35.2% more. It rose by 0.86 ppm a year on average in the 1960s, 1.22, 1.64, 1.53, 1.91 and 2.40 in the next decades, and 2.62 a year from 2020 to 2025, three times the 1960s rate.
- **How closely they line up.** Temperature rises 2.74 degrees C per doubling of carbon dioxide over 1959 to 2025 (2.49 to 2.99, R squared 0.93); 2.84 over 1979 to 2025 (2.50 to 3.19); 3.24 over 2000 to 2025 (2.58 to 3.91).

## Not done, and caveats

- One temperature dataset. Other groups publish global series, and the annual measurement uncertainty that GISTEMP provides is not carried into the intervals here; it is small beside the year-to-year variation, but not zero.
- The intervals describe the scatter of years around a smooth trend, allowing for dependence by a Newey-West covariance, which is an approximation; windows of 15 to 26 years give wide intervals that move with the window. The windows
  are round-number choices and all are reported; the 1970-2022 window was added to see what the record years do, and it changes the answer.
- The quadratic is a convenient shape, not a physical model; it cannot tell a steady bend from a step. The record years of 2023 to 2025 include the 2023-24 El Niño, and El Niño, like the 2015-16 event, adds warmth to a year; it is not removed
  here, so part of the recent rise in the rate may be that and not a change in the underlying trend.
- Degrees per doubling is an association between two series that both trend, not climate sensitivity: other greenhouse gases, particles in the air and ocean cycles changed over the same years, and the number depends on the window.
- Mauna Loa is one remote site, and the last year's value can be revised. 2026 is left out because it has 8 months.
