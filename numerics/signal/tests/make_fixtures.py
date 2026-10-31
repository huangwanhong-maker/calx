"""Reference values for calx's signal checks, from scipy.signal.

Writes tests/fixtures/spectra-reference.txt: one line per array, a key and then its numbers,

    key v0 v1 v2 …

read by check_tapers.cpp and check_spectra.cpp. Comment lines start with '#'.

- dpss/<n>/<nw>/taper<k>: the stride s, then h_0, h_s, h_2s, … of scipy.signal.windows.dpss(n, nw, K)
  (norm 2: unit energy; scipy's signs), K = 2 nw − 1; every sample for n ≤ 64, every 4th for n = 1000,
  every 16th for n = 4096, every 1000th for n = 100000.
- dpss/<n>/<nw>/ratios: the concentrations, dpss(…, return_ratios=True).
- series/<name>: input series (numpy's PCG64, seeded), so that C++ reads exactly what scipy read.
- periodogram/<series>/<variant>: scipy.signal.periodogram(x, fs = 10, window, detrend, scaling='density')
  — boxcar, boxcar without detrending, the periodic Hann window (calx's hann(n)) — for n = 1000 and 999.
- csd/<L>/{pxx, pyy, re, im, coh}: scipy.signal.welch, csd and coherence of series x4000 and y4000
  (y a lagged, filtered copy of x plus noise), fs = 2, Hann, nperseg = L, noverlap = L // 2
  (calx's step round(L/2) for even L, and for L = 255 its step 128 = 255 − 127).

    .venv/Scripts/python.exe -I foundations/calx/numerics/signal/tests/make_fixtures.py
"""

import datetime
from pathlib import Path

import numpy as np
import scipy
from scipy import signal
from scipy.signal.windows import dpss, get_window

OUT = Path(__file__).resolve().parent / "fixtures" / "spectra-reference.txt"


def row(key, values):
    return key + " " + " ".join(f"{float(v):.17g}" for v in values)


def main() -> None:
    lines = [
        f"# calx signal: reference values from scipy {scipy.__version__} (numpy {np.__version__}), {datetime.date.today().isoformat()}.",
        "# Written by make_fixtures.py; each line is a key and its numbers.",
    ]
    # ---- discrete prolate spheroidal sequences ----
    for n in (64, 1000, 4096, 100000):
        for nw in (2, 4):
            if n == 100000 and nw == 2:
                continue
            K = 2 * nw - 1
            h, ratios = dpss(n, nw, K, return_ratios=True)
            stride = 1 if n <= 64 else 4 if n == 1000 else 16 if n == 4096 else 1000
            for k in range(K):
                lines.append(row(f"dpss/{n}/{nw}/taper{k}", [stride] + list(h[k, ::stride])))
            lines.append(row(f"dpss/{n}/{nw}/ratios", ratios))

    # ---- series ----
    rng = np.random.default_rng(20261030)
    fs = 10.0
    series = {}
    for n in (1000, 999):
        t = np.arange(n) / fs
        series[f"x{n}"] = 1.5 + 2.0 * np.sin(2 * np.pi * 1.25 * t + 0.4) + 0.7 * np.cos(2 * np.pi * 3.1 * t) + rng.standard_normal(n)
    n = 4000
    x = rng.standard_normal(n + 50)
    xf = np.convolve(x, [0.5, 1.0, 0.5], mode="same")  # a low-pass colouring
    y = 0.8 * np.roll(xf, 3) + 0.6 * rng.standard_normal(n + 50)  # y lags x by 3 samples
    series["x4000"] = xf[50:]
    series["y4000"] = y[50:]
    for name, v in series.items():
        lines.append(row(f"series/{name}", v))

    # ---- periodograms ----
    for name in ("x1000", "x999"):
        v = series[name]
        m = len(v)
        _, p = signal.periodogram(v, fs, window="boxcar", detrend="constant", scaling="density")
        lines.append(row(f"periodogram/{name}/boxcar", p))
        _, p = signal.periodogram(v, fs, window="boxcar", detrend=False, scaling="density")
        lines.append(row(f"periodogram/{name}/nodemean", p))
        w = get_window("hann", m)  # periodic (DFT-even), as calx's hann
        _, p = signal.periodogram(v, fs, window=w, detrend="constant", scaling="density")
        lines.append(row(f"periodogram/{name}/hann", p))

    # ---- cross spectra and coherence ----
    fs2 = 2.0
    for L in (256, 255):
        kw = dict(fs=fs2, window="hann", nperseg=L, noverlap=L // 2, detrend="constant", scaling="density")
        _, pxx = signal.welch(series["x4000"], **kw)
        _, pyy = signal.welch(series["y4000"], **kw)
        _, pxy = signal.csd(series["x4000"], series["y4000"], **kw)
        _, coh = signal.coherence(series["x4000"], series["y4000"], fs=fs2, window="hann", nperseg=L, noverlap=L // 2)
        lines.append(row(f"csd/{L}/pxx", pxx))
        lines.append(row(f"csd/{L}/pyy", pyy))
        lines.append(row(f"csd/{L}/re", pxy.real))
        lines.append(row(f"csd/{L}/im", pxy.imag))
        lines.append(row(f"csd/{L}/coh", coh))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({len(lines)} lines, {OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
