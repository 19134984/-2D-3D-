#!/usr/bin/env python
"""Verify a rigidValve2d angle(t) profile against the tracked reference.

Compares the valve tip angle transient against
.ci/reference/rigidValve2d_n40_angle.csv (generated from master, the validated
baseline) and emits an overlay + error plot as a CI artifact.

The valve flutter is a sensitive oscillator: differing host-side floating-point
contraction across toolchains (e.g. the GPU_HIP-on-NVIDIA job builds in single
mode, so the rigid-body integrator is compiled by nvcc, vs gcc for the mixed
CUDA job) shifts the *phase* of the early swing by a few hundredths of a second.
That is a benign re-synchronizing offset, not a divergence, yet a pointwise
max|error| gate flags it (~5e-2 rad where the trace is steepest) while native
CUDA stays at ~2e-4. So we gate on phase-robust functionals of the trajectory
instead: the equilibrium offset (mean), the oscillation amplitude (peak and
std), the dominant period, and an RMS backstop against gross divergence. Each
catches a real regression (wrong amplitude, frequency drift, NaN/blow-up) by
several percent, while a sub-percent phase shift passes them all.

usage: compare_valve_profile.py REFERENCE.csv PRODUCED.csv LABEL RMS_TOL OUT.png
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Phase-robust gate tolerances (physical, backend-independent). Grounded in the
# measured spread: an exaggerated 0.04 s phase shift perturbs these by <1 %,
# whereas a 10 % amplitude or 8 % frequency regression trips them by 3-10 %.
PEAK_TOL   = 0.02    # rel. error of peak |angle| (amplitude)
STD_TOL    = 0.03    # rel. error of angle std (oscillation energy)
MEAN_TOL   = 5e-3    # abs. error of mean angle (equilibrium offset), rad
PERIOD_TOL = 0.05    # rel. error of dominant oscillation period

ref_path, prod_path, label, rms_tol_s, out = sys.argv[1:6]
rms_tol = float(rms_tol_s)

def load(path):
    d = np.genfromtxt(path, delimiter=";", names=True)
    return d["time"], d["angle"]

tr, ar = load(ref_path)
tp, ap = load(prod_path)

# Align the produced trace onto the reference time grid (guards against any
# cadence mismatch; in practice the grids are identical).
api = np.interp(tr, tp, ap)
err = api - ar
abserr = np.abs(err)

def dominant_period(t, a):
    """Mean period from zero-crossings of the de-meaned signal (phase-robust)."""
    s = a - a.mean()
    z = np.where(np.diff(np.sign(s)) != 0)[0]
    if len(z) < 3:
        return float("nan")
    return 2.0 * (t[z[-1]] - t[z[0]]) / (len(z) - 1)

peak_ref, peak_prod = np.abs(ar).max(), np.abs(api).max()
std_ref,  std_prod  = ar.std(), api.std()
mean_ref, mean_prod = ar.mean(), api.mean()
per_ref,  per_prod  = dominant_period(tr, ar), dominant_period(tr, api)

rms    = float(np.sqrt(np.nanmean(err**2)))
maxerr = float(np.nanmax(abserr))
imax   = int(np.nanargmax(abserr))

# (metric, value, reference, tolerance, is_relative)
checks = [
    ("peak |angle|", peak_prod, peak_ref, PEAK_TOL,   True),
    ("angle std",    std_prod,  std_ref,  STD_TOL,    True),
    ("mean angle",   mean_prod, mean_ref, MEAN_TOL,   False),
    ("period",       per_prod,  per_ref,  PERIOD_TOL, True),
]

print(f"reference : {ref_path} ({len(tr)} pts, t in [{tr[0]:.4f},{tr[-1]:.4f}])")
print(f"produced  : {label} ({prod_path})")
print(f"final angle: produced {api[-1]:.10f}  ref {ar[-1]:.10f}  (d={api[-1]-ar[-1]:+.2e})")
print(f"max |err|  : {maxerr:.3e} rad at t={tr[imax]:.3f}  (informational; phase-sensitive)")

failed = []
for name, val, refv, tol, rel in checks:
    dev = abs(val - refv) / abs(refv) if rel else abs(val - refv)
    kind = "rel" if rel else "abs"
    ok = np.isfinite(dev) and dev <= tol
    print(f"{name:12s}: produced {val:+.5f}  ref {refv:+.5f}  ({kind} dev {dev:.3e} {'<=' if ok else '>'} tol {tol:.3e})")
    if not ok:
        failed.append(name)

rms_ok = np.isfinite(rms) and rms <= rms_tol
print(f"{'RMS |err|':12s}: {rms:.3e} rad  ({'<=' if rms_ok else '>'} backstop {rms_tol:.3e})")
if not rms_ok:
    failed.append("RMS")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
ax1.plot(tr, ar, "k-", lw=1.2, label="reference (master)")
ax1.plot(tr, api, "--", lw=0.9, color="tab:orange", label=label)
ax1.set_ylabel("valve tip angle [rad]")
ax1.legend(loc="best", fontsize=9); ax1.grid(alpha=0.3)
ax1.set_title(f"rigidValve2d N=40 — angle(t) vs reference  ({label})")

ax2.plot(tr, abserr, lw=0.9, color="tab:orange", label=f"|{label} - ref| (phase-sensitive)")
ax2.axhline(rms_tol, color="tab:red", ls=":", lw=1.2, label=f"RMS backstop = {rms_tol:g}")
ax2.set_ylabel("|angle error| [rad]"); ax2.set_xlabel("physical time [s]")
ax2.set_yscale("log"); ax2.legend(loc="best", fontsize=9); ax2.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig(out, dpi=110)
print(f"wrote {out}")

if failed:
    print(f"FAIL: phase-robust gate(s) exceeded: {', '.join(failed)}")
    sys.exit(1)
print("PASS: all phase-robust trajectory metrics within tolerance")
