#!/usr/bin/env python3
"""VETTAM unit test suite.

Builds and runs the VETTAM verification problems, compares each against its
analytic solution, and exits non-zero if any tolerance is exceeded -- so a
change to the module can be checked with one command.

    ./run_tests.py                     # quick suite (~2 min incl. builds)
    ./run_tests.py --plots             # also write diagnostic figures
    ./run_tests.py --full              # adds the expensive advection test
    ./run_tests.py --test diffusion    # one test
    ./run_tests.py --compile-only      # build, do not run
    ./run_tests.py --list              # show the registry

FLASH is located by walking up from this file until a `setup` script is found,
or via --flash-dir / the FLASH_DIR environment variable. Parfiles come from
parfiles/ next to this script.

Only the VETTAM-based tests are covered. Marshak_Wave, Radiation_Shock_LE,
RHDThinShell, Rosen_SinkTest and Rosen_ShadowBlobs all select the legacy
physics/RadTrans/RadTransMain/VET or Cool/CoolMain/RadTransVET unit instead, so
they are not registered here and currently have no driver -- porting them means
rewriting their simulation-level cool_vet.F90 / cool_hydro_update.F90 overrides
against VETTAM's rt_* interface.

Requires numpy and h5py. Reference: Menon et al. 2022, MNRAS 512, 401.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import vettam_io as vio  # noqa: E402

C_LIGHT = 2.99792458e10          # cm/s
A_RAD = 7.565733250033928e-15    # erg cm^-3 K^-4
R_GAS = 8.314462618e7            # erg K^-1 mol^-1

SETUP_COMMON = ["-1d", "-auto", "-nxb=8", "+bouchut", "-maxblocks=500"]


# --------------------------------------------------------------------------
# Checks. Each returns (passed, {metric: value}, [message, ...]).
# --------------------------------------------------------------------------

def _check(results, label, value, limit, fmt="{:.3e}"):
    ok = bool(value <= limit)
    results.append(
        ("PASS" if ok else "FAIL")
        + f"  {label}: " + fmt.format(value) + " (limit " + fmt.format(limit) + ")"
    )
    return ok


def check_streaming(case_dir, tol):
    """Gaussian pulse must advect at exactly c with the energy conserved."""
    files = vio.plotfiles(case_dir)
    t0, x0, dx0, d0 = vio.read(files[0])
    e0 = vio.integral(d0["erad"], dx0)
    box = x0.max() - x0.min() + dx0[0]
    nu = 20.0

    metrics, msgs, ok = {}, [], True
    for idx in (2, 5, 10):
        t, x, dx, d = vio.read(files[idx])
        erad = d["erad"]
        # analytic: the initial Gaussian rigidly shifted by c*t, periodic
        shift = (x - C_LIGHT * t + box / 2.0) % box - box / 2.0
        analytic = np.exp(-(nu ** 2) * shift ** 2) + erad.min()

        drift = abs(vio.integral(erad, dx) - e0) / e0
        # position error, measured in cells, accounting for periodicity
        x_num = x[np.argmax(erad)]
        x_an = (C_LIGHT * t + box / 2.0) % box - box / 2.0
        dpos = abs((x_num - x_an + box / 2.0) % box - box / 2.0) / dx[0]
        l1 = vio.l1_error(erad, analytic, dx)

        ncross = t * C_LIGHT / box
        metrics[f"t{idx}"] = {"crossings": ncross, "energy_drift": drift,
                              "position_err_cells": dpos, "L1": l1}
        ok &= _check(msgs, f"[{ncross:.1f} crossings] energy drift", drift, tol["energy"])
        ok &= _check(msgs, f"[{ncross:.1f} crossings] pulse position (cells)",
                     dpos, tol["position_cells"], "{:.2f}")
        ok &= _check(msgs, f"[{ncross:.1f} crossings] L1 vs analytic",
                     l1, tol["l1"][idx])
    return ok, metrics, msgs


def check_diffusion(case_dir, tol, sigma_t):
    """Pulse spreading must follow the analytic diffusion solution."""
    files = vio.plotfiles(case_dir)
    diff = C_LIGHT / (3.0 * sigma_t)
    nu = 20.0
    _, _, dx0, d0 = vio.read(files[0])
    e0 = vio.integral(d0["erad"], dx0)

    metrics, msgs, ok = {}, [], True
    for idx in (4, 16, 50, 100):
        t, x, dx, d = vio.read(files[idx])
        erad = d["erad"]
        spread = 4.0 * diff * t * nu ** 2 + 1.0
        analytic = np.exp(-(nu ** 2) * x ** 2 / spread) / np.sqrt(spread) + erad.min()

        drift = abs(vio.integral(erad, dx) - e0) / e0
        l1 = vio.l1_error(erad, analytic, dx)
        sqrt_dt = np.sqrt(diff * t)
        metrics[f"t{idx}"] = {"sqrt_Dt": sqrt_dt, "energy_drift": drift, "L1": l1}
        ok &= _check(msgs, f"[sqrt(Dt)={sqrt_dt:.2f}] energy drift", drift, tol["energy"])
        ok &= _check(msgs, f"[sqrt(Dt)={sqrt_dt:.2f}] L1 vs analytic", l1, tol["l1"])
    return ok, metrics, msgs


def check_diffusion_width(case_dir, tol, sigma_t):
    """Pulse *width* must track the analytic diffusion rate.

    The second moment measures how much physical time the radiation field
    actually saw, so this is what catches the field being advanced by the wrong
    interval when VETTAM's update gate skips hydro steps (BUG-10). It is also
    insensitive to the peak clipping that dominates the L1 error at large steps,
    which is why the tolerance can be tight.
    """
    files = vio.plotfiles(case_dir)
    diff = C_LIGHT / (3.0 * sigma_t)
    nu = 20.0
    _, _, dx0, d0 = vio.read(files[0])
    e0 = vio.integral(d0["erad"], dx0)

    metrics, msgs, ok = {}, [], True
    for idx in (50, 100):
        t, x, dx, d = vio.read(files[idx])
        erad = d["erad"]
        # the uniform floor (a_r T_min^4) carries no information and would
        # dominate the second moment over a 10 cm box, so take it out
        signal = erad - erad.min()
        rms = float(np.sqrt(np.sum(signal * x ** 2 * dx) / np.sum(signal * dx)))
        rms_analytic = float(np.sqrt(1.0 / (2 * nu ** 2) + 2.0 * diff * t))
        err = abs(rms - rms_analytic) / rms_analytic
        drift = abs(vio.integral(erad, dx) - e0) / e0

        sqrt_dt = np.sqrt(diff * t)
        metrics[f"t{idx}"] = {"sqrt_Dt": sqrt_dt, "rms": rms,
                              "rms_analytic": rms_analytic, "rms_err": err,
                              "energy_drift": drift}
        ok &= _check(msgs, f"[sqrt(Dt)={sqrt_dt:.2f}] pulse width "
                           f"({rms:.5f} vs {rms_analytic:.5f})", err, tol["width"])
        ok &= _check(msgs, f"[sqrt(Dt)={sqrt_dt:.2f}] energy drift", drift, tol["energy"])
    return ok, metrics, msgs


def _lowrie_edwards_reference(sim_dir):
    """Dimensional Lowrie & Edwards semi-analytic shock structure.

    M03.0_soln.out holds six dimensionless rows (x, rho, v, T_gas, T_rad, F_r);
    the normalisation below is the one Simulation_init applies.
    """
    sigma_a, kappa, p0, gamma, box = 1.0e6, 1.0, 1.0e-4, 5.0 / 3.0, 1.0
    c_s = np.sqrt(C_LIGHT ** 2 / (3.0 * sigma_a * kappa))
    t_norm = c_s ** 2 / (gamma * R_GAS)
    x, _, _, t_gas, t_rad, _ = np.loadtxt(
        os.path.join(sim_dir, "M03.0_soln.out"))
    return x * box, t_gas * t_norm, t_rad * t_norm


def check_marshak(case_dir, tol, sim_dir):
    """Marshak wave against the Su & Olson (1996) tables.

    The dimensionless time is tau = epsilon * c * rho*kappa * t. Output lands
    where the time stepping puts it rather than exactly on the tabulated taus,
    so the closest dump to each is used and its real tau reported. The tau=0.1
    table is measured but not asserted on: the earliest dump sits at tau~0.117
    and there is no closer one, so the comparison there is dominated by that
    17% time offset rather than by solver error.
    """
    epsilon, sigma_t = 0.1, 40.0
    files = vio.plotfiles(case_dir)
    times = np.array([vio.read(f)[0] for f in files])
    taus = epsilon * C_LIGHT * sigma_t * times

    metrics, msgs, ok = {}, [], True
    for target in (0.1, 1, 10, 100):
        idx = int(np.argmin(np.abs(taus - target)))
        t, x, dx, d = vio.read(files[idx], fields=("erad", "temp"))
        # the tables hold U = E_r and V = a_r T^4, both already dimensionless
        v_sim = A_RAD * d["temp"] ** 4
        table = np.loadtxt(os.path.join(sim_dir, f"tau_{target}.dat"), skiprows=1)
        x_ref, u_ref, v_ref = table[:, 1], table[:, 2], table[:, 3]

        inside = (x >= x_ref.min()) & (x <= x_ref.max())
        l1_u = vio.l1_error(d["erad"][inside], np.interp(x, x_ref, u_ref)[inside], dx[inside])
        l1_v = vio.l1_error(v_sim[inside], np.interp(x, x_ref, v_ref)[inside], dx[inside])
        offset = abs(taus[idx] - target) / target

        metrics[f"tau{target}"] = {"tau_actual": taus[idx], "tau_offset": offset,
                                   "L1_U": l1_u, "L1_V": l1_v}
        if offset > tol["tau_offset"]:
            msgs.append(f"SKIP  [tau={target}] nearest dump is tau={taus[idx]:.3f}"
                        f" ({offset:.0%} off); L1 U={l1_u:.3f} V={l1_v:.3f} not asserted")
            continue
        ok &= _check(msgs, f"[tau={taus[idx]:.2f}] L1 radiation energy", l1_u, tol["l1"])
        ok &= _check(msgs, f"[tau={taus[idx]:.2f}] L1 material energy", l1_v, tol["l1"])
    return ok, metrics, msgs


def check_radshock(case_dir, tol, sim_dir):
    """Gas and radiation temperature against Lowrie & Edwards (2008), M = 3.

    The radiation moments are stored in the lab frame, so they are transformed
    to the comoving frame before a radiation temperature is taken from them.
    """
    files = vio.plotfiles(case_dir)
    x_an, tgas_an, trad_an = _lowrie_edwards_reference(sim_dir)

    metrics, msgs, ok = {}, [], True
    for idx in (len(files) // 2, len(files) - 1):
        t, x, dx, d = vio.read(files[idx],
                               fields=("erad", "dens", "temp", "mohx", "velx"))
        erad = d["erad"] - 2.0 * (d["velx"] / C_LIGHT ** 2) * d["mohx"] \
            * np.sign(d["velx"]) * np.sign(d["mohx"])
        t_rad = (erad / A_RAD) ** 0.25

        inside = (x >= x_an.min()) & (x <= x_an.max())
        gas_ref = np.interp(x, x_an, tgas_an)
        rad_ref = np.interp(x, x_an, trad_an)
        l1_gas = vio.l1_error(d["temp"][inside], gas_ref[inside], dx[inside])
        l1_rad = vio.l1_error(t_rad[inside], rad_ref[inside], dx[inside])

        metrics[f"t{idx}"] = {"time": t, "L1_Tgas": l1_gas, "L1_Trad": l1_rad}
        ok &= _check(msgs, f"[t={t:.2e}] L1 gas temperature", l1_gas, tol["l1"])
        ok &= _check(msgs, f"[t={t:.2e}] L1 radiation temperature", l1_rad, tol["l1"])
    return ok, metrics, msgs


def check_radtempeq(case_dir, tol, rho, mu=0.6, gamma=5.0 / 3.0,
                    t_init=100.0, e_init=7.565733250033928e-13):
    """Gas and radiation must relax to the energy-conserving equilibrium.

    a_r T^4 + rho R T / ((gamma - 1) mu) = e_total, solved for T.
    """
    files = vio.plotfiles(case_dir)
    t, x, dx, d = vio.read(files[-1], fields=("erad", "temp"))

    c_v = R_GAS * rho / ((gamma - 1.0) * mu)
    e_total = c_v * t_init + e_init
    roots = np.roots([A_RAD, 0.0, 0.0, c_v, -e_total])
    t_eq = float(next(r.real for r in roots if abs(r.imag) < 1e-8 and r.real > 0))

    t_gas = float(d["temp"].mean())
    t_rad = float((d["erad"].mean() / A_RAD) ** 0.25)
    err_gas = abs(t_gas - t_eq) / t_eq
    err_rad = abs(t_rad - t_eq) / t_eq

    msgs, ok = [], True
    metrics = {"T_equilibrium": t_eq, "T_gas": t_gas, "T_rad": t_rad,
               "err_gas": err_gas, "err_rad": err_rad, "t_final": t}
    ok &= _check(msgs, f"gas temperature ({t_gas:.4f} K vs {t_eq:.4f} K)",
                 err_gas, tol["temperature"])
    ok &= _check(msgs, f"radiation temperature ({t_rad:.4f} K vs {t_eq:.4f} K)",
                 err_rad, tol["temperature"])
    return ok, metrics, msgs


def check_advection(case_dir, tol):
    """Energy conservation only -- the reference solution is semi-analytic."""
    files = vio.plotfiles(case_dir)
    _, _, dx0, d0 = vio.read(files[0])
    e0 = vio.integral(d0["erad"], dx0)
    _, _, dx, d = vio.read(files[-1])
    drift = abs(vio.integral(d["erad"], dx) - e0) / e0
    msgs = []
    ok = _check(msgs, "energy drift", drift, tol["energy"])
    return ok, {"energy_drift": drift}, msgs


# --------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------

TESTS = {
    "streaming": {
        "simulation": "Streaming_Pulse",
        "objdir": "objVETTAM_Streaming",
        "quick": True,
        "cases": {"streaming": "RadPulseStreaming.par"},
        "check": {"streaming": lambda d: check_streaming(
            d, {"energy": 1e-6, "position_cells": 2.0,
                "l1": {2: 0.12, 5: 0.25, 10: 0.40}})},
        "blurb": "free-streaming pulse; transport limit (Fig. 1 left)",
    },
    "diffusion": {
        "simulation": "Diffusion_Pulse",
        "objdir": "objVETTAM_Diffusion",
        "quick": True,
        "cases": {"weak": "RadPulseWeakDiffusion.par",
                  "strong": "RadPulseStrongDiffusion.par",
                  "skipgate": "RadPulseSkipGate.par"},
        "check": {
            "weak": lambda d: check_diffusion(d, {"energy": 1e-6, "l1": 0.06}, 40.0),
            "strong": lambda d: check_diffusion(d, {"energy": 1e-6, "l1": 0.06}, 4.0e4),
            # Skipping ~4 of every 5 hydro steps must not change the diffusion
            # rate. Observed 2-3% with the fix, ~53% without it (BUG-10).
            "skipgate": lambda d: check_diffusion_width(
                d, {"energy": 2e-6, "width": 0.10}, 40.0),
        },
        "blurb": "diffusing pulse: weak, strong, and update-gate skipping",
    },
    "marshak": {
        # Su & Olson Marshak wave. Ported off the legacy Cool unit: its private
        # cool_hydro_update.F90 existed only to add a closed-form temperature
        # update for c_v = alpha*T^3, which is now rt_temp_type=5 with
        # rt_cv_alpha. Simulation_initBlock supplies alpha = 4*a_r/sim_epsilon.
        "simulation": "Marshak_Wave",
        "objdir": "objVETTAM_Marshak",
        "quick": True,
        "cases": {"suolson": "MarshakWave.par"},
        "check": {"suolson": lambda d: check_marshak(
            d, {"l1": 0.06, "tau_offset": 0.05}, str(HERE / "Marshak_Wave"))},
        "blurb": "Marshak wave vs Su & Olson (Fig. 4)",
    },
    "radshock": {
        # Lowrie & Edwards (2008) M=3 radiation shock. Ported off the legacy
        # Cool unit: it needed per-face Dirichlet radiation temperatures and a
        # constant total opacity, both of which VETTAM now supports, so the
        # 2841-line cool_vet.F90 override it used to carry is gone.
        "simulation": "Radiation_Shock_LE",
        "objdir": "objVETTAM_RadShock",
        "quick": True,
        "cases": {"m3": "RadnShock.par"},
        "check": {"m3": lambda d: check_radshock(
            d, {"l1": 0.03},
            str(HERE / "Radiation_Shock_LE"))},
        "blurb": "radiation shock vs Lowrie & Edwards (Fig. 5)",
    },
    "radtempeq": {
        "simulation": "RadTempEq",
        "objdir": "objVETTAM_RadTempEq",
        "quick": True,
        # One 8-cell block, so extra ranks own no blocks and the PETSc solve
        # stalls. This problem is 0-D; a single rank is the right way to run it.
        "max_nprocs": 1,
        "cases": {"p1e-9": "RadTempEq_P1e-9.par",
                  "p10": "RadTempEq_P10.par",
                  "p1e11": "RadTempEq_P1e11.par"},
        "check": {
            "p1e-9": lambda d: check_radtempeq(d, {"temperature": 1e-3}, 1e-7),
            "p10": lambda d: check_radtempeq(d, {"temperature": 1e-3}, 1e-17),
            "p1e11": lambda d: check_radtempeq(d, {"temperature": 1e-3}, 1e-27),
        },
        "blurb": "radiation-gas energy equilibration (Fig. 3)",
    },
    "radtempeq_lin": {
        # Same problem and same parfiles as "radtempeq", but built with the
        # LinTempUpdate overrides, which replace the outer Picard/Newton-Raphson
        # temperature iteration with a fully coupled linearised solve: gas
        # temperature becomes an extra matrix unknown (NDIM+2 per cell) with the
        # a_r T^4 emission linearised about T^n. rt_hydro_type and rt_temp_type
        # are inert in that build -- rt_dustTemperature only re-runs the EOS.
        #
        # Both schemes must land on the same equilibrium, since that is fixed by
        # energy conservation; they differ during the relaxation (~1% at the
        # first dump for P = 10), so this is an independent discretisation
        # checked against the same analytic answer, not a rerun.
        "simulation": "RadTempEq",
        "objdir": "objVETTAM_RadTempEqLin",
        "setup_extra": ["-unit=physics/RadTrans/RadTransMain/VETTAM/LinTempUpdate"],
        "quick": True,
        "max_nprocs": 1,
        "cases": {"p1e-9": "RadTempEq_P1e-9.par",
                  "p10": "RadTempEq_P10.par",
                  "p1e11": "RadTempEq_P1e11.par"},
        "check": {
            "p1e-9": lambda d: check_radtempeq(d, {"temperature": 1e-3}, 1e-7),
            "p10": lambda d: check_radtempeq(d, {"temperature": 1e-3}, 1e-17),
            "p1e11": lambda d: check_radtempeq(d, {"temperature": 1e-3}, 1e-27),
        },
        "blurb": "rad-gas equilibration, linearised temperature update",
    },
    "advection": {
        "simulation": "Advection_Pulse",
        "objdir": "objVETTAM_Advection",
        "quick": False,
        "cases": {"vel0": "Advection_Pulse_vel0.par",
                  "vel10": "Advection_Pulse_vel10.par"},
        "check": {"vel0": lambda d: check_advection(d, {"energy": 1e-4}),
                  "vel10": lambda d: check_advection(d, {"energy": 1e-4})},
        "blurb": "advecting pulse, dynamic diffusion (Fig. 6) -- SLOW, hours",
    },
}


# --------------------------------------------------------------------------
# Build / run
# --------------------------------------------------------------------------

def find_flash(explicit=None):
    if explicit:
        p = Path(explicit).resolve()
    elif os.environ.get("FLASH_DIR"):
        p = Path(os.environ["FLASH_DIR"]).resolve()
    else:
        p = HERE
        while p != p.parent and not (p / "setup").is_file():
            p = p.parent
    if not (p / "setup").is_file():
        sys.exit(f"error: no FLASH 'setup' script found (looked at {p}).\n"
                 "       Pass --flash-dir or set FLASH_DIR.")
    return p


def run(cmd, cwd, log):
    with open(log, "w") as fh:
        return subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)


def missing_units(flash, test):
    """Units named in setup_extra that are absent from this FLASH tree.

    A test may depend on an optional sub-unit. Where that is not part of the
    installation -- VETTAM is distributed without LinTempUpdate, for instance --
    the test is skipped rather than reported as a failure.
    """
    missing = []
    for arg in test.get("setup_extra", []):
        if arg.startswith("-unit="):
            unit = arg.split("=", 1)[1]
            if not (flash / "source" / unit).is_dir():
                missing.append(unit)
    return missing


def build(flash, test, jobs, verbose):
    objdir = flash / test["objdir"]
    cmd = ["./setup", f"unitTest/VETTAM/{test['simulation']}"] + SETUP_COMMON + \
          test.get("setup_extra", []) + [f"-objdir={test['objdir']}"]
    print(f"    setup: {' '.join(cmd)}")
    if run(cmd, flash, flash / f"setup_{test['objdir']}.log") != 0:
        return None, f"setup failed (see setup_{test['objdir']}.log)"
    print(f"    make -j{jobs}")
    if run(["make", f"-j{jobs}"], objdir, objdir / "make.log") != 0:
        return None, f"build failed (see {test['objdir']}/make.log)"
    return objdir, None


def stale(case_dir, objdir, parfile):
    """Why the stored output for this case cannot be trusted, or None.

    Reusing output produced by an older binary or an older parfile is the single
    easiest way to get a confidently wrong answer out of this suite -- it is how
    a stale objdir once made the diffusion test look like it had a 1500% error.
    So compare mtimes rather than trusting the caller to remember --overwrite.
    """
    outputs = vio.plotfiles(str(case_dir))
    if not outputs:
        return "no output"
    newest = max(Path(f).stat().st_mtime for f in outputs)
    binary = objdir / "flash4"
    if binary.exists() and binary.stat().st_mtime > newest:
        return "flash4 is newer than the stored output"
    if Path(parfile).stat().st_mtime > newest:
        return "the parfile is newer than the stored output"
    return None


def simulate(objdir, parfile, case, nprocs, overwrite):
    case_dir = objdir / case
    if case_dir.exists() and not overwrite:
        why = stale(case_dir, objdir, parfile)
        if why is None:
            print(f"    {case}: reusing existing output (--overwrite to rerun)")
            return case_dir, None
        if why != "no output":
            print(f"    {case}: rerunning -- {why}")
    shutil.rmtree(case_dir, ignore_errors=True)
    case_dir.mkdir(parents=True)

    shutil.copy2(parfile, objdir / "flash.par")
    for leftover in list(objdir.glob("*_hdf5_plt_cnt_*")) + list(objdir.glob("*_hdf5_chk_*")):
        leftover.unlink()

    t0 = time.time()
    rc = run(["./flash4"], objdir, objdir / f"run_{case}.log")
    if rc != 0:
        return case_dir, f"FLASH exited {rc} (see {objdir.name}/run_{case}.log)"
    for f in list(objdir.glob("*_hdf5_plt_cnt_*")) + list(objdir.glob("*_hdf5_chk_*")):
        shutil.move(str(f), case_dir / f.name)
    print(f"    {case}: ran in {time.time() - t0:.1f}s on {nprocs} rank(s)")
    return case_dir, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flash-dir")
    ap.add_argument("--test", action="append", choices=sorted(TESTS),
                    help="run only this test (repeatable)")
    ap.add_argument("--full", action="store_true", help="include the slow tests")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--compile-only", action="store_true")
    ap.add_argument("--no-compile", action="store_true")
    ap.add_argument("--overwrite", action="store_true", help="rerun even if output exists")
    ap.add_argument("-np", "--nprocs", type=int, default=4)
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--metrics", default=str(HERE / "unitTest_metrics.json"))
    ap.add_argument("--plots", action="store_true",
                    help="write diagnostic figures (needs matplotlib)")
    ap.add_argument("--plot-dir", default=str(HERE / "figures"),
                    help="where --plots writes its PDFs (default: ./figures)")
    args = ap.parse_args()

    if args.list:
        for name, t in sorted(TESTS.items()):
            print(f"  {name:15s} {'quick' if t['quick'] else 'slow ':5s}  {t['blurb']}")
        return 0

    selected = args.test or [n for n, t in TESTS.items() if t["quick"] or args.full]
    flash = find_flash(args.flash_dir)
    parfiles = HERE / "parfiles"
    print(f"FLASH tree: {flash}")
    print(f"Tests:      {', '.join(selected)}\n")

    metrics, failures, skipped = {}, [], []
    for name in selected:
        test = TESTS[name]
        print(f"[{name}] {test['blurb']}")
        absent = missing_units(flash, test)
        if absent:
            print(f"    SKIP  not installed: {', '.join(absent)}\n")
            skipped.append(f"{name}: {', '.join(absent)} not installed")
            continue
        objdir = flash / test["objdir"]
        if not args.no_compile:
            objdir, err = build(flash, test, args.jobs, False)
            if err:
                print(f"    ERROR {err}\n"); failures.append(f"{name}: {err}"); continue
        elif not (objdir / "flash4").exists():
            err = f"{objdir} has no flash4 binary; drop --no-compile"
            print(f"    ERROR {err}\n"); failures.append(f"{name}: {err}"); continue
        if args.compile_only:
            print("    built (compile-only)\n"); continue

        metrics[name] = {}
        produced = {}
        for case, par in test["cases"].items():
            nprocs = min(args.nprocs, test.get("max_nprocs", args.nprocs))
            case_dir, err = simulate(objdir, parfiles / par, case, nprocs, args.overwrite)
            if err:
                print(f"    ERROR {case}: {err}"); failures.append(f"{name}/{case}: {err}")
                continue
            ok, m, msgs = test["check"][case](str(case_dir))
            metrics[name][case] = m
            produced[case] = case_dir
            for line in msgs:
                print("      " + line)
            if not ok:
                failures.append(f"{name}/{case}: tolerance exceeded")

        if args.plots and produced:
            try:
                import vettam_plots
                out = Path(args.plot_dir) / f"{name}.pdf"
                vettam_plots.PLOTS[name](produced, out)
                print(f"    figure: {out}")
            except ImportError as exc:
                print(f"    figure skipped: {exc} (matplotlib not installed)")
            except Exception as exc:                      # noqa: BLE001
                # A broken plot must not turn a passing suite red.
                print(f"    figure FAILED: {type(exc).__name__}: {exc}")
        print()

    if metrics:
        with open(args.metrics, "w") as fh:
            json.dump(metrics, fh, indent=2, default=float)
        print(f"metrics written to {args.metrics}")

    print("\n" + "=" * 66)
    for s_ in skipped:
        print(f"SKIPPED  {s_}")
    if failures:
        print(f"FAILED ({len(failures)}):")
        for f in failures:
            print("  - " + f)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
