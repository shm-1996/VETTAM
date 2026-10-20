"""Diagnostic figures for the VETTAM unit tests.

Ported from the plotting content of Menon21_Plots.py (Menon et al. 2022, MNRAS
512, 401) onto the numpy + h5py reader in vettam_io, so no yt, scipy or
flashplotlib is needed. matplotlib is the only extra dependency and is imported
lazily -- run_tests.py only touches this module when --plots is given.

Each figure overlays the numerical solution on its analytic counterpart and is
annotated with the measured L1 error next to the value published in the method
paper, so agreement (or drift) is visible in the figure itself.
"""

import numpy as np

import vettam_io as vio

C_LIGHT = 2.99792458e10
A_RAD = 7.565733250033928e-15
R_GAS = 8.314462618e7

# Fixed categorical order, never cycled. These are the method paper's colours,
# kept for continuity with the published figures; they also clear the
# colour-vision checks comfortably (worst CVD dE = 19.1 against a target of 8,
# worst normal-vision dE = 27.9 against a floor of 15). Series identity is
# additionally carried by line style (numerical solid, analytic dashed), so it
# never rests on colour alone.
SERIES = ["#a2191f", "#009d9a", "#6929c4"]

NUM = dict(ls="-", lw=1.4, alpha=0.85)        # numerical
ANA = dict(ls="--", lw=1.1, alpha=0.9)        # analytic

# L1 errors reported in Menon et al. 2022, for context in the annotations.
PUBLISHED = {
    "streaming": [0.037, 0.16, 0.277],
    "weak": [0.05, 0.037, 0.033],
    "strong": [0.024, 0.03, 0.032],
}


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "xtick.top": False,
        "ytick.right": False,
        "figure.dpi": 150,
        # mathtext, not a local LaTeX install: keeps the figures reproducible
        # and stops stray characters being treated as LaTeX specials
        "text.usetex": False,
    })
    return plt


# Annotations sit on top of the curves, so give them a surface to rest on.
NOTE = dict(color="0.35", fontsize=7.5,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.5))


def _finish(fig, outfile):
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return outfile


# --------------------------------------------------------------------------
# Analytic solutions (same expressions the checkers in run_tests.py use)
# --------------------------------------------------------------------------

def analytic_streaming(x, t, box, nu=20.0, floor=0.0):
    shift = (x - C_LIGHT * t + box / 2.0) % box - box / 2.0
    return np.exp(-(nu ** 2) * shift ** 2) + floor


def analytic_diffusion(x, t, sigma_t, nu=20.0, floor=0.0):
    spread = 4.0 * (C_LIGHT / (3.0 * sigma_t)) * t * nu ** 2 + 1.0
    return np.exp(-(nu ** 2) * x ** 2 / spread) / np.sqrt(spread) + floor


def equilibrium_temperature(rho, mu=0.6, gamma=5.0 / 3.0,
                            t_init=100.0, e_init=7.565733250033928e-13):
    c_v = R_GAS * rho / ((gamma - 1.0) * mu)
    roots = np.roots([A_RAD, 0.0, 0.0, c_v, -(c_v * t_init + e_init)])
    return float(next(r.real for r in roots if abs(r.imag) < 1e-8 and r.real > 0))


# --------------------------------------------------------------------------
# Figures, one per registered test
# --------------------------------------------------------------------------

def plot_streaming(case_dirs, outfile):
    """E_r(x) at three crossing times, against the rigidly advected profile."""
    plt = _mpl()
    files = vio.plotfiles(str(case_dirs["streaming"]))
    _, x0, dx0, _ = vio.read(files[0])
    box = x0.max() - x0.min() + dx0[0]

    fig, ax = plt.subplots(figsize=(5.0, 3.6), tight_layout=True)
    errs = []
    for colour, idx in zip(SERIES, (2, 5, 10)):
        t, x, dx, d = vio.read(files[idx])
        erad = d["erad"]
        an = analytic_streaming(x, t, box, floor=erad.min())
        ncross = t * C_LIGHT / box
        ax.plot(x, erad, color=colour, label=rf"$ct/L = {ncross:.1f}$", **NUM)
        ax.plot(x, an, color=colour, **ANA)
        errs.append(vio.l1_error(erad, an, dx))

    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel(r"$x\ (\mathrm{cm})$")
    ax.set_ylabel(r"$E_r\ (\mathrm{erg\ cm^{-3}})$")
    ax.set_title("Free-streaming pulse", loc="left")
    ax.text(0.02, 0.97, r"$\rho\kappa = 10^{-20}\ \mathrm{cm^{-1}}$",
            transform=ax.transAxes, va="top", color="0.35")
    ax.text(0.02, 0.06,
            "L1  " + ", ".join(f"{e:.3f}" for e in errs) +
            "\nMenon+22  " + ", ".join(f"{e:.3f}" for e in PUBLISHED["streaming"]),
            transform=ax.transAxes, **NOTE)
    ax.plot([], [], color="0.4", label="numerical", **NUM)
    ax.plot([], [], color="0.4", label="analytic", **ANA)
    ax.legend(loc="upper right", fontsize=7.5)
    return _finish(fig, outfile)


def plot_diffusion(case_dirs, outfile):
    """E_r(x) for each diffusion case against the spreading Gaussian."""
    plt = _mpl()
    panels = [("weak", 40.0, r"$\rho\kappa = 40$"),
              ("strong", 4.0e4, r"$\rho\kappa = 4\times10^{4}$"),
              ("skipgate", 40.0, r"$\rho\kappa = 40$, gate skipping")]
    panels = [p for p in panels if p[0] in case_dirs]

    fig, axs = plt.subplots(ncols=len(panels), figsize=(4.4 * len(panels), 3.6),
                            sharey=True, tight_layout=True)
    axs = np.atleast_1d(axs)

    for ax, (case, sigma_t, blurb) in zip(axs, panels):
        files = vio.plotfiles(str(case_dirs[case]))
        errs = []
        # For the skip-gate case the first gated radiation update lands after
        # plotfile 4, so that output is still essentially the initial condition
        # and comparing it to the analytic solution says nothing. Start later.
        indices = (16, 50, 100) if case == "skipgate" else (4, 16, 100)
        for colour, idx in zip(SERIES, indices):
            t, x, dx, d = vio.read(files[idx])
            erad = d["erad"]
            an = analytic_diffusion(x, t, sigma_t, floor=erad.min())
            sqrt_dt = np.sqrt(C_LIGHT / (3.0 * sigma_t) * t)
            ax.plot(x, erad, color=colour, label=rf"$\sqrt{{Dt}} = {sqrt_dt:.2f}$", **NUM)
            ax.plot(x, an, color=colour, **ANA)
            errs.append(vio.l1_error(erad, an, dx))

        ax.set_xlim(-5, 5)
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 1.5)
        ax.set_xlabel(r"$x\ (\mathrm{cm})$")
        ax.set_title(blurb, loc="left")
        note = "L1  " + ", ".join(f"{e:.3f}" for e in errs)
        if case in PUBLISHED:
            note += "\nMenon+22  " + ", ".join(f"{e:.3f}" for e in PUBLISHED[case])
        ax.text(0.02, 0.06, note, transform=ax.transAxes, **NOTE)
        ax.legend(loc="upper right", fontsize=7.5)

    axs[0].set_ylabel(r"$E_r\ (\mathrm{erg\ cm^{-3}})$")
    return _finish(fig, outfile)


def plot_radtempeq(case_dirs, outfile):
    """Gas and radiation temperature relaxing to the equilibrium value.

    The x axis is logarithmic: all three cases relax within a couple of thermal
    times and then sit flat out to t/t_therm ~ 12, so a linear axis spends most
    of its width on straight lines.
    """
    plt = _mpl()
    cases = [("p1e-9", 1e-7, r"$\mathcal{P} = 10^{-9}$"),
             ("p10", 1e-17, r"$\mathcal{P} = 10$"),
             ("p1e11", 1e-27, r"$\mathcal{P} = 10^{11}$")]
    cases = [c for c in cases if c[0] in case_dirs]
    t_therm = 1.0 / C_LIGHT

    fig, ax = plt.subplots(figsize=(5.4, 3.8), tight_layout=True)
    t_end = 1.0
    for colour, (case, rho, label) in zip(SERIES, cases):
        times, t_gas, t_rad = [], [], []
        for f in vio.plotfiles(str(case_dirs[case])):
            t, _, _, d = vio.read(f, fields=("erad", "temp"))
            times.append(t)
            t_gas.append(d["temp"].mean())
            t_rad.append((d["erad"].mean() / A_RAD) ** 0.25)
        times = np.array(times) / t_therm
        t_gas, t_rad = np.array(t_gas), np.array(t_rad)

        # t = 0 cannot go on a log axis; the first dump after it is the start.
        keep = times > 0
        times, t_gas, t_rad = times[keep], t_gas[keep], t_rad[keep]
        t_end = max(t_end, times[-1])

        t_eq = equilibrium_temperature(rho)
        ax.axhline(t_eq, color=colour, lw=0.8, alpha=0.4)
        ax.plot(times, t_gas, color=colour, **NUM)
        ax.plot(times, t_rad, color=colour, **ANA)
        # Direct label on the plateau, inside the axes.
        ax.text(times[-1], t_eq * 1.12, label, color=colour, fontsize=8,
                ha="right", va="bottom")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(right=t_end)
    ax.set_xlabel(r"$t/t_\mathrm{therm}$")
    ax.set_ylabel(r"$T\ (\mathrm{K})$")
    ax.set_title("Radiation-gas energy equilibration", loc="left")

    # Line-style key only; identity is carried by the direct labels above.
    ax.plot([], [], color="0.4", label=r"$T_\mathrm{gas}$", **NUM)
    ax.plot([], [], color="0.4", label=r"$T_\mathrm{rad}$", **ANA)
    ax.plot([], [], color="0.4", lw=0.8, alpha=0.4, label="equilibrium")
    ax.legend(loc="center left", fontsize=7.5)
    return _finish(fig, outfile)


def plot_advection(case_dirs, outfile):
    """Static vs advected pulse, with the moving case shifted back into frame."""
    plt = _mpl()
    v_adv = 1.0e6  # cm/s, the vel10 case

    f0 = vio.plotfiles(str(case_dirs["vel0"]))[-1]
    f1 = vio.plotfiles(str(case_dirs["vel10"]))[-1]
    _, x, dx, a = vio.read(f0, fields=("dens", "temp", "velx"))
    t1, _, _, b = vio.read(f1, fields=("dens", "temp", "velx"))

    box = x.max() - x.min() + dx[0]
    shift_cells = int(round((v_adv * t1 % box) / dx[0]))
    b = {k: np.roll(v, -shift_cells) for k, v in b.items()}

    fig, axs = plt.subplots(nrows=3, figsize=(5.4, 8.0), sharex=True, tight_layout=True)
    for ax, key, ylab in ((axs[0], "dens", r"$\rho\ (\mathrm{g\ cm^{-3}})$"),
                          (axs[1], "temp", r"$T\ (\mathrm{K})$")):
        ax.plot(x, a[key], color=SERIES[0], label=r"$v = 0$", **NUM)
        ax.plot(x, b[key], color=SERIES[1], label=r"$v = 10\ \mathrm{km\ s^{-1}}$", **ANA)
        ax.set_ylabel(ylab)
    axs[0].legend(loc="upper right", fontsize=7.5)
    axs[0].set_title("Advecting pulse (dynamic diffusion)", loc="left")

    for colour, key, lab in ((SERIES[0], "dens", r"$\rho$"), (SERIES[1], "temp", r"$T$")):
        axs[2].plot(x, (a[key] - b[key]) / a[key], color=colour, label=lab, **NUM)
    axs[2].set_ylabel("relative error")
    axs[2].set_xlabel(r"$x\ (\mathrm{cm})$")
    axs[2].legend(loc="upper right", fontsize=7.5)

    l1 = {k: vio.l1_error(b[k], a[k], dx) for k in ("dens", "temp")}
    axs[2].text(0.02, 0.06,
                f"L1  rho {l1['dens']:.4f}, T {l1['temp']:.4f}"
                "\nMenon+22  rho 0.0014, T 0.0017",
                transform=axs[2].transAxes, **NOTE)
    return _finish(fig, outfile)


def plot_radshock(case_dirs, outfile):
    """Gas and radiation temperature against the Lowrie & Edwards profile."""
    plt = _mpl()
    import os
    sigma_a, kappa, p0, gamma, box = 1.0e6, 1.0, 1.0e-4, 5.0 / 3.0, 1.0
    c_s = np.sqrt(C_LIGHT ** 2 / (3.0 * sigma_a * kappa))
    t_norm = c_s ** 2 / (gamma * R_GAS)
    ref = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "Radiation_Shock_LE", "M03.0_soln.out")
    x_an, _, _, tgas_an, trad_an, _ = np.loadtxt(ref)
    x_an, tgas_an, trad_an = x_an * box, tgas_an * t_norm, trad_an * t_norm

    files = vio.plotfiles(str(case_dirs["m3"]))
    t, x, dx, d = vio.read(files[-1],
                           fields=("erad", "dens", "temp", "mohx", "velx"))
    erad = d["erad"] - 2.0 * (d["velx"] / C_LIGHT ** 2) * d["mohx"] \
        * np.sign(d["velx"]) * np.sign(d["mohx"])
    t_rad = (erad / A_RAD) ** 0.25

    fig, ax = plt.subplots(figsize=(5.4, 3.8), tight_layout=True)
    ax.plot(x_an, tgas_an, color="0.35", **ANA)
    ax.plot(x_an, trad_an, color="0.35", **ANA)
    ax.plot(x, d["temp"], color=SERIES[0], label=r"$T_\mathrm{gas}$", **NUM)
    ax.plot(x, t_rad, color=SERIES[1], label=r"$T_\mathrm{rad}$", **NUM)
    ax.plot([], [], color="0.35", label="Lowrie and Edwards", **ANA)

    inside = (x >= x_an.min()) & (x <= x_an.max())
    l1g = vio.l1_error(d["temp"][inside], np.interp(x, x_an, tgas_an)[inside], dx[inside])
    l1r = vio.l1_error(t_rad[inside], np.interp(x, x_an, trad_an)[inside], dx[inside])
    ax.text(0.02, 0.92, f"L1  gas {l1g:.4f}, rad {l1r:.4f}",
            transform=ax.transAxes, **NOTE)

    ax.set_xlim(x_an.min(), x_an.max())
    ax.set_xlabel(r"$x\ (\mathrm{cm})$")
    ax.set_ylabel(r"$T\ (\mathrm{K})$")
    ax.set_title("Radiation shock, $M = 3$", loc="left")
    ax.legend(loc="lower right", fontsize=7.5)
    return _finish(fig, outfile)


def plot_marshak(case_dirs, outfile):
    """Marshak wave at four dimensionless times against Su & Olson."""
    plt = _mpl()
    import os
    epsilon, sigma_t = 0.1, 40.0
    ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Marshak_Wave")

    files = vio.plotfiles(str(case_dirs["suolson"]))
    times = np.array([vio.read(f)[0] for f in files])
    taus = epsilon * C_LIGHT * sigma_t * times

    fig, axs = plt.subplots(ncols=2, figsize=(9.0, 3.8), sharey=False, tight_layout=True)
    colours = SERIES + ["#0043ce"]
    for colour, target in zip(colours, (0.1, 1, 10, 100)):
        idx = int(np.argmin(np.abs(taus - target)))
        t, x, dx, d = vio.read(files[idx], fields=("erad", "temp"))
        v_sim = A_RAD * d["temp"] ** 4
        table = np.loadtxt(os.path.join(ref_dir, f"tau_{target}.dat"), skiprows=1)
        x_ref, u_ref, v_ref = table[:, 1], table[:, 2], table[:, 3]

        axs[0].plot(x, d["erad"], color=colour, label=rf"$\tau = {target:g}$", **NUM)
        axs[0].plot(x_ref, u_ref, color=colour, **ANA)
        axs[1].plot(x, v_sim, color=colour, **NUM)
        axs[1].plot(x_ref, v_ref, color=colour, **ANA)

    for ax, lab, lo in ((axs[0], r"$E_r$", 1e-8), (axs[1], r"$a_R T^4$", 1e-5)):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylim(lo, 1.8)
        ax.set_xlabel(r"$x\ (\mathrm{cm})$")
        ax.set_ylabel(lab + r"$\ (\mathrm{erg\ cm^{-3}})$")
    axs[0].plot([], [], color="0.4", label="numerical", **NUM)
    axs[0].plot([], [], color="0.4", label="Su and Olson", **ANA)
    axs[0].legend(loc="lower left", fontsize=7.5)
    axs[0].set_title("Marshak wave", loc="left")
    return _finish(fig, outfile)


PLOTS = {
    "streaming": plot_streaming,
    "diffusion": plot_diffusion,
    "marshak": plot_marshak,
    "radshock": plot_radshock,
    "radtempeq": plot_radtempeq,
    # The linearised-temperature build solves the same problem, so it gets the
    # same figure; comparing the two PDFs shows where the schemes differ.
    "radtempeq_lin": plot_radtempeq,
    "advection": plot_advection,
}
