"""Minimal FLASH plotfile reader for the VETTAM unit tests.

Only numpy and h5py are required -- deliberately no yt, so the suite can run
on a bare CI worker. The unit tests all use a uniform (lrefine_max = 1) 1-D
grid, which keeps this reader to a few lines.
"""

import glob
import os

import h5py
import numpy as np


def plotfiles(directory, basename=None):
    """Every plotfile in *directory*, ordered by dump number."""
    pattern = f"{basename}hdf5_plt_cnt_*" if basename else "*_hdf5_plt_cnt_*"
    return sorted(glob.glob(os.path.join(directory, pattern)))


def read(path, fields=("erad",)):
    """Return (time, x, dx, {field: values}) for a 1-D FLASH plotfile.

    Cells are returned sorted by x, guard cells excluded (FLASH does not write
    them to plotfiles), and only from leaf blocks.
    """
    with h5py.File(path, "r") as h:
        leaf = h["node type"][:] == 1
        size = h["block size"][leaf][:, 0]
        centre = h["coordinates"][leaf][:, 0]
        ncell = h[fields[0]].shape[-1]

        dx = size / ncell
        offsets = (np.arange(ncell) + 0.5)[None, :]
        x = (centre - size / 2.0)[:, None] + offsets * dx[:, None]
        x = x.ravel()
        cell_dx = np.repeat(dx, ncell)

        data = {}
        for f in fields:
            data[f] = h[f][leaf][:, 0, 0, :].ravel()

        time = float(h["real scalars"][:]["value"][
            [n.decode().strip() for n in h["real scalars"][:]["name"]].index("time")
        ])

    order = np.argsort(x)
    return time, x[order], cell_dx[order], {k: v[order] for k, v in data.items()}


def l1_error(numerical, analytic, dx):
    """Volume-weighted L1 error, normalised by the analytic solution."""
    return float(np.sum(np.abs(numerical - analytic) * dx) / np.sum(analytic * dx))


def integral(values, dx):
    return float(np.sum(values * dx))
