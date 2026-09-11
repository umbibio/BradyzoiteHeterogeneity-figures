"""Three-dimensional UMAP scatter, rendered exactly as in the published figures.

The panels are not 2-D scatters of a 2-D embedding: they are matplotlib 3-D axes
holding a 3-D UMAP, viewed from ``elev=60, azim=0``. Matplotlib projects the
*unit-cube-normalised* coordinates, so the aspect of each panel depends on the
range of the particular subset being drawn. Projecting to 2-D by hand and
scattering that gives a visibly different picture for some subsets, which is why
this keeps the original 3-D rendering path.

Vendored from ``notebooks/movies.py`` of the analysis repository (function
``plot_3d_preview``), unchanged apart from documentation and typing.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

# Viewing angle used for every published 3-D panel.
ELEVATION = 60
AZIMUTH = 0


def plot_3d_preview(
    points: NDArray,
    colors: Sequence | None = None,
    point_size: float = 5,
    cmap: str | None = None,
    alpha: float = 0.7,
    background_color: str = "white",
    dpi: int = 100,
    padding_factor: float = 0.1,
    elevation: float = 6,
    angle: float = 0,
    vmin: float | None = None,
    vmax: float | None = None,
    remove_panels: bool = False,
    show_colorscale: bool = False,
    tick_step: float = 1.0,
):
    fig = plt.figure(figsize=(12, 9), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")

    if colors is None:
        colors = points[:, 2]
        cmap = "plasma"

    ax.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2],
        c=colors,
        cmap=cmap,
        s=point_size,
        alpha=alpha,
        vmin=vmin,
        vmax=vmax,
    )

    ax.set_facecolor(background_color)
    fig.patch.set_facecolor(background_color)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    # Tight bounds with a small padding, computed from this subset alone.
    for axis, setter in enumerate(("set_xlim", "set_ylim", "set_zlim")):
        low, high = points[:, axis].min(), points[:, axis].max()
        pad = (high - low) * padding_factor
        getattr(ax, setter)(low - pad, high + pad)

    ax.view_init(elev=elevation, azim=angle)

    if show_colorscale:
        mappable = cm.ScalarMappable(norm=mcolors.Normalize(vmin=vmin, vmax=vmax), cmap=plt.get_cmap(cmap))
        mappable.set_array([])
        bar = plt.colorbar(mappable, ax=ax, shrink=0.75)
        bar.set_label("Normalized Expression")
        bar.set_ticks(ticks=np.arange(0, vmax, tick_step))

    if remove_panels:
        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pane.fill = False
            pane.set_edgecolor("none")
        ax.grid(False)
        ax.set_axis_off()

    return ax


def matplotlib_projection(elevation_deg: float, azimuth_deg: float) -> NDArray:
    """The 3x2 projection matplotlib's 3-D axes applies for ``view_init``.

    Only needed to reproduce the ``2d_projection`` embedding stored alongside the
    3-D one; the panels themselves are rendered by :func:`plot_3d_preview`.
    """
    elev = np.radians(elevation_deg)
    azim = np.radians(azimuth_deg - 90)
    cos_e, sin_e = np.cos(elev), np.sin(elev)
    cos_a, sin_a = np.cos(azim), np.sin(azim)
    return np.array([[cos_a, -sin_a, 0], [sin_e * sin_a, sin_e * cos_a, cos_e]])
