# Reproduz todas as figuras e numeros do artigo.
# Uso:  python run_experiments.py  ->  gera ../figures/ e ../results.json

from __future__ import annotations
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from seir_ode import solve_seir, final_size_and_peak, MEASLES
from seir_ca import SEIR_CA
from similarity_dimension import (box_count, infection_dimension_over_time,
                                  koch_curve, koch_similarity_dimension,
                                  rasterize_curve)

HERE = os.path.dirname(__file__)
FIGDIR = os.path.abspath(os.path.join(HERE, "..", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold",
    "legend.frameon": False, "font.family": "DejaVu Sans",
})

# cores para S, E, I, R
CMAP = ListedColormap(["#dfe7ef", "#f4b942", "#d1495b", "#4a5859"])
COL = {"S": "#8aa1b4", "E": "#f4b942", "I": "#d1495b", "R": "#4a5859"}

results = {}


def fig1_ode_vs_ca():
    t, Y = solve_seir(N=10_000, I0=10, days=120, dt=0.1)
    ar, peak, day = final_size_and_peak(t, Y)
    results["ode"] = {"attack_rate": ar, "peak_prev": peak, "peak_day": day,
                      **{k: MEASLES[k] for k in ("R0", "beta", "sigma", "gamma")}}

    tau = 0.30
    ca = SEIR_CA(L=200, tau=tau, seed=3)
    ts, _ = ca.run(steps=800)
    ca_peak = ts[:, 2].max()
    results["ca_headline"] = {"tau": tau, "R0_proxy": ca.local_R0_proxy(),
                              "attack_rate": float(ts[-1, 3]),
                              "peak_prev": float(ca_peak),
                              "peak_day": int(ts[:, 2].argmax()),
                              "duration": len(ts) - 1}

    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    for j, lab in enumerate(["S", "E", "I", "R"]):
        ax[0].plot(t, Y[:, j], color=COL[lab], lw=2, label=lab)
    ax[0].set(title="(a) SEIR de campo médio (EDO)", xlabel="tempo (dias)",
              ylabel="fração da população", xlim=(0, 120))
    ax[0].legend(ncol=4, loc="upper right", fontsize=8)

    steps = np.arange(len(ts))
    for j, lab in enumerate(["S", "E", "I", "R"]):
        ax[1].plot(steps, ts[:, j], color=COL[lab], lw=2, label=lab)
    ax[1].set(title="(b) SEIR espacial (autômato celular)",
              xlabel="tempo (passos = dias)", ylabel="fração da população")
    ax[1].legend(ncol=4, loc="center right", fontsize=8)

    ax[2].plot(t, Y[:, 2], color="black", lw=2, label=f"EDO  (pico {peak*100:.1f}%)")
    ax[2].plot(steps, ts[:, 2], color=COL["I"], lw=2,
               label=f"AC  (pico {ca_peak*100:.1f}%)")
    ax[2].set(title="(c) Prevalência de infectados I(t)", xlabel="tempo (dias)",
              ylabel="fração infecciosa", xlim=(0, 800))
    ax[2].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig1_ode_vs_ca.png"), bbox_inches="tight")
    plt.close(fig)
    print("fig1 ok:", results["ca_headline"])


def fig2_ca_snapshots():
    frames_at = [5, 60, 150, 300, 500, 700]
    ca = SEIR_CA(L=200, tau=0.30, seed=3)
    ts, frames = ca.run(steps=max(frames_at), record_frames=frames_at)

    fig, axes = plt.subplots(1, len(frames_at), figsize=(13, 2.6))
    for ax, k in zip(axes, frames_at):
        g = frames.get(k, ca.grid)
        ax.imshow(g, cmap=CMAP, vmin=0, vmax=3, interpolation="nearest")
        ax.set_title(f"passo {k}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    handles = [Patch(color=CMAP.colors[i], label=l)
               for i, l in enumerate(["S", "E", "I", "R"])]
    fig.legend(handles=handles, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.08), fontsize=9)
    fig.suptitle("Propagação espacial do autômato celular SEIR "
                 "(núcleo infeccioso central único)", y=1.02, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_ca_snapshots.png"), bbox_inches="tight")
    plt.close(fig)
    print("fig2 ok")


def fig3_dimension():
    x, y = koch_curve(order=6)
    img = rasterize_curve(x, y, grid=1024)
    d_box, sizes, counts = box_count(img)
    d_exact = koch_similarity_dimension()
    results["koch"] = {"D_similarity_exact": d_exact, "D_boxcount": d_box}

    record = list(range(10, 720, 20))
    ca = SEIR_CA(L=200, tau=0.30, seed=3)
    _, frames = ca.run(steps=max(record), record_frames=record)
    dim_t = infection_dimension_over_time(frames, states=(1, 2))
    ks = sorted(dim_t)
    dims = [dim_t[k] for k in ks]
    results["front_dimension"] = {"mean": float(np.nanmean(dims)),
                                  "min": float(np.nanmin(dims)),
                                  "max": float(np.nanmax(dims))}

    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.2))
    ax[0].plot(x, y, color="#2b3a67", lw=0.8)
    ax[0].set_title("(a) Curva de Koch"); ax[0].axis("equal"); ax[0].axis("off")
    ax[0].text(0.5, -0.15, f"$D_s=\\log4/\\log3={d_exact:.4f}$",
               ha="center", transform=ax[0].transAxes, fontsize=9)

    m = counts > 0
    xx, yy = np.log(1.0 / sizes[m]), np.log(counts[m])
    slope, intercept = np.polyfit(xx, yy, 1)
    ax[1].plot(xx, yy, "o", color="#2b3a67", ms=5)
    ax[1].plot(xx, slope * xx + intercept, "-", color="#d1495b",
               label=f"inclinação $D_{{cx}}={slope:.3f}$")
    ax[1].set(title="(b) Contagem de caixas da curva de Koch",
              xlabel=r"$\log(1/\varepsilon)$", ylabel=r"$\log N(\varepsilon)$")
    ax[1].legend(fontsize=9)

    ax[2].plot(ks, dims, "o-", color=COL["I"], ms=4)
    ax[2].axhline(2.0, ls=":", color="gray"); ax[2].axhline(1.0, ls=":", color="gray")
    ax[2].set(title="(c) Dimensão da frente ativa (E+I)",
              xlabel="tempo (passos)", ylabel=r"$D_{cx}$", ylim=(0.8, 2.05))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig3_dimension.png"), bbox_inches="tight")
    plt.close(fig)
    print("fig3 ok:", results["koch"], results["front_dimension"])


def fig4_phase_transition():
    taus = np.array([0.03, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12,
                     0.15, 0.20, 0.30, 0.50, 0.70])
    seeds = range(5)
    attack_mean, attack_sd, clust_dim = [], [], []
    for tau in taus:
        ars, dims = [], []
        for s in seeds:
            ca = SEIR_CA(L=150, tau=float(tau), seed=s)
            ts, _ = ca.run(steps=1200)
            ars.append(ts[-1, 3])
            if ts[-1, 3] > 0.02:                       # apenas execuções que decolaram
                d, _, _ = box_count(np.isin(ca.grid, (1, 2, 3)))
                dims.append(d)
        attack_mean.append(np.mean(ars))
        attack_sd.append(np.std(ars))
        clust_dim.append(np.nan if not dims else np.mean(dims))
    results["phase"] = {"tau": taus.tolist(),
                        "attack_mean": [float(a) for a in attack_mean],
                        "cluster_dim": [None if np.isnan(d) else float(d)
                                        for d in clust_dim]}

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    ax[0].errorbar(taus, np.array(attack_mean) * 100, yerr=np.array(attack_sd) * 100,
                   fmt="o-", color="#2b3a67", capsize=3)
    ax[0].axvspan(0.06, 0.08, color="#d1495b", alpha=0.12)
    ax[0].set(title="(a) Limiar de invasão",
              xlabel=r"transmissibilidade por contato $\tau$",
              ylabel="taxa de ataque final (%)")
    ax[0].annotate("limiar\nepidêmico", xy=(0.07, 40), xytext=(0.15, 30),
                   fontsize=8, arrowprops=dict(arrowstyle="->", color="gray"))

    dvals = np.array([np.nan if d is None else d for d in results["phase"]["cluster_dim"]])
    mm = ~np.isnan(dvals)
    ax[1].plot(taus[mm], dvals[mm], "s-", color=COL["I"])
    ax[1].axhline(2.0, ls=":", color="gray")
    ax[1].axhline(1.896, ls="--", color="#2b3a67", lw=1, label="percolação 2D $\\approx1{,}90$")
    ax[1].set(title="(b) Dimensão do cluster infectado",
              xlabel=r"transmissibilidade por contato $\tau$",
              ylabel=r"$D_{cx}$", ylim=(1.6, 2.05))
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig4_phase_transition.png"), bbox_inches="tight")
    plt.close(fig)
    print("fig4 ok")


if __name__ == "__main__":
    fig1_ode_vs_ca()
    fig2_ca_snapshots()
    fig3_dimension()
    fig4_phase_transition()
    with open(os.path.join(HERE, "..", "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("\nFiguras geradas em", FIGDIR)
