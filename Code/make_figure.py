#!/usr/bin/env python
"""Two-panel figure for the audit; reads ../Data/exhaustive.json."""
import os, json
from math import gcd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
recs = json.load(open(os.path.join(HERE, "..", "Data", "exhaustive.json")))

plt.rcParams.update({"font.size": 8, "font.family": "serif", "mathtext.fontset": "cm",
                     "axes.linewidth": 0.6, "pdf.fonttype": 42})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 3.05))

# ---- panel (a): exhaustive counts
rows = []
for r in recs:
    for t in r["tuples"]:
        rows.append((t["normal_per_normal_tuple"], t["gcd"], t["pn_min"], t["pn_max"],
                     t["max_normal_on_nonnormal_tuple"]))
rows.sort()
xs = list(range(1, len(rows) + 1))
C1, C2, C3, C4 = "#1f4e79", "#c0504d", "#7f7f7f", "#2e7d32"
ax.plot(xs, [r[0] for r in rows], "_", color=C3, ms=9, mew=1.2,
        label=r"predicted $\Phi(x^m-1)/\Phi(g)$", zorder=1)
for g1, mk, lab in ((True, "o", r"normal $\mathbf{a}$, $\gcd(d_1,d_2)=1$"),
                    (False, "s", r"normal $\mathbf{a}$, $\gcd(d_1,d_2)>1$")):
    sel = [(x, r) for x, r in zip(xs, rows) if (r[1] == 1) == g1]
    ax.plot([x for x, _ in sel], [r[0] for _, r in sel], mk, color=C1, ms=3.2,
            mfc=C1 if g1 else "white", mew=0.8, label=lab, zorder=3)
ax.vlines(xs, [r[2] for r in rows], [r[3] for r in rows], color=C4, lw=2.2,
          label="primitive normal (min\u2013max)", zorder=2)
ax.plot(xs, [r[4] for r in rows], "x", color=C2, ms=3.5, mew=0.9,
        label=r"non-normal admissible $\mathbf{a}$", zorder=4)
ax.set_yscale("symlog", linthresh=1, linscale=0.5)
ax.set_ylim(-0.6, 3e6)
ax.set_yticks([0, 1, 10, 100, 1e3, 1e4, 1e5, 1e6])
ax.set_xlim(0, len(rows) + 1)
ax.set_xlabel(r"test case $(q,m,\mathbf{d})$, ordered by predicted count")
ax.set_ylabel(r"elements with $\mathrm{Tr}_{m/d_i}(\alpha)=a_i$")
ax.text(0.03, 0.93, "(a)", transform=ax.transAxes, fontweight="bold")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False,
          fontsize=6.6, handletextpad=0.4, columnspacing=1.0, borderaxespad=0)

# ---- panel (b): exponents for k = 2, m = d1 d2
cols = {3: C2, 5: C4, 7: C1}
for d1 in (3, 5, 7):
    d2s = [d2 for d2 in range(d1 + 1, 27) if gcd(d1, d2) == 1]
    bx.plot(d2s, [d1 * d2 / 2 - 2 * (d1 + d2) + 1 for d2 in d2s], "--", marker="o", ms=2.2,
            lw=0.8, color=cols[d1], mfc="white", mew=0.6,
            label=rf"$d_1={d1}$: printed, $m/2-\lambda-D$")
    bx.plot(d2s, [d1 * d2 / 2 - (d1 + d2) + 1 for d2 in d2s], "-", marker="o", ms=2.2, lw=0.9,
            color=cols[d1], label=rf"$d_1={d1}$: corrected, $m/2-\lambda$")
bx.axhline(0, color="k", lw=0.6)
bx.set_ylim(-14, 30)
bx.set_xlim(3.5, 26.5)
bx.set_xlabel(r"$d_2$ (with $\gcd(d_1,d_2)=1$, $m=d_1d_2$)")
bx.set_ylabel(r"exponent of $q$ in the sufficient condition")
bx.text(0.03, 0.93, "(b)", transform=bx.transAxes, fontweight="bold")
_h, _l = bx.get_legend_handles_labels()
_o = [0, 2, 4, 1, 3, 5]
bx.legend([_h[i] for i in _o], [_l[i] for i in _o], loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False,
          fontsize=6.6, handletextpad=0.4, columnspacing=1.0, borderaxespad=0)
for a in (ax, bx):
    a.tick_params(width=0.6, length=2.5)
fig.tight_layout(w_pad=1.6)
fig.savefig(os.path.join(HERE, "..", "fig_intermediate_traces_audit.pdf"))
print("cases:", len(rows))
