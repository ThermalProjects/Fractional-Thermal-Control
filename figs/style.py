import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "vplant"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("FIGOUT", os.path.join(HERE, "out"))
os.makedirs(OUT, exist_ok=True)
W = 4.6          # text width in inches (117 mm)
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["TeX Gyre Termes", "Times"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "lines.linewidth": 1.0, "grid.linewidth": 0.4,
    "grid.alpha": 0.35, "axes.grid": True, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02, "figure.dpi": 150,
})
# grayscale-safe palette: black, dark grey, mid grey + linestyles
K, G1, G2, G3 = "0.0", "0.35", "0.55", "0.75"
C = ["#1f4e79", "#c0504d", "#4f8f3a", "#7f6000"]   # used sparingly, always with a linestyle too

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    plt.close(fig)
    print("saved", name)

def simlabel(ax, txt="simulated"):
    ax.text(0.99, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6, style="italic", color="0.4")
