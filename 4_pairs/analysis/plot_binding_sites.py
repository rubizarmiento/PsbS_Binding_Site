"""
Plot a top view (XY plane) of the PSII supercomplex with selected proteins
colored and labeled by name.

All PSII beads are drawn as a light grey background. The highlighted chains
(CP24, CP26, CP29, CP43, the S- and M-LHCII trimers and the small subunits
PsbE, PsbF, PsbJ, PsbK, PsbP and PsbZ) are drawn on top with the colors from
``chain_colors.yaml``, and the protein name from ``chain_labels.yaml`` is
written at the center of each chain. The figure is saved as a transparent PNG.

Text uses Arial when it is installed (it is not shipped with matplotlib, but
macOS provides it). Otherwise the closest metric-compatible font available is
used: Helvetica, Liberation Sans, Arimo, and finally DejaVu Sans.

Arguments
---------
-f, --file : str
    Input PSII PDB file (default: 5_psii/base_dir/rotated.pdb).
-colors, --colors : str
    YAML file mapping chainID -> color (default: definitions_yaml/chain_colors.yaml).
-labels, --labels : str
    YAML file mapping chainID -> protein name (default: definitions_yaml/chain_labels.yaml).
-o, --output : str
    Output figure (default: 4_pairs/analysis/figures/binding_sites/binding_sites.png).
-fontsize, --fontsize : float
    Font size of the protein labels in points (default: 9).

Defaults are relative paths, so run from the repository root:

    python 4_pairs/analysis/plot_binding_sites.py
    python 4_pairs/analysis/plot_binding_sites.py -fontsize 12 -o figure.png
"""

import argparse
import os
import warnings

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import MDAnalysis as mda
import yaml

# Ignores warnings from MDAnalysis
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

# Arial first; the rest are fallbacks with (nearly) the same glyph metrics
FONT_FAMILY = ["Arial", "Helvetica", "Liberation Sans", "Arimo", "DejaVu Sans"]

# Proteins that get a name written at their center
LABELED_PROTEINS = ["LHCBM", "LHCB3", "CP24", "CP26", "CP29", "CP43",
                    "PsbE", "PsbF", "PsbJ", "PsbK", "PsbP", "PsbZ"]

# Chains to highlight (both monomers); colors come from chain_colors.yaml
# Drawn in this order: PsbP (lumenal, covers its neighbours in top view) goes first.
CHAINS = [
    "P", "p",                                         # PsbP (grey scale)
    "4", "s", "r", "c", "8", "S", "R", "C",           # CP24, CP26, CP29, CP43
    "G", "N", "Y", "g", "n", "y",                     # S-LHCII trimers
    "1", "2", "3", "5", "6", "7",                     # M-LHCII trimers
    "E", "e", "F", "f", "J", "j",                     # PsbE, PsbF, PsbJ (grey scale)
    "K", "k", "Z", "z",                               # PsbK, PsbZ (grey scale)
]

# Label offsets in Angstrom (x, y). PsbE and PsbF (cyt b559) share the same
# center in top view, so their labels are stacked above/below each other.
LABEL_OFFSETS = {"PsbE": (0, 4), "PsbF": (0, -4)}


def parse_args():
    """
    Parse the command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments: ``file``, ``colors``, ``labels``, ``output`` and
        ``fontsize``. See the module docstring for their meaning and defaults.
    """
    parser = argparse.ArgumentParser(
        description='Plot a top view of PSII with the highlighted chains '
                    'colored and labeled by protein name')
    parser.add_argument('-f', '--file', type=str,
                        default='5_psii/base_dir/rotated.pdb',
                        help='Input PSII PDB file (default: %(default)s)')
    parser.add_argument('-colors', '--colors', type=str,
                        default='definitions_yaml/chain_colors.yaml',
                        help='YAML file with chainID -> color (default: %(default)s)')
    parser.add_argument('-labels', '--labels', type=str,
                        default='definitions_yaml/chain_labels.yaml',
                        help='YAML file with chainID -> protein name (default: %(default)s)')
    parser.add_argument('-o', '--output', type=str,
                        default='4_pairs/analysis/figures/binding_sites/binding_sites.png',
                        help='Output figure (default: %(default)s)')
    parser.add_argument('-fontsize', '--fontsize', type=float, default=9,
                        help='Font size of the protein labels in points (default: %(default)s)')
    return parser.parse_args()


def load_chain_colors(path):
    """
    Load the chainID -> color mapping.

    Parameters
    ----------
    path : str
        YAML file with entries such as ``"s": "#fdb462"``.

    Returns
    -------
    dict of str to str
        Color (any matplotlib color string) for each chainID.
    """
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_chain_labels(path):
    """
    Load the chainID -> protein name mapping.

    Names starting with ``PSB`` are rewritten in the usual notation,
    e.g. ``PSBE`` -> ``PsbE``. Other names (CP24, LHCBM, ...) are kept as is.

    Parameters
    ----------
    path : str
        YAML file with entries such as ``"e": PSBE``.

    Returns
    -------
    dict of str to str
        Protein name for each chainID.
    """
    with open(path, "r", encoding="utf-8") as f:
        labels = yaml.safe_load(f)
    return {chain: (f"Psb{name[3:]}" if name.startswith("PSB") else name)
            for chain, name in labels.items()}


def text_color_for(background):
    """
    Choose black or white text for readability on a given background.

    Parameters
    ----------
    background : str
        Background color (any matplotlib color string).

    Returns
    -------
    str
        ``"white"`` for dark backgrounds, ``"black"`` otherwise.
    """
    r, g, b = mcolors.to_rgb(background)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if luminance < 0.5 else "black"


def plot_background(ax, universe):
    """
    Draw every bead of the structure as a light grey background.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    universe : MDAnalysis.Universe
        PSII structure.
    """
    positions = universe.atoms.positions
    ax.scatter(positions[:, 0], positions[:, 1], color="lightgrey", s=70)


def plot_chains(ax, universe, chains, colors):
    """
    Draw the selected chains on top of the background, in the given order.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    universe : MDAnalysis.Universe
        PSII structure.
    chains : list of str
        ChainIDs to draw; later chains are drawn on top of earlier ones.
    colors : dict of str to str
        Color for each chainID.
    """
    for chain in chains:
        positions = universe.select_atoms(f"chainID {chain}").positions
        ax.scatter(positions[:, 0], positions[:, 1], color=colors[chain], s=100)


def plot_labels(ax, universe, chains, colors, labels, labeled_proteins,  # pylint: disable=too-many-arguments,too-many-positional-arguments
                offsets, fontsize):
    """
    Write the protein name at the XY center of each selected chain.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    universe : MDAnalysis.Universe
        PSII structure.
    chains : list of str
        Candidate chainIDs to label.
    colors : dict of str to str
        Color for each chainID, used to pick black or white text.
    labels : dict of str to str
        Protein name for each chainID.
    labeled_proteins : list of str
        Only chains whose protein name is in this list are labeled.
    offsets : dict of str to tuple of float
        Optional (x, y) shift in Angstrom of the label, per protein name.
    fontsize : float
        Font size in points.
    """
    for chain in chains:
        name = labels[chain]
        if name not in labeled_proteins:
            continue
        center = universe.select_atoms(f"chainID {chain}").positions[:, :2].mean(axis=0)
        center = center + offsets.get(name, (0, 0))
        ax.text(center[0], center[1], name, ha="center", va="center",
                fontsize=fontsize, fontweight="bold",
                color=text_color_for(colors[chain]))


def clean_axes(ax):
    """
    Use an equal aspect ratio and remove ticks and borders.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to clean.
    """
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def save_figure(fig, output):
    """
    Save the figure as a transparent PNG, creating the folder if needed.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    output : str
        Output file path.
    """
    odir = os.path.dirname(output)
    if odir:
        os.makedirs(odir, exist_ok=True)
    fig.savefig(output, transparent=True, bbox_inches="tight", dpi=600)


def main():
    """
    Build the PSII top-view figure from the command line arguments and save it.
    """
    args = parse_args()

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = FONT_FAMILY

    colors = load_chain_colors(args.colors)
    labels = load_chain_labels(args.labels)
    universe = mda.Universe(args.file)

    fig, ax = plt.subplots(figsize=(10, 10))
    plot_background(ax, universe)
    plot_chains(ax, universe, CHAINS, colors)
    # Labels after all chains so they stay on top
    plot_labels(ax, universe, CHAINS, colors, labels, LABELED_PROTEINS,
                LABEL_OFFSETS, args.fontsize)
    clean_axes(ax)

    save_figure(fig, args.output)

    print(f"Figure saved to {args.output}")

if __name__ == "__main__":
    main()
