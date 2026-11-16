"""Generate raw spatial layout PDF — FFPE block photos vs. Xenium cell centroids.

One row per multi-tissue output: left = FFPE block photo (matched by block id),
right = cell centroids with sample-id labels.

Usage:
    python xenium/scripts/plot_raw_spatial.py
    python xenium/scripts/plot_raw_spatial.py --out figures/raw_spatial.pdf
    python xenium/scripts/plot_raw_spatial.py \\
        --xlsx xenium/data/xenium\\ data\\ tracking.xlsx \\
        --cells-dir xenium/data/multitissue_cells \\
        --blocks-zip xenium/data/Fred\\ Hutch\\ FFPE\\ blocks.zip \\
        --out figures/raw_spatial.pdf
"""

import argparse
import io
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

# ---------------------------------------------------------------------------
# Defaults (relative to repo root, i.e. the cwd when called as shown above)
# ---------------------------------------------------------------------------
REPO = Path("/home/workspace/PCD_data/xenium")
DEFAULT_XLSX = REPO / "data" / "xenium data tracking.xlsx"
DEFAULT_CELLS_DIR = REPO / "data" / "multitissue_cells"
DEFAULT_BLOCKS_ZIP = REPO / "data" / "Fred Hutch FFPE blocks.zip"
DEFAULT_OUT = REPO / "figures" / "raw_spatial.pdf"

# Manually-placed sample-id label positions (axes fraction).
# (0,0) = bottom-left, (1,1) = top-right of the centroid panel.
SAMPLE_LABEL_POS = {
    "TMA00238": {"TIS02600-001": (0.2, 0.25), "TIS02435-001": (0.8, 0.6)},
    "TMA00240": {"TIS06483-001": (0.25, 0.5), "TIS02619-001": (0.88, 0.85)},
    "TMA00244": {"TIS04965-001": (0.25, 0.25), "TIS02618-001": (0.75, 0.9)},
    "TMA00246": {"TIS02599-001": (0.18, 0.35), "TIS02429-001": (0.8, 0.525)},
    "TMA00142": {"TIS04966-001": (0.9, 0.48), "TIS04397-001": (0.2, 0.55)},
    "TMA00145": {"TIS02439-001": (0.25, 0.5), "TIS04436-001": (0.75, 0.5)},
    "TMA00161": {"TIS02608-001": (0.25, 0.3), "TIS02424-001": (0.75, 0.8)},
    "TMA00157": {"TIS04436-001": (0.25, 0.2), "TIS02439-001": (0.75, 0.7)},
    "TMA00159": {"TIS02609-001": (0.2, 0.15), "TIS02594-001": (0.75, 0.5)},
}


def load_groups(xlsx: Path) -> pd.DataFrame:
    main = pd.read_excel(xlsx, sheet_name="main")
    main.columns = [str(c).strip() for c in main.columns]
    sep = main[main["multiple samples"].astype(str).str.upper().str.startswith("Y")].copy()
    groups = (
        sep.groupby("file name")
        .agg(
            slide=("slide id", "first"),
            exp=("experiment id", "first"),
            section=("section id", "first"),
            block=("block id", "first"),
            samples=("potential sample ids (left, right)", "first"),
            subjects=("cohort", lambda s: sorted(s.dropna().astype(str).unique())),
        )
        .reset_index()
    )
    groups["label"] = groups["section"].astype(str)
    return groups


def load_cells(label: str, cells_dir: Path) -> pd.DataFrame | None:
    p = cells_dir / f"{label}.parquet"
    if not p.exists():
        return None
    return pd.read_parquet(p, columns=["cell_id", "x_centroid", "y_centroid"])


def load_block(block_id: str, blocks_zip: Path) -> Image.Image | None:
    if not blocks_zip.exists():
        return None
    with zipfile.ZipFile(blocks_zip) as z:
        name = next(
            (n for n in z.namelist() if Path(n).stem == block_id and n.lower().endswith(".jpg")),
            None,
        )
        if name is None:
            return None
        with z.open(name) as f:
            return Image.open(io.BytesIO(f.read())).copy()


def build_figure(groups: pd.DataFrame, cells_dir: Path, blocks_zip: Path) -> plt.Figure:
    labels_list = groups["label"].tolist()
    fig, axes = plt.subplots(len(labels_list), 2, figsize=(11, 4.5 * len(labels_list)))

    for (ax_img, ax_pts), lab in zip(axes, labels_list):
        row = groups[groups.label == lab].iloc[0]
        subs = " + ".join(row.subjects)
        lr = [s.strip() for s in str(row.samples).split(",")] if pd.notna(row.samples) else []

        # left: FFPE block photo
        img = load_block(row.block, blocks_zip)
        if img is None:
            ax_img.set_title(f"{row.block} — block photo MISSING")
        else:
            ax_img.imshow(img)
            ax_img.set_title(f"{row.block} — FFPE block")
            if len(lr) == 2:
                box = dict(boxstyle="round", fc="black", alpha=0.55, ec="none")
                ax_img.text(0.15, 0.75, f"L: {lr[0]}", transform=ax_img.transAxes,
                            ha="left", va="top", fontsize=8, color="white", bbox=box)
                ax_img.text(0.85, 0.75, f"R: {lr[1]}", transform=ax_img.transAxes,
                            ha="right", va="top", fontsize=8, color="white", bbox=box)
        ax_img.axis("off")

        # right: cell centroids
        df = load_cells(lab, cells_dir)
        if df is None:
            ax_pts.set_title(f"{lab} · {row.block} — MISSING (not downloaded)")
            ax_pts.axis("off")
            continue
        s = df.sample(60_000, random_state=0) if len(df) > 60_000 else df
        ax_pts.scatter(s.x_centroid, s.y_centroid, s=1, alpha=0.3, linewidths=0)
        ax_pts.set_aspect("equal")
        ax_pts.invert_yaxis()
        ax_pts.set_title(f"{lab} · {row.block}  ({subs})  n={len(df):,}")
        ax_pts.set_xlabel("x (µm)")
        ax_pts.set_ylabel("y (µm)")

        default_xy = [(0.25, 0.5), (0.75, 0.5)]
        lblbox = dict(boxstyle="round", fc="white", alpha=0.7, ec="0.4")
        for i, sid in enumerate(lr):
            x, y = SAMPLE_LABEL_POS.get(lab, {}).get(sid, default_xy[i] if i < 2 else (0.5, 0.5))
            ax_pts.text(x, y, sid, transform=ax_pts.transAxes, ha="center", va="center",
                        fontsize=8, color="black", bbox=lblbox)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.05)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX,
                        help="Path to 'xenium data tracking.xlsx'")
    parser.add_argument("--cells-dir", type=Path, default=DEFAULT_CELLS_DIR,
                        help="Directory containing <section>.parquet files")
    parser.add_argument("--blocks-zip", type=Path, default=DEFAULT_BLOCKS_ZIP,
                        help="Path to 'Fred Hutch FFPE blocks.zip'")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT,
                        help="Output PDF path")
    args = parser.parse_args()

    print(f"Loading tracking sheet: {args.xlsx}")
    groups = load_groups(args.xlsx)
    print(f"  {len(groups)} multi-tissue outputs\n")
    for _, r in groups.iterrows():
        subs = " + ".join(r.subjects)
        print(f"  {r.label:10s}  block {r.block}  {r.exp}  ->  {subs}")

    print(f"\nBuilding figure ...")
    fig = build_figure(groups, args.cells_dir, args.blocks_zip)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
