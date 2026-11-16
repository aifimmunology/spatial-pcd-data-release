#!/usr/bin/env python
"""Generate a content tree for each zip in the PCD_Spatial project store
WITHOUT keeping the zips on disk.

The full zips are large (tens of GB each), so instead of holding them all we
process them ONE AT A TIME:

    download zip -> read its internal structure -> write a .md tree file
    -> delete the zip -> move on to the next.

Peak disk usage is therefore a single zip at any moment, never all of them.
hisepy only offers a full download (no byte-range read), so each zip must be
fully fetched once before we can list its contents.

If a zip is already present under TARGET_DIR (e.g. left over from
visium_data_download.py), it is read in place instead of re-downloaded, and
then deleted after its tree is written.

Output: one flat Markdown file per zip under TREE_DIR, named after the zip's
basename (e.g. 'visium-zip-trees/EXP00875_Visium_BM.zip.md'). A few zips share
the same basename under different AIFI ids; only those are disambiguated by
prefixing the AIFI id so nothing gets overwritten.

Usage:
    python visium_zip_tree.py
"""

import glob
import os
import zipfile
from collections import Counter

import hisepy as hp

# Project store holding the data (the store is named "PCD_Spatial").
STORE_NAME = "PCD_Spatial"

# Where downloaded zips are staged/kept (matches visium_data_download.py).
TARGET_DIR = "/home/workspace/PCD_data/visium/data/visium-data"

# Where the generated tree files are written.
TREE_DIR = "/home/workspace/PCD_data/visium/data/visium-zip-trees"

# Where hise stages downloaded files before we can read them.
INPUT_ROOT = os.path.expanduser("~/input")

# Exact file names in the project store, as listed by
# list_files_in_project_store().
FILES = [
    "AIFI-PSP-2025-11-11T18:48:44.995148662Z/TIS05393-001-010-EXP-00971-SB005.zip",
    "AIFI-2024-07-03T19:50:28.897069523Z/SB008_EXP-01092_TIS06377-001_corrected_sr_outs.zip",
    "AIFI-2024-07-03T19:51:46.907032375Z/SB008_EXP-01092_TIS06392-002_corrected_sr_outs.zip",
    "AIFI-2024-07-03T19:52:34.884583329Z/SB008_EXP-01095_TIS05681-001_sr_outs.zip",
    "AIFI-2024-07-03T19:53:08.243222988Z/SB008_EXP-01095_TIS05683-001_sr_outs.zip",
    "AIFI-2024-07-03T20:01:18.179156501Z/SB008_EXP-01106_TIS05684-001-003.zip",
    "AIFI-2024-07-03T20:10:25.871592784Z/SB008_EXP-01106_TIS06402-001-004.zip",
    "AIFI-2024-06-10T23:11:33.601005705Z/TIS05672_EXP00950_SB5.zip",
    "AIFI-2024-06-10T23:12:12.363347736Z/TIS05671_EXP00949_SB5.zip",
    "AIFI-2024-06-10T23:18:19.262091389Z/TIS04789_EXP00948_SB5.zip",
    "AIFI-2024-06-10T23:19:27.648275352Z/TIS04788_EXP00947_SB5.zip",
    "AIFI-2024-03-14T17:52:07.667273266Z/TIS05398-001-007-EXP-00971-SB005.zip",
    "AIFI-2024-03-14T17:51:49.757643408Z/EXP-00971 AACR 3+4-20240229T010204Z-001.zip",
    "AIFI-2024-03-14T17:52:01.813837457Z/TIS05394-001-010-EXP-00971-SB005.zip",
    "AIFI-2024-06-10T23:08:48.250692659Z/TIS05393_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:09:28.900356701Z/TIS05395_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:12:49.667054782Z/TIS05398_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:15:03.199377095Z/TIS05395_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:16:03.199656687Z/TIS05394_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:17:13.523763688Z/TIS05393_EXP00971_AACR_SB5.zip",
    "AIFI-2024-03-14T17:51:55.870095031Z/SB005_EXP00963_rerun_outs.zip",
    "AIFI-2024-06-10T23:13:27.824091334Z/TIS05397_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:14:13.664453937Z/TIS05396_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:17:45.568843181Z/TIS05392_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:18:52.069714292Z/TIS04788_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:20:01.122865514Z/EXP00937_BM.zip",
    "AIFI-2024-06-10T23:22:06.740736449Z/EXP00898_BM.zip",
    "AIFI-2024-06-10T23:22:51.980946983Z/EXP00897_BM.zip",
    "AIFI-2024-06-10T23:21:35.580329354Z/EXP00934_BM.zip",
    "AIFI-2024-06-10T23:20:26.232014326Z/EXP00936_BM.zip",
    "AIFI-2024-06-10T23:21:01.685829928Z/EXP00935_BM.zip",
    "AIFI-2024-06-10T23:23:30.477514464Z/EXP00875_Visium_BM.zip",
    "AIFI-2024-07-03T20:20:20.619605327Z/SB008_EXP-01106_TMA00144.zip",
    "AIFI-2024-07-03T20:30:34.747740525Z/SB008_EXP-01106_TMA00147_full.zip",
]


# Basenames that occur more than once in FILES need disambiguating on output.
_BASENAME_COUNTS = Counter(os.path.basename(f) for f in FILES)


def tree_path(file_name):
    """Return the flat output path for file_name's tree.

    Normally 'TREE_DIR/<basename>.md'. If the basename is shared by more than
    one zip, prefix it with the AIFI id so the files don't collide.
    """
    base = os.path.basename(file_name)
    if _BASENAME_COUNTS[base] > 1:
        prefix = file_name.split("/")[0]
        return os.path.join(TREE_DIR, f"{prefix}__{base}.md")
    return os.path.join(TREE_DIR, base + ".md")


def locate_local(file_name):
    """Return a path to an already-present copy of file_name, or None.

    Checks TARGET_DIR first, then anywhere hise may have staged it.
    """
    candidate = os.path.join(TARGET_DIR, file_name)
    if os.path.isfile(candidate):
        return candidate
    matches = glob.glob(os.path.join(INPUT_ROOT, "**", file_name), recursive=True)
    if matches:
        return max(matches, key=os.path.getmtime)
    return None


def human_size(num_bytes):
    """Format a byte count as a short human-readable string."""
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.1f}{unit}" if unit != "B" else f"{int(size)}B"
        size /= 1024


def build_tree(infolist):
    """Build an ASCII tree (list of lines) from a zip's ZipInfo list.

    Directories are shown with their entries nested underneath; each file is
    annotated with its uncompressed size.
    """
    # Nested dict: name -> child dict; files store a special key for size.
    root = {}
    for info in infolist:
        parts = [p for p in info.filename.split("/") if p]
        if not parts:
            continue
        node = root
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        last = parts[-1]
        if info.is_dir():
            node.setdefault(last, {})
        else:
            child = node.setdefault(last, {})
            child["__size__"] = info.file_size

    lines = []

    def walk(node, prefix):
        # Sort: directories first, then files, alphabetically.
        entries = sorted(
            (k for k in node if k != "__size__"),
            key=lambda k: ("__size__" in node[k], k.lower()),
        )
        for i, name in enumerate(entries):
            last = i == len(entries) - 1
            connector = "└── " if last else "├── "
            child = node[name]
            is_file = "__size__" in child
            label = name
            if is_file:
                label = f"{name}  ({human_size(child['__size__'])})"
            else:
                label = f"{name}/"
            lines.append(f"{prefix}{connector}{label}")
            if not is_file:
                extension = "    " if last else "│   "
                walk(child, prefix + extension)

    walk(root, "")
    return lines


def write_tree_file(file_name, infolist):
    """Write a Markdown tree file for file_name and return its path."""
    dst = tree_path(file_name)

    n_files = sum(1 for i in infolist if not i.is_dir())
    n_dirs = sum(1 for i in infolist if i.is_dir())
    total = sum(i.file_size for i in infolist if not i.is_dir())
    tree_lines = build_tree(infolist)

    with open(dst, "w") as fh:
        fh.write(f"# {os.path.basename(file_name)}\n\n")
        fh.write(f"- Project store: `{STORE_NAME}`\n")
        fh.write(f"- Store path: `{file_name}`\n")
        fh.write(f"- Entries: {n_files} files, {n_dirs} directories\n")
        fh.write(f"- Uncompressed total: {human_size(total)}\n\n")
        fh.write("```\n")
        fh.write(os.path.basename(file_name) + "\n")
        fh.write("\n".join(tree_lines))
        fh.write("\n```\n")
    return dst


def process(file_name):
    """Ensure file_name is available locally, write its tree, delete the zip.

    Returns True on success.
    """
    local = locate_local(file_name)
    downloaded = False

    # A pre-existing copy may be a truncated/incomplete download (valid PK
    # header but no central directory). Discard it and re-download if so.
    if local is not None and not zipfile.is_zipfile(local):
        print(f"  existing copy is incomplete/corrupt, removing: {local}")
        try:
            os.remove(local)
        except OSError as err:
            print(f"  WARNING: could not delete corrupt copy: {err}")
        local = None

    if local is None:
        print("  downloading ...")
        ok = hp.download_from_project_store(
            store_name=STORE_NAME,
            file_name=file_name,
        )
        if not ok:
            print("  download returned False")
            return False
        local = locate_local(file_name)
        downloaded = True
        if local is None:
            print("  ERROR: download succeeded but staged file not found")
            return False
    else:
        print(f"  using existing copy: {local}")

    try:
        with zipfile.ZipFile(local) as zf:
            infolist = zf.infolist()
        dst = write_tree_file(file_name, infolist)
        print(f"  -> {dst}")
        return True
    except zipfile.BadZipFile as err:
        print(f"  ERROR: not a valid zip ({err})")
        return False
    finally:
        # Delete the zip to reclaim space before the next one, whether it was
        # just downloaded or already present.
        try:
            os.remove(local)
            print(f"  deleted {local} ({'downloaded' if downloaded else 'pre-existing'})")
        except OSError as err:
            print(f"  WARNING: could not delete {local}: {err}")


def main():
    os.makedirs(TREE_DIR, exist_ok=True)
    succeeded, failed = [], []
    for i, file_name in enumerate(FILES, start=1):
        print(f"[{i}/{len(FILES)}] {file_name}")
        try:
            if process(file_name):
                succeeded.append(file_name)
            else:
                failed.append(file_name)
        except Exception as err:  # keep going if one file fails
            print(f"  ERROR: {err}")
            failed.append(file_name)

    print(f"\nDone. {len(succeeded)} tree files in {TREE_DIR}, {len(failed)} failed.")
    if failed:
        print("Failed files:")
        for name in failed:
            print(f"  {name}")


if __name__ == "__main__":
    main()
