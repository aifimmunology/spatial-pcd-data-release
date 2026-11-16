"""One-zip probe for the Space Ranger version extraction logic.

Downloads EXP00875_Visium_BM.zip from PCD_Spatial, reads only the small
housekeeping members from its central directory, and prints:

  1. The raw content of each housekeeping file (_version, _cmdline, _invocation)
     so you can read the version string with your own eyes.
  2. What the extraction logic parsed out and which file it came from.

Verification: compare (2) against what you can read in (1). If the extracted
version matches the string you can see in the raw text, the logic is correct.

Run from the repo root (requires an authenticated hisepy session):
    python visium/scripts/probe_spaceranger_version.py
"""

import glob
import os
import re
import zipfile

import hisepy as hp

STORE = "PCD_Spatial"
FILE_NAME = "AIFI-2024-06-10T23:23:30.477514464Z/EXP00875_Visium_BM.zip"
ZIP_BASE = "EXP00875_Visium_BM.zip"

RE_PATH_VER = re.compile(r"spaceranger[-_](\d+\.\d+\.\d+)", re.I)
RE_INV_VER = re.compile(r'pipelines_version\s*=\s*"([^"]+)"')
RE_WEB_VER = re.compile(r'"Pipeline Version","spaceranger-(\d+\.\d+\.\d+)"', re.I)


def read_member(zf, name):
    """Return decoded text of a zip member, or None if not present."""
    try:
        return zf.read(name).decode("utf-8", "replace")
    except KeyError:
        return None


def find_member_ending(names, suffix):
    """Return the shortest member path ending with `suffix`, or None."""
    hits = sorted((n for n in names if n.endswith(suffix)), key=len)
    return hits[0] if hits else None


def extract_version(zf):
    """Try four extraction strategies and return (version, source_file) or (None, None)."""
    names = zf.namelist()

    # Strategy 1: _version -- full content is the version string.
    target = find_member_ending(names, "_version")
    if target:
        txt = read_member(zf, target)
        if txt:
            ver = txt.strip()
            if re.match(r"^\d+\.\d+[\.\d]*$", ver):
                return ver, target

    # Strategy 2: _cmdline -- binary path contains "spaceranger-X.Y.Z".
    target = find_member_ending(names, "_cmdline")
    if target:
        txt = read_member(zf, target)
        if txt:
            m = RE_PATH_VER.search(txt)
            if m:
                return m.group(1), target

    # Strategy 3: _invocation -- pipelines_version = "X.Y.Z".
    target = find_member_ending(names, "_invocation")
    if target:
        txt = read_member(zf, target)
        if txt:
            m = RE_INV_VER.search(txt)
            if m:
                return m.group(1), target

    # Strategy 4: web_summary.html -- "Pipeline Version","spaceranger-X.Y.Z"
    # Fallback for outputs-only zips that have no housekeeping files.
    target = find_member_ending(names, "web_summary.html")
    if target:
        txt = read_member(zf, target)
        if txt:
            m = RE_WEB_VER.search(txt)
            if m:
                return m.group(1), target

    return None, None


def main():
    cands = glob.glob(os.path.expanduser(f"~/input/**/{ZIP_BASE}"), recursive=True)
    if cands:
        zp = cands[0]
        print(f"Already in FUSE mount, skipping download.")
    else:
        print(f"Downloading: {FILE_NAME}")
        hp.download_from_project_store(store_name=STORE, file_name=FILE_NAME)
        cands = glob.glob(os.path.expanduser(f"~/input/**/{ZIP_BASE}"), recursive=True)
        if not cands:
            raise FileNotFoundError(f"{ZIP_BASE} not found in FUSE mount after download")
        zp = cands[0]
    print(f"Found at: {zp}\n")

    with zipfile.ZipFile(zp) as zf:
        # Show every housekeeping member that the extraction logic looks at,
        # with its raw content -- so you can see exactly what was parsed.
        names = zf.namelist()
        for suffix in ("_version", "_cmdline", "_invocation", "web_summary.html"):
            target = find_member_ending(names, suffix)
            if target:
                txt = read_member(zf, target)
                preview = txt[:300].replace("\n", " ↵ ") if txt else "(empty)"
                print(f"--- {target} ---")
                print(preview)
                print()
            else:
                print(f"--- (no {suffix} member found) ---\n")

        ver, source = extract_version(zf)

    print("=" * 50)
    if ver:
        print(f"Extracted version : {ver}")
        print(f"Source file       : {source}")
        print()
        print("Verify: does the version above appear literally in the raw content printed above?")
        print("If yes, the extraction logic is correct.")
    else:
        print("Version NOT found. Review the member previews above to see what's there.")
        print("You may need to adjust the extraction logic in extract_spaceranger_version.py.")


if __name__ == "__main__":
    main()
