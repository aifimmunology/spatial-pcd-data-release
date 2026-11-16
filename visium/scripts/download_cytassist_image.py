"""
Download a single zip from the HISE project store `PCD_Spatial` and extract
ONLY the CytAssist image (`cytassist_image.tiff`) from it (without unpacking
the whole archive).

Target zip (exact key/name):
    AIFI-2024-06-10T23:08:48.250692659Z/TIS05393_EXP00971_AACR_SB5.zip

The zip contains several copies of `cytassist_image.tiff` (intermediate
SPATIAL_REPORTER pipeline copies plus the canonical one). We want only the
canonical Space Ranger output at `outs/spatial/cytassist_image.tiff`.

zipfile reads the archive's central directory and lets us pull out just the
one member we want, so the rest of the zip is never extracted.

Note on storage: `download_from_project_store` does NOT copy the zip onto
writable disk. It exposes the file through HISE's read-only FUSE cache mount
under `~/input/<job_id>/<store>/<uuid>/<file_name>`, so the large archive
never consumes local disk and cannot (and need not) be deleted afterwards.
"""

import glob
import os
import zipfile

import hisepy as hp

STORE_NAME = "PCD_Spatial"
FILE_NAME = "AIFI-2024-06-10T23:08:48.250692659Z/TIS05393_EXP00971_AACR_SB5.zip"
# Match the canonical Space Ranger output, not the SPATIAL_REPORTER copies.
TARGET_SUFFIX = "outs/spatial/cytassist_image.tiff"
TARGET_BASENAME = "cytassist_image.tiff"
# Extracted image is written next to this script.
OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def locate_downloaded_zip():
    """Find where hisepy placed the zip.

    The download lands in the read-only FUSE cache at
    `~/input/<job_id>/<store>/<uuid>/<FILE_NAME>`; the job id and uuid are not
    knowable ahead of time, so glob for the file and take the newest match.
    """
    pattern = os.path.expanduser(os.path.join("~", "input", "**", FILE_NAME))
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        return None
    return max(matches, key=os.path.getmtime)


def main():
    # 1. Download the single zip from the project store.
    print(f"Downloading: {FILE_NAME}")
    ok = hp.download_from_project_store(store_name=STORE_NAME, file_name=FILE_NAME)
    if not ok:
        print("Download reported failure.")
        return

    zip_path = locate_downloaded_zip()
    if not zip_path or not os.path.exists(zip_path):
        print("Could not locate the downloaded zip under ~/input/.")
        return
    print(f"Located zip: {zip_path}")

    # 2. Find the canonical cytassist_image.tiff member inside the zip.
    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if m.endswith(TARGET_SUFFIX)]

        if not members:
            print(f"No '{TARGET_SUFFIX}' inside the zip. Candidates:")
            for m in zf.namelist():
                if m.endswith(TARGET_BASENAME):
                    print(f"  {m}")
            return
        if len(members) > 1:
            print("Multiple matches; review:")
            for m in members:
                print(f"  {m}")
            return

        member = members[0]
        # 3. Extract just that one member, writing it next to this script with
        #    a sample-prefixed name so it is unambiguous (e.g.
        #    TIS05393_EXP00971_AACR_SB5_cytassist_image.tiff).
        sample = os.path.splitext(os.path.basename(FILE_NAME))[0]
        out_path = os.path.join(OUT_DIR, f"{sample}_{TARGET_BASENAME}")
        with zf.open(member) as src, open(out_path, "wb") as dst:
            dst.write(src.read())

    print(f"Extracted '{member}' -> {out_path}")

    # 4. No cleanup needed: the zip lives in HISE's read-only FUSE cache and
    #    never consumed writable disk, so there is nothing to delete.


if __name__ == "__main__":
    main()
