"""
Download a single zip from the HISE project store `PCD_Spatial` and extract
ONLY `web_summary.html` from it (without unpacking the whole archive).

Target file (exact key/name):
    AIFI-2024-06-10T23:22:51.980946983Z/EXP00897_BM.zip

zipfile reads the archive's central directory and lets us pull out just the
one member we want, so the rest of the zip is never extracted.
"""

import os
import zipfile

import hisepy as hp

STORE_NAME = "PCD_Spatial"
FILE_NAME = "AIFI-2024-06-10T23:22:51.980946983Z/EXP00897_BM.zip"
TARGET_BASENAME = "web_summary.html"


def main():
    # 1. Download the single zip from the project store.
    print(f"Downloading: {FILE_NAME}")
    ok = hp.download_from_project_store(store_name=STORE_NAME, file_name=FILE_NAME)
    if not ok:
        print("Download reported failure.")
        return

    zip_path = os.path.expanduser(os.path.join("~", STORE_NAME, FILE_NAME))
    if not os.path.exists(zip_path):
        print(f"Expected zip not found at: {zip_path}")
        return

    # 2. Find the web_summary.html member inside the zip (path inside the
    #    archive is unknown, so match on basename).
    out_dir = os.path.dirname(zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if m.endswith(TARGET_BASENAME)]

        if not members:
            print(f"No '{TARGET_BASENAME}' inside the zip. Members:")
            for m in zf.namelist():
                print(f"  {m}")
            return
        if len(members) > 1:
            print("Multiple matches; review:")
            for m in members:
                print(f"  {m}")
            return

        member = members[0]
        # 3. Extract just that one member, writing it flat into out_dir.
        out_path = os.path.join(out_dir, TARGET_BASENAME)
        with zf.open(member) as src, open(out_path, "wb") as dst:
            dst.write(src.read())

    print(f"Extracted '{member}' -> {out_path}")


if __name__ == "__main__":
    main()
