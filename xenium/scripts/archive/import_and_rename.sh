#!/bin/bash
# For each sample: download its xenium bundle from HISE, rewrite the wrong sample/slide
# IDs with `xeniumranger rename`, upload the corrected bundle to the staging bucket, and
# move any wrong-named copy already in the bucket to wrong_id/. A failed sample is skipped
# and reported.

set -uo pipefail

# "WRONG_SAMPLE WRONG_SLIDE CORRECT_SAMPLE CORRECT_SLIDE HISE_UUID"
SAMPLES=(
    "TMA00145 0031088 TMA00142 0031088 56afab01-d275-420e-8e25-41a900098ce0"
    "TMA00142 0031088 TMA00145 0031088 f2a09064-72be-42ef-b8a6-99fc592c1fef"
)

BUCKET="gs://temp-xenium-hise-transfer"
BUNDLE_BASE="/home/workspace/xenium/data/bundles"
OUTPUT_BASE="/home/workspace/xenium/data/renamed"
XENIUMRANGER="/home/workspace/xeniumranger-xenium4.0/xeniumranger"
PYTHON="/home/workspace/environment/scanpyv3/bin/python"

mkdir -p "$BUNDLE_BASE" "$OUTPUT_BASE"
bucket=$(gcloud storage ls "$BUCKET/")

echo "Plan:"
rows=()
for row in "${SAMPLES[@]}"; do
    read -r ws wl cs cl uuid <<< "$row"
    echo "  $ws/$wl → $cs/$cl   (HISE UUID: $uuid)"
    rows+=("$ws|$wl|$cs|$cl|$uuid")
done

read -r -p "Proceed? This mutates $BUCKET. [y/N] " ans
[[ "$ans" =~ ^[Yy] ]] || { echo "Aborted."; exit 0; }

ok=(); fail=(); uploaded_finalnames=()
for r in "${rows[@]}"; do
    IFS='|' read -r ws wl cs cl uuid <<< "$r"
    bundle="$BUNDLE_BASE/${ws}_${wl}"
    run="$OUTPUT_BASE/output_corrected_${cs}"

    echo; echo "=== $ws → $cs ==="

    # Download from HISE by UUID; redirect logging to stderr, take last stdout line as path
    zippath=$("$PYTHON" -c "
import unittest.mock as mock, sys, logging
logging.basicConfig(stream=sys.stderr)
sys.modules['hisepy.training_job'] = mock.MagicMock()
sys.modules['hisepy.ai_ray_transformer'] = mock.MagicMock()
import hisepy as hp
paths = hp.cache_files(file_ids=['$uuid'])
if not paths: sys.exit('cache_files returned no paths')
print(paths[0])
" | tail -1) || { echo "✗ HISE download failed for $ws"; fail+=("$ws → $cs"); continue; }

    srcname=$(basename "${zippath%.zip}")
    finalname="${srcname//__${wl}__${ws}__/__${cl}__${cs}__}"
    final="$OUTPUT_BASE/$finalname"

    if mkdir -p "$bundle" &&
       unzip -qo "$zippath" -d "$bundle/" &&
       rm -rf "$run" &&
       ( cd "$OUTPUT_BASE" && "$XENIUMRANGER" rename \
             --id="output_corrected_${cs}" \
             --xenium-bundle="$bundle/$srcname" \
             --region-name="$cs" ) &&
       rm -rf "$final" && mv "$run/outs" "$final" && rm -rf "$run" &&
       gcloud storage rsync -r "$final" "$BUCKET/$finalname/" &&
       rm -rf "$bundle" "$final"; then

        uploaded_finalnames+=("$finalname")

        # Move any pre-existing wrong-named copy in the bucket to wrong_id/,
        # but skip folders we just uploaded in this run (same name, new content).
        wrong_match=$(awk -v p="__${wl}__${ws}__" 'index($0,p) && !index($0,"wrong_id/")' <<< "$bucket")
        if [[ -n "$wrong_match" ]]; then
            wrongsrc=$(basename "${wrong_match%/}")
            skip=false
            for fn in "${uploaded_finalnames[@]}"; do [[ "$fn" == "$wrongsrc" ]] && skip=true && break; done
            if $skip; then
                echo "Skipping wrong_id move for $wrongsrc (just uploaded as corrected output)"
            else
                echo "Moving $wrongsrc → wrong_id/"
                gcloud storage mv --no-clobber "$BUCKET/$wrongsrc" "$BUCKET/wrong_id/$wrongsrc" || \
                    echo "Warning: could not move $wrongsrc to wrong_id/"
            fi
        fi
        ok+=("$cs"); echo "✓ $cs"
    else
        fail+=("$ws → $cs"); echo "✗ $ws failed (bundle kept: $bundle)"
    fi
done

echo; echo "Done: ${#ok[@]} succeeded, ${#fail[@]} failed."
if ((${#fail[@]})); then printf '  ✗ %s\n' "${fail[@]}"; exit 1; fi
