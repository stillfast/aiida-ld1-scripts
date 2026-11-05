#!/usr/bin/env bash
# Export the results group (and its role subgroups) to a self-contained AiiDA archive.
#
# The archive holds every calculation, its inputs (the recipe as a Dict node), the
# generated UPF and the retrieved ld1.x output, so the pseudopotential can be audited
# or re-analysed from the file alone:
#
#     verdi archive import data/aiida-ld1-6p25.aiida
#     verdi group list
#
# Note the argument order: `verdi archive create` takes the OUTPUT FILE positionally and
# `-G/--groups` takes one or more labels, so the file has to come first.  The role
# subgroups are passed explicitly -- they are what records which calculation is the
# generation, which are the rlderiv points and which is the transferability test, and
# without them the archive would lose that information.
set -euo pipefail

GROUP="${1:-aiida-ld1-6p25}"
OUT="${2:-$(dirname "$0")/../data/${GROUP}.aiida}"

echo "exporting '${GROUP}' (+ .generation, .rlderiv_scan, .configts) -> ${OUT}"
verdi archive create --force "${OUT}" -G \
    "${GROUP}" "${GROUP}.generation" "${GROUP}.rlderiv_scan" "${GROUP}.configts"
ls -lh "${OUT}"
echo
echo "verify with:  verdi archive info ${OUT}"
