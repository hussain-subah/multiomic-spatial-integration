#!/bin/bash
#
# Download commonly used Enrichr GMT gene-set libraries.
#
# Usage:
#   bash resources/download_gmt_resources.sh
#

set -euo pipefail

OUTDIR="resources/gmt"
mkdir -p "${OUTDIR}"

BASE_URL="https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="

download_gmt () {
    LIBRARY="$1"
    OUTFILE="$2"

    echo "Downloading ${LIBRARY}..."
    curl -L \
        "${BASE_URL}${LIBRARY}" \
        -o "${OUTDIR}/${OUTFILE}"
}

download_gmt "KEGG_2021_Human"        "KEGG_2021_Human.gmt"
download_gmt "Reactome_2022"          "Reactome_2022.gmt"
download_gmt "WikiPathway_2023_Human" "WikiPathways_2023_Human.gmt"

echo
echo "Downloaded GMT files:"
ls -lh "${OUTDIR}"
