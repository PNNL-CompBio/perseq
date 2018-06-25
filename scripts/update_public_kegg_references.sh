#!/usr/bin/env bash
set -euo pipefail

outdir=${1:-"tests/refdata"}

curl "http://rest.kegg.jp/list/{ko,pathway}" \
    --create-dirs --output "${outdir}/#1_list.txt"
curl "http://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json" \
    --location --create-dirs --output "${outdir}/hierarchy.json"
