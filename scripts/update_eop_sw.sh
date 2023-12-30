#!/usr/bin/env bash
set -e
cd $(dirname ${0})/..

RATE=256K

rm -f external/EOP-All.csv external/SW-All.csv

curl --limit-rate ${RATE} \
    "https://celestrak.org/SpaceData/EOP-All.csv" \
    > external/EOP-All.csv

curl --limit-rate ${RATE} \
    "https://celestrak.org/SpaceData/SW-All.csv" \
    > external/SW-All.csv
