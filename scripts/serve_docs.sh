#!/usr/bin/env bash
set -e

cd $(dirname ${0})/../doc/api

python3 -m http.server 8080
