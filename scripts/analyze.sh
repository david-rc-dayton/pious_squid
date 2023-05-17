#!/usr/bin/env bash
set -e

cd $(dirname ${0})/..

dart fix --apply .
dart format --fix .
dart analyze .
