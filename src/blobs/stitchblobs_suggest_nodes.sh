#!/bin/bash
# Estimates the minimum Perlmutter CPU nodes and the total MPI ranks for StitchBlobs.
# <in_list> is the same list passed to StitchBlobs --in_list. Opens every file in it
# and samples tagged points from some of them. Runs serially; a login node is fine.
# Assumes a StitchBlobs build that includes PR #111.
#
# Usage: stitchblobs_suggest_nodes.sh <in_list> <var> [StitchBlobsNodeEstimator options]
# Run it from the directory StitchBlobs will run in, since relative paths in the
# list resolve against the current directory.

set -eo pipefail

if [[ "$1" == "-h" || "$1" == "--help" || $# -lt 2 ]]; then
	sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'
	[[ $# -ge 2 || "$1" == "-h" || "$1" == "--help" ]] && exit 0
	exit 2
fi
if [[ "$NERSC_HOST" != "perlmutter" ]]; then
	echo "ERROR: this script uses Perlmutter defaults; run it on Perlmutter." >&2
	exit 1
fi

in_list=$1
var=$2
shift 2

type module > /dev/null 2>&1 || source /etc/profile
module -q load cpu cray-hdf5 cray-netcdf

dir=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
src=$dir/StitchBlobsNodeEstimator.cpp
if [[ -f "$src" ]]; then
	# Running from a source checkout: rebuild when the helper changed.
	bin=$dir/../../bin/StitchBlobsNodeEstimator
	if [[ ! -x "$bin" || "$src" -nt "$bin" ]]; then
		echo "Building StitchBlobsNodeEstimator..." >&2
		mkdir -p "$(dirname "$bin")"
		nc=$dir/../netcdf-cxx-4.2
		CC -O2 -std=c++11 -I"$nc" "$src" "$nc/netcdf.cpp" "$nc/ncvalues.cpp" -o "$bin"
	fi
else
	# Installed: the helper sits next to this script.
	bin=$dir/StitchBlobsNodeEstimator
fi

exec "$bin" --in_list "$in_list" --var "$var" "$@"
