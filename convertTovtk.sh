#!/bin/bash

# Convert glisX octbin results into VTK (.vtu) files for ParaView.
#
# Run from the build directory after a simulation, e.g.
#   ../convertTovtk.sh 4 . single
#   ../convertTovtk.sh 4 . two

if [ $# -ne 3 ]; then
    echo "Usage: $0 <NPROC> <folder_with_results_dir> <single|two>"
    echo "Example: $0 4 . single"
    exit 1
fi

NPROC=$1
CURRENT_WD=$2
MODEL=$3

ORIGINAL_WD=$(pwd)

cd "$CURRENT_WD" || { echo "Failed to cd into $CURRENT_WD"; exit 1; }

case "$MODEL" in
  single)
    FIELD_BASENAMES="{'swe_h_%4.4d_%4.4d', 'swe_Ux_%4.4d_%4.4d', 'swe_Uy_%4.4d_%4.4d', 'swe_Z_%4.4d_%4.4d'}"
    FIELD_NAMES="{'h', 'Ux', 'Uy', 'Z'}"
    MESH_BASENAME="swe_h_%4.4d_%4.4d"
    ;;
  two)
    FIELD_BASENAMES="{'swe_hw_%4.4d_%4.4d', 'swe_hs_%4.4d_%4.4d', 'swe_Uxw_%4.4d_%4.4d', 'swe_Uyw_%4.4d_%4.4d', 'swe_Uxs_%4.4d_%4.4d', 'swe_Uys_%4.4d_%4.4d', 'swe_Z_%4.4d_%4.4d'}"
    FIELD_NAMES="{'hw', 'hs', 'Uxw', 'Uyw', 'Uxs', 'Uys', 'Z'}"
    MESH_BASENAME="swe_hw_%4.4d_%4.4d"
    ;;
  *)
    echo "Unknown model '$MODEL' (expected 'single' or 'two')"
    exit 1
    ;;
esac

"${ORIGINAL_WD}/_deps/octave/install/bin/octave" --no-gui <<EOF
addpath(canonicalize_file_name('${ORIGINAL_WD}/../scripts/m/'))
export_tmesh_data('${MESH_BASENAME}', ...
    ${FIELD_BASENAMES}, ...
    ${FIELD_NAMES}, ...
    {}, {}, 'glisX', 0:1000, $NPROC)
EOF
