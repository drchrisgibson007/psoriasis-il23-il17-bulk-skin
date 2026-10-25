# Source this file from the packet folder:  source config.sh
# Sets GRG_ROOT to this packet's top folder. All scripts resolve paths from it.
export GRG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Large public inputs that are not in the zip go under $GRG_ROOT/external (see README).
mkdir -p "$GRG_ROOT/external"
