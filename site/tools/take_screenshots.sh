#!/usr/bin/env bash
# Screenshots for the project page (site/screenshots/), from the public
# 1000 Genomes example data only.
#
# Runs a disposable copy of the app with its own user-data folder, so the
# owner's real workspace (which may hold private data) is never opened,
# changed or captured. Needs a display (the renderer runs).
#
#   site/tools/take_screenshots.sh
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
GODOT=${GODOT:-$HOME/Applications/Godot_v4.7.1-stable_linux.x86_64}
COPY=$(mktemp -d)
trap 'rm -rf "$COPY"; rm -rf "$HOME/.local/share/godot/app_userdata/MitoWebPageShots"' EXIT
cp -a "$ROOT/Godot_MitoWeb/." "$COPY/"
# A separate user-data folder for the copy.
python3 - "$COPY/project.godot" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
s = s.replace('config/name="Mito Web"', 'config/name="Mito Web"\nconfig/use_custom_user_dir=true\nconfig/custom_user_dir_name="MitoWebPageShots"', 1)
# No tooltips: the real mouse pointer may rest over the window.
s += "\n[gui]\n\ntimers/tooltip_delay_sec=100000.0\n" if "[gui]" not in s else ""
s = s.replace("[gui]\n", "[gui]\n\ntimers/tooltip_delay_sec=100000.0\n", 1) if "timers/tooltip_delay_sec" not in s else s
open(p, "w").write(s)
PY
rm -rf "$HOME/.local/share/godot/app_userdata/MitoWebPageShots"
"$GODOT" --headless --path "$COPY" --import >/dev/null 2>&1 || true
OUT="$ROOT/site/screenshots" "$GODOT" --path "$COPY" --resolution 1600x900 --script "$ROOT/site/tools/take_screenshots.gd"
