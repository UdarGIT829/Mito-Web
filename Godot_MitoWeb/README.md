# Mito Web

Mitochondrial variant analysis that runs entirely in the browser. It is built
with Godot 4.7 (GDScript, single-threaded web export) and hosted on GitHub
Pages. There is no server: uploaded files and projects stay in the user's
browser (`user://`, which is backed by IndexedDB). The exported
`*.mitoweb.json` workspace file is the durable save.

The Python code elsewhere in this repository (`mito_viewer/`, `main.py`,
`web_viewer.py`) is the reference implementation being ported.

## Layout

| Path | Purpose |
|---|---|
| `scripts/workspace_model.gd` | Workspace data, validation, JSON format (versioned) |
| `scripts/autoload/workspace.gd` | Open workspace, autosave to `user://` |
| `scripts/autoload/web_files.gd` | Browser upload/download bridge; desktop file-dialog fallback |
| `scripts/main.gd` | Main screen |
| `tests/run_tests.gd` | Headless tests |

## Commands

From this folder, with `godot` pointing at Godot 4.7.1:

```bash
godot --headless --script res://tests/run_tests.gd               # tests
godot --headless --export-release "Web" ../build/web/index.html  # web build
python3 -m http.server -d ../build/web 8000                      # try it locally
```

Pushing changes under `Godot_MitoWeb/` to `main` runs
`.github/workflows/pages.yml`, which tests, exports and deploys to Pages.
