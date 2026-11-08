extends Node
## Autoload "WebFiles": opens files from, and saves files to, the user's
## computer. In a web export this goes through a small JavaScript bridge
## (<input type=file> and browser downloads); on desktop it falls back to
## native file dialogs so the app can be run from the editor.

## Emitted once per pick. `files` is an Array of
## {"name": String, "bytes": PackedByteArray}; `purpose` is the caller's tag.
signal files_picked(purpose: String, files: Array)
## usage/quota are bytes, or -1 when the browser does not report them.
signal storage_estimated(usage: int, quota: int, persisted: bool)

const JS_SETUP := """
window.mitoWeb = window.mitoWeb || (function () {
	const pending = [];
	return {
		pick: function (accept, multiple, onDone) {
			const input = document.createElement('input');
			input.type = 'file';
			input.accept = accept;
			input.multiple = multiple;
			input.addEventListener('change', async function () {
				const files = Array.from(input.files || []);
				for (const f of files) {
					pending.push({ name: f.name, bytes: new Uint8Array(await f.arrayBuffer()) });
				}
				onDone(files.length);
			});
			input.click();
		},
		nextName: function () { return pending.length ? pending[0].name : ''; },
		takeBytes: function () { const f = pending.shift(); return f ? f.bytes : new Uint8Array(0); },
		estimate: function (requestPersist, onDone) {
			const s = navigator.storage;
			if (!s || !s.estimate) { onDone(-1, -1, false); return; }
			const persist = requestPersist && s.persist ? s.persist()
				: (s.persisted ? s.persisted() : Promise.resolve(false));
			Promise.all([persist, s.estimate()])
				.then(function (r) { onDone(r[1].usage || 0, r[1].quota || 0, !!r[0]); })
				.catch(function () { onDone(-1, -1, false); });
		}
	};
})();
"""

var _is_web := OS.has_feature("web")
var _js: JavaScriptObject
# JavaScriptBridge callbacks are garbage-collected unless a reference is kept.
var _pick_callback: JavaScriptObject
var _estimate_callback: JavaScriptObject
var _pick_purpose := ""
var _dialog: FileDialog


func _ready() -> void:
	if _is_web:
		JavaScriptBridge.eval(JS_SETUP, true)
		_js = JavaScriptBridge.get_interface("mitoWeb")
		_pick_callback = JavaScriptBridge.create_callback(_on_web_pick_done)
		_estimate_callback = JavaScriptBridge.create_callback(_on_web_estimate_done)


## `accept` is a comma-separated list of extensions, e.g. ".vcf,.csv".
func pick_files(purpose: String, accept: String, multiple: bool) -> void:
	_pick_purpose = purpose
	if _is_web:
		_js.pick(accept, multiple, _pick_callback)
	else:
		_open_desktop_dialog(accept, multiple)


func save_file(file_name: String, bytes: PackedByteArray, mime: String = "application/json") -> void:
	if _is_web:
		JavaScriptBridge.download_buffer(bytes, file_name, mime)
		return
	var dialog := _new_dialog()
	dialog.file_mode = FileDialog.FILE_MODE_SAVE_FILE
	dialog.current_file = file_name
	dialog.file_selected.connect(func(path: String) -> void:
		var file := FileAccess.open(path, FileAccess.WRITE)
		if file:
			file.store_buffer(bytes)
	)
	dialog.popup_centered_ratio(0.6)


## Reports browser storage use. Pass request_persist=true after the user adds
## data, to ask the browser not to evict it (some browsers show a prompt).
func estimate_storage(request_persist: bool = false) -> void:
	if _is_web:
		_js.estimate(request_persist, _estimate_callback)
	else:
		storage_estimated.emit(-1, -1, true)


func _on_web_pick_done(args: Array) -> void:
	var files: Array = []
	for i in int(args[0]):
		var file_name := str(JavaScriptBridge.eval("window.mitoWeb.nextName()", true))
		var bytes = JavaScriptBridge.eval("window.mitoWeb.takeBytes()", true)
		if typeof(bytes) != TYPE_PACKED_BYTE_ARRAY:
			bytes = PackedByteArray()
		files.append({"name": file_name, "bytes": bytes})
	files_picked.emit(_pick_purpose, files)


func _on_web_estimate_done(args: Array) -> void:
	storage_estimated.emit(int(args[0]), int(args[1]), bool(args[2]))


func _open_desktop_dialog(accept: String, multiple: bool) -> void:
	var dialog := _new_dialog()
	dialog.file_mode = FileDialog.FILE_MODE_OPEN_FILES if multiple else FileDialog.FILE_MODE_OPEN_FILE
	var patterns := PackedStringArray()
	for ext in accept.split(",", false):
		patterns.append("*" + ext.strip_edges())
	if not patterns.is_empty():
		dialog.add_filter(", ".join(patterns), "Supported files")
	var deliver := func(paths: PackedStringArray) -> void:
		var files: Array = []
		for path in paths:
			files.append({"name": path.get_file(), "bytes": FileAccess.get_file_as_bytes(path)})
		files_picked.emit(_pick_purpose, files)
	dialog.files_selected.connect(deliver)
	dialog.file_selected.connect(func(path: String) -> void: deliver.call(PackedStringArray([path])))
	dialog.popup_centered_ratio(0.6)


func _new_dialog() -> FileDialog:
	if _dialog:
		_dialog.queue_free()
	_dialog = FileDialog.new()
	_dialog.access = FileDialog.ACCESS_FILESYSTEM
	_dialog.use_native_dialog = true
	get_tree().root.add_child(_dialog)
	return _dialog
