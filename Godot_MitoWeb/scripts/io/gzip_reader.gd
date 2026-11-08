class_name GzipReader
extends RefCounted
## Decompresses gzip and BGZF (blocked gzip, used by .vcf.gz from bgzip/htslib).
## zlib stops at the end of the first gzip member, so BGZF files are split
## into their blocks using the BSIZE field of each block header.


static func is_gzip(bytes: PackedByteArray) -> bool:
	return bytes.size() >= 2 and bytes[0] == 0x1f and bytes[1] == 0x8b


## Returns {"ok": true, "bytes": PackedByteArray} or {"ok": false, "error": String}.
static func decompress(bytes: PackedByteArray) -> Dictionary:
	if not is_gzip(bytes):
		return {"ok": true, "bytes": bytes}
	if _bgzf_block_size(bytes, 0) <= 0:
		var single := bytes.decompress_dynamic(-1, FileAccess.COMPRESSION_GZIP)
		if single.is_empty() and bytes.size() > 20:
			return {"ok": false, "error": "The gzip data could not be decompressed."}
		return {"ok": true, "bytes": single}
	var out := PackedByteArray()
	var offset := 0
	while offset < bytes.size():
		var block_size := _bgzf_block_size(bytes, offset)
		if block_size <= 0 or offset + block_size > bytes.size():
			return {"ok": false, "error": "Corrupt BGZF block at byte %d." % offset}
		var block := bytes.slice(offset, offset + block_size)
		var isize := block.decode_u32(block_size - 4)
		if isize > 0:
			var data := block.decompress_dynamic(-1, FileAccess.COMPRESSION_GZIP)
			if data.size() != isize:
				return {"ok": false, "error": "BGZF block at byte %d did not decompress correctly." % offset}
			out.append_array(data)
		offset += block_size
	return {"ok": true, "bytes": out}


## Total size of the BGZF block starting at `offset`, or -1 if it is not one.
static func _bgzf_block_size(bytes: PackedByteArray, offset: int) -> int:
	if offset + 18 > bytes.size() or bytes[offset] != 0x1f or bytes[offset + 1] != 0x8b:
		return -1
	if not (bytes[offset + 3] & 0x04):  # FEXTRA
		return -1
	var xlen := bytes.decode_u16(offset + 10)
	var pos := offset + 12
	var end := pos + xlen
	while pos + 4 <= end:
		var subfield_len := bytes.decode_u16(pos + 2)
		if bytes[pos] == 66 and bytes[pos + 1] == 67 and subfield_len == 2:  # "BC"
			return bytes.decode_u16(pos + 4) + 1
		pos += 4 + subfield_len
	return -1
