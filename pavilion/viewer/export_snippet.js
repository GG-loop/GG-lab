// ---- 3D download (shared by every viewer page; build scripts inline this file)
// Writes binary STL in metres with Z up (the generator's frame; scenes are Y up),
// packs it with any extra files into a .zip and offers it through the viewer's
// `downloads` capability. The viewer confirms every save; the button stays hidden
// where downloads are not available.
const Export = (() => {
  let ns = null;
  const ready = (window.claude && typeof window.claude.use === "function")
    ? window.claude.use("downloads").then((d) => { ns = d; return !!d; }, () => false)
    : Promise.resolve(false);

  // items: [{ geometry: THREE.BufferGeometry, matrix?: THREE.Matrix4 }]
  function stl(items, header) {
    let n = 0;
    for (const { geometry: g } of items) n += (g.index ? g.index.count : g.attributes.position.count) / 3;
    const buf = new ArrayBuffer(84 + 50 * n), dv = new DataView(buf);
    const head = (header || "binary STL, metres, Z up").slice(0, 79);
    for (let i = 0; i < head.length; i++) dv.setUint8(i, head.charCodeAt(i) & 0x7f);
    dv.setUint32(80, n, true);
    const a = new THREE.Vector3(), b = new THREE.Vector3(), c = new THREE.Vector3();
    const u = new THREE.Vector3(), w = new THREE.Vector3();
    let o = 84;
    const put = (v) => { dv.setFloat32(o, v.x, true); dv.setFloat32(o + 4, -v.z, true); dv.setFloat32(o + 8, v.y, true); o += 12; };
    for (const { geometry: g, matrix: m } of items) {
      const p = g.attributes.position, ix = g.index, cnt = ix ? ix.count : p.count;
      for (let t = 0; t < cnt; t += 3) {
        a.fromBufferAttribute(p, ix ? ix.getX(t) : t);
        b.fromBufferAttribute(p, ix ? ix.getX(t + 1) : t + 1);
        c.fromBufferAttribute(p, ix ? ix.getX(t + 2) : t + 2);
        if (m) { a.applyMatrix4(m); b.applyMatrix4(m); c.applyMatrix4(m); }
        u.subVectors(b, a); w.subVectors(c, a); u.cross(w).normalize();
        put(u); put(a); put(b); put(c);
        dv.setUint16(o, 0, true); o += 2;
      }
    }
    return buf;
  }

  // ---- minimal zip writer: deflate through the browser's CompressionStream, else stored
  const CRC = (() => {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; }
    return t;
  })();
  function crc32(u8) {
    let c = 0xffffffff;
    for (let i = 0; i < u8.length; i++) c = CRC[(c ^ u8[i]) & 0xff] ^ (c >>> 8);
    return (c ^ 0xffffffff) >>> 0;
  }
  async function deflate(u8) {
    if (typeof CompressionStream !== "function") return null;
    try {
      const stream = new Blob([u8]).stream().pipeThrough(new CompressionStream("deflate-raw"));
      return new Uint8Array(await new Response(stream).arrayBuffer());
    } catch (e) { return null; }
  }
  async function zip(files) {
    const parts = [], central = [];
    let offset = 0;
    for (const [name, value] of Object.entries(files)) {
      const data = typeof value === "string" ? new TextEncoder().encode(value) : new Uint8Array(value);
      const fname = new TextEncoder().encode(name);
      const packed = await deflate(data);
      const method = packed && packed.length < data.length ? 8 : 0;
      const body = method ? packed : data;
      const crc = crc32(data);
      const local = new DataView(new ArrayBuffer(30));
      local.setUint32(0, 0x04034b50, true); local.setUint16(4, 20, true); local.setUint16(6, 0x0800, true);
      local.setUint16(8, method, true); local.setUint16(10, 0, true); local.setUint16(12, 0x21, true);
      local.setUint32(14, crc, true); local.setUint32(18, body.length, true); local.setUint32(22, data.length, true);
      local.setUint16(26, fname.length, true); local.setUint16(28, 0, true);
      const cen = new DataView(new ArrayBuffer(46));
      cen.setUint32(0, 0x02014b50, true); cen.setUint16(4, 20, true); cen.setUint16(6, 20, true); cen.setUint16(8, 0x0800, true);
      cen.setUint16(10, method, true); cen.setUint16(12, 0, true); cen.setUint16(14, 0x21, true);
      cen.setUint32(16, crc, true); cen.setUint32(20, body.length, true); cen.setUint32(24, data.length, true);
      cen.setUint16(28, fname.length, true); cen.setUint32(42, offset, true);
      parts.push(local.buffer, fname, body);
      central.push(cen.buffer, fname);
      offset += 30 + fname.length + body.length;
    }
    const size = central.reduce((a, b) => a + b.byteLength, 0);
    const end = new DataView(new ArrayBuffer(22));
    end.setUint32(0, 0x06054b50, true);
    end.setUint16(8, Object.keys(files).length, true); end.setUint16(10, Object.keys(files).length, true);
    end.setUint32(12, size, true); end.setUint32(16, offset, true);
    return new Blob([...parts, ...central, end.buffer], { type: "application/zip" });
  }

  // files: { "name.ext": ArrayBuffer | string }
  async function save(name, files) {
    if (!(await ready) || !ns) throw { code: "unavailable", message: "Downloads are not available here." };
    return ns.save({ filename: name + ".zip", data: await zip(files) });
  }

  // wire a button + status element; build() returns { name, files }
  function bind(button, status, build) {
    ready.then((ok) => { button.hidden = !ok; });
    button.addEventListener("click", async () => {
      button.disabled = true;
      status.textContent = "Preparing the 3D file…";
      try {
        const { name, files } = await build();
        const r = await save(name, files);
        status.textContent = r.status === "saved" ? "Saved." : "Sent.";
      } catch (e) {
        const code = e && e.code;
        status.textContent = code === "declined" ? "Download cancelled."
          : code === "rate_limited" ? "Another download is waiting for confirmation."
          : code === "unavailable" || code === "not_granted" ? "Downloads are not available in this view."
          : "The download failed: " + ((e && e.message) || "unknown error") + ".";
        if (code === "unavailable" || code === "not_granted") button.hidden = true;
      } finally {
        button.disabled = false;
      }
    });
  }

  const readme = (title, lines) => [title, "", "Units: metres. Up axis: Z. Binary STL.", ...lines, ""].join("\r\n");
  return { ready, stl, save, bind, readme };
})();
