#!/usr/bin/env python3
"""
Build the interactive 3D viewer (viewer/index.html).

The page lays modules out in the browser, so it only needs one mesh per
surface: the bare minimal surface for the assembly (light enough to instance
hundreds of times) and the thickened, printable module for the detail view.
Both are quantised to 16 bit and inlined into viewer/template.html.

    python build_viewer.py
"""

import base64
import json
from pathlib import Path

import numpy as np

import pavilion as pv

HERE = Path(__file__).parent
SURFACES = ["gyroid", "schwarz_p", "diamond", "neovius"]


def pack(mesh, size):
    """Positions as uint16 over [0, size], indices as uint16/uint32, all base64."""
    q = np.clip(np.round(mesh.vertices / size * 65535), 0, 65535).astype("<u2")
    idx_type = "<u2" if len(mesh.vertices) < 65536 else "<u4"
    idx = mesh.faces.astype(idx_type)
    return {
        "pos": base64.b64encode(q.tobytes()).decode(),
        "idx": base64.b64encode(idx.tobytes()).decode(),
        "idx32": idx_type == "<u4",
    }


def main():
    cfg = pv.DEFAULT_CONFIG
    M = cfg["module_size_m"]
    t_ref = 0.04  # connector pattern and solid preview at a mid-range sheet
    data = {"module_size_m": M, "defaults": {
        "sheet_thickness_m": cfg["sheet_thickness_m"],
        "density_kg_m3": cfg["printer"]["density_kg_m3"],
        "material": cfg["printer"]["material"],
        "build_volume_m": cfg["printer"]["build_volume_m"],
        "max_roof_span_modules": cfg["max_roof_span_modules"],
        "envelope_depth_modules": cfg["envelope_depth_modules"],
    }, "surfaces": {}}

    for s in SURFACES:
        c = pv.merge(cfg, {"surface": s})
        mid = pv.mid_surface(c, 32)
        solid = pv.sheet_solid(c, cfg["sheet_thickness_m"]["base"], [0, 0, 0], [M] * 3, 48)
        conn = {ax: pv.face_connectors(c, t_ref, a) for a, ax in enumerate("xyz")}
        data["surfaces"][s] = {
            "area_m2": round(float(mid.area), 4),  # mid-surface area per module
            "mid": pack(mid, M),
            "solid": pack(solid, M),
            "solid_thickness_m": cfg["sheet_thickness_m"]["base"],
            "connectors": conn,
        }
        print(f"{s}: mid {len(mid.faces)} tris, solid {len(solid.faces)} tris, "
              f"{sum(len(v) for v in conn.values())} connector points on 3 faces")

    html = (HERE / "viewer" / "template.html").read_text()
    html = html.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
    (HERE / "viewer" / "index.html").write_text(html)
    print(f"wrote viewer/index.html ({len(html) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
