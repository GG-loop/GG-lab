#!/usr/bin/env python3
"""
Build the feasibility report page (viewer/feasibility.html) from out_analysis/.

The page shows the analysed pavilion in 3D (shell with the stair openings,
the printed stairs, utilisation colours, critical zones and the worst joints)
next to the stair compliance and the structural results.

    python analysis.py && python build_report.py
"""

import base64
import json
from pathlib import Path

import numpy as np
import trimesh
from scipy.spatial import cKDTree

HERE = Path(__file__).parent
SRC = HERE / "out_analysis"


def pack(mesh, lo, size, extra=None):
    q = np.clip(np.round((mesh.vertices - lo) / size * 65535), 0, 65535).astype("<u2")
    out = {"pos": base64.b64encode(q.tobytes()).decode(),
           "idx": base64.b64encode(mesh.faces.astype("<u4").tobytes()).decode()}
    if extra is not None:
        out["val"] = base64.b64encode(extra.astype("<u1").tobytes()).decode()
    return out


def decimate(mesh, faces):
    try:
        return mesh.simplify_quadric_decimation(face_count=min(faces, len(mesh.faces)))
    except ImportError:
        return mesh


def main():
    report = json.loads((SRC / "report.json").read_text())
    shell = trimesh.load(SRC / "pavilion_with_stairs.stl")
    util = np.load(SRC / "vertex_utilisation.npy") if (SRC / "vertex_utilisation.npy").exists() else None
    stairs = trimesh.load(SRC / "stairs.stl") if (SRC / "stairs.stl").exists() else None

    small = decimate(shell, 160000)
    if util is not None:
        _, nearest = cKDTree(shell.vertices).query(small.vertices)
        u = util[nearest]
        val = np.clip(np.round(u / 1.5 * 255), 0, 255)   # 0..1.5 utilisation -> 0..255
    else:
        val = np.zeros(len(small.vertices))
    lo = shell.bounds[0]
    size = float((shell.bounds[1] - shell.bounds[0]).max())
    data = {"lo": lo.round(4).tolist(), "size": size, "shell": pack(small, lo, size, val)}
    if stairs is not None:
        data["stairs"] = pack(decimate(stairs, 30000), lo, size)

    st = report.get("structure", {})
    keep = {k: report[k] for k in ("rules", "mesh") if k in report}
    keep["config"] = {k: report["config"][k] for k in ("surface", "cell_m", "pods", "chamber_centre_m", "sheet",
                                                        "organic", "printer")}
    if "stairs" in report:
        s = dict(report["stairs"])
        s.pop("route", None); s.pop("route_kind", None)
        keep["stairs"] = s
    if st:
        keep["structure"] = st
    data["report"] = keep

    html = (HERE / "viewer" / "feasibility_template.html").read_text()
    html = html.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":"), default=float))
    (HERE / "viewer" / "feasibility.html").write_text(html)
    print(f"wrote viewer/feasibility.html ({len(html) / 1e6:.1f} MB), shell {len(small.faces)} tris")


if __name__ == "__main__":
    main()
