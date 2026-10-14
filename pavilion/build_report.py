#!/usr/bin/env python3
"""
Build the feasibility report page (viewer/feasibility.html) from out_analysis/.

The page shows the analysed pavilion in 3D (shell with the stair openings,
the printed stairs, utilisation colours, critical zones and the worst joints)
next to the stair compliance and the structural results.

    python analysis.py && python build_report.py
    python analysis.py --art && python build_report.py --sculpture   # the climbable sculpture
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


def main_sculpture(src=HERE / "out_sculpture"):
    report = json.loads((src / "report.json").read_text())
    mesh = trimesh.load(src / "sculpture.stl")
    small = decimate(mesh, 260000)
    _, nearest = cKDTree(mesh.vertices).query(small.vertices)
    util = np.load(src / "vertex_utilisation.npy")[nearest] if (src / "vertex_utilisation.npy").exists() \
        else np.zeros(len(small.vertices))
    cls = np.load(src / "vertex_class.npy")[nearest]
    if "--reclassify" in __import__("sys").argv:
        # recompute the ledge classes without rerunning the structural model
        import analysis as an
        import organic_pavilion as op
        from pavilion import merge
        cfg = report["config"]
        F0, axes = op.field(merge(cfg, {"terraces": {"enabled": False}}))
        cls = an.ledge_classes(small, op.to_mesh(F0, axes, cfg["voxel_m"]))
    lo = mesh.bounds[0]
    size = float((mesh.bounds[1] - mesh.bounds[0]).max())
    shell = pack(small, lo, size, np.clip(np.round(util / 1.5 * 255), 0, 255))
    shell["cls"] = base64.b64encode(cls.astype("<u1").tobytes()).decode()
    data = {"lo": lo.round(4).tolist(), "size": size, "shell": shell}
    if "membranes" in report and (src / "membranes.stl").exists():
        mem = trimesh.load(src / "membranes.stl")
        data["membranes"] = pack(decimate(mem, 60000), lo, size)
    keep = {k: report[k] for k in ("climb", "mesh", "print", "structure", "membranes") if k in report}
    keep["config"] = {k: report["config"][k] for k in ("surface", "cell_m", "chamber_centre_m", "sheet", "organic",
                                                        "printer", "terraces")}
    data["report"] = keep
    variants = []
    for name, d in [("as drawn", src), ("thicker sheet", HERE / "out_sculpture_thick"),
                    ("edge beams", HERE / "out_sculpture_rim")]:
        f = d / "report.json"
        if f.exists():
            r = json.loads(f.read_text())
            st = r.get("structure", {})
            variants.append({"name": name, "sheet": r["config"]["sheet"], "mass_kg": r["mesh"]["mass_kg"],
                             "mass_hollow_ledges_kg": r["mesh"]["mass_hollow_ledges_kg"],
                             "disp_mm": st.get("max_displacement_mm"), "limit_mm": st.get("displacement_limit_mm"),
                             "util_max": st.get("utilisation", {}).get("max")})
    keep["variants"] = variants
    html = (HERE / "viewer" / "sculpture_template.html").read_text()
    html = html.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":"), default=float))
    (HERE / "viewer" / "sculpture.html").write_text(html)
    print(f"wrote viewer/sculpture.html ({len(html) / 1e6:.1f} MB), {len(small.faces)} tris")


def main():
    import sys
    if "--sculpture" in sys.argv:
        return main_sculpture()
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
