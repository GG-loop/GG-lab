#!/usr/bin/env python3
"""
Modular minimal-surface pavilion generator.

The pavilion is a voxel grid of cubic modules. Every module is exactly one
period of a triply periodic minimal surface (TPMS), thickened into a printable
sheet. Because the surface repeats with the module size, every module of the
same sheet thickness is geometrically identical, and neighbouring modules meet
on flat faces where the surface continues without a seam. The pavilion grows
by adding modules; it never needs new geometry.

    python pavilion.py                      # default 6 x 6 x 6 m pavilion
    python pavilion.py -c my_config.json    # custom extent / surface / printer
    python pavilion.py --size 9 6 6         # quick override of the extent (m)

Outputs (in out/ by default):
    modules/<type>.stl      one watertight print file per unique module type
    parts/<type>_pNN.stl    the same type split to fit the printer build volume
    assembly_preview.stl    low resolution assembly for viewing
    manifest.json           placement, roles, connectors, quantities, mass
    preview.png             axonometric render + plan of the module layout
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import ndimage
from skimage import measure
import trimesh

MIN_EXTENT_M = 6.0

DEFAULT_CONFIG = {
    "surface": "gyroid",            # gyroid | schwarz_p | diamond | neovius
    "module_size_m": 1.5,           # edge of one cubic module = one TPMS period
    "size_m": [6.0, 6.0, 6.0],      # requested extent, rounded up to whole modules
    "envelope_depth_modules": 1,    # thickness of the walls / roof in modules
    "sheet_thickness_m": {          # graded from ground level to roof level
        "base": 0.050,
        "top": 0.030,
        "levels": "per_level"       # "per_level" or "uniform" (uses base everywhere)
    },
    "entrances": [                  # openings cut out of the ground-level walls
        {"side": "south", "width_modules": 2, "height_modules": 2},
        {"side": "north", "width_modules": 2, "height_modules": 2}
    ],
    "max_roof_span_modules": 4,
    "connectors": {"spacing_m": 0.30, "edge_distance_m": 0.06},     # add interior columns when the clear span is larger
    "resolution": {"module": 64, "preview": 20},  # samples per module edge
    "printer": {
        "build_volume_m": [1.0, 1.0, 1.0],
        "material": "rPETG",
        "density_kg_m3": 1270
    },
    "output_dir": "out"
}

SIDES = ("south", "north", "west", "east")


# --------------------------------------------------------------------------- #
# Minimal surfaces
# --------------------------------------------------------------------------- #

def tpms(name, x, y, z):
    """Implicit TPMS value; x, y, z already scaled to one period = 2*pi."""
    s, c = np.sin, np.cos
    if name == "gyroid":
        return s(x) * c(y) + s(y) * c(z) + s(z) * c(x)
    if name == "schwarz_p":
        return c(x) + c(y) + c(z)
    if name == "diamond":
        return (s(x) * s(y) * s(z) + s(x) * c(y) * c(z)
                + c(x) * s(y) * c(z) + c(x) * c(y) * s(z))
    if name == "neovius":
        return 3 * (c(x) + c(y) + c(z)) + 4 * c(x) * c(y) * c(z)
    raise ValueError(f"unknown surface '{name}'")


def sheet_distance(name, period, xs, ys, zs):
    """Approximate distance (m) to the minimal surface: f / |grad f|."""
    k = 2 * math.pi / period
    X, Y, Z = np.meshgrid(xs * k, ys * k, zs * k, indexing="ij")
    f = tpms(name, X, Y, Z)
    h = xs[1] - xs[0]
    gx, gy, gz = np.gradient(f, h)
    grad = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    return np.abs(f) / np.maximum(grad, 1e-6)


def box_sdf(xs, ys, zs, lo, hi):
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    c = (np.asarray(lo) + np.asarray(hi)) / 2
    r = (np.asarray(hi) - np.asarray(lo)) / 2
    q = np.stack([np.abs(X - c[0]) - r[0], np.abs(Y - c[1]) - r[1], np.abs(Z - c[2]) - r[2]])
    outside = np.linalg.norm(np.maximum(q, 0), axis=0)
    inside = np.minimum(q.max(axis=0), 0)
    return outside + inside


def sheet_solid(cfg, thickness, lo, hi, samples_per_module):
    """Watertight mesh of the thickened TPMS clipped to the box [lo, hi]."""
    M = cfg["module_size_m"]
    h = M / samples_per_module
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    # Samples are offset half a step from the box faces so no sample lies exactly
    # on the cut plane (that would produce degenerate marching-cubes triangles).
    axes = [np.arange(lo[a] - 1.5 * h, hi[a] + 2 * h, h) for a in range(3)]
    d = sheet_distance(cfg["surface"], M, *axes)
    field = np.maximum(d - thickness / 2, box_sdf(*axes, lo, hi))
    verts, faces, _, _ = measure.marching_cubes(field, 0.0, spacing=(h, h, h))
    verts += np.array([axes[0][0], axes[1][0], axes[2][0]])
    mesh = trimesh.Trimesh(verts, faces, process=True)
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    # Clipping a periodic sheet to a box can leave crumbs (~1 cm3) at the box
    # corners that belong to the neighbouring module; drop them so each print
    # file is a single solid.
    bodies = mesh.split(only_watertight=False)
    if len(bodies) > 1:
        keep = max(b.volume for b in bodies) * 1e-3
        mesh = trimesh.util.concatenate([b for b in bodies if b.volume > keep])
    mesh.fix_normals()
    return mesh


def mid_surface(cfg, samples_per_module):
    """The bare minimal surface of one module (no thickness), for previews."""
    M = cfg["module_size_m"]
    h = M / samples_per_module
    xs = np.linspace(0, M, samples_per_module + 1)
    k = 2 * math.pi / M
    X, Y, Z = np.meshgrid(xs * k, xs * k, xs * k, indexing="ij")
    verts, faces, _, _ = measure.marching_cubes(tpms(cfg["surface"], X, Y, Z), 0.0,
                                                spacing=(h, h, h))
    return trimesh.Trimesh(verts, faces, process=True)


# --------------------------------------------------------------------------- #
# Layout: which grid cells carry a module and what they do
# --------------------------------------------------------------------------- #

def grid_extent(cfg):
    M = cfg["module_size_m"]
    size = cfg["size_m"]
    for a, s in zip("xyz", size):
        if s < MIN_EXTENT_M:
            raise ValueError(f"size along {a} is {s} m; the minimum boundary is {MIN_EXTENT_M} m")
    n = [math.ceil(round(s / M, 6)) for s in size]
    return n, [ni * M for ni in n]


def entrance_cells(cfg, n):
    nx, ny, _ = n
    cut = set()
    for e in cfg.get("entrances", []):
        side = e["side"]
        if side not in SIDES:
            raise ValueError(f"entrance side must be one of {SIDES}")
        along = nx if side in ("south", "north") else ny
        w = min(e["width_modules"], along - 2)  # always keep the two corner columns
        start = (along - w) // 2
        for t in range(start, start + w):
            for k in range(min(e["height_modules"], n[2] - 1)):  # never cut the roof
                if side == "south":
                    cut.add((t, 0, k))
                elif side == "north":
                    cut.add((t, ny - 1, k))
                elif side == "west":
                    cut.add((0, t, k))
                else:
                    cut.add((nx - 1, t, k))
    return cut


def column_cells(cfg, n):
    """Interior columns on a regular grid when the roof span gets too large."""
    nx, ny, nz = n
    d = cfg["envelope_depth_modules"]
    span = cfg["max_roof_span_modules"]
    cols = set()

    def lines(count):
        clear = count - 2 * d
        if clear <= span:
            return []
        bays = math.ceil((clear + 1) / (span + 1))
        step = (clear + 1) / bays
        return [d - 1 + round(step * b) for b in range(1, bays)]

    xi, yj = lines(nx), lines(ny)
    if xi or yj:
        # A long direction alone still needs columns: put them on the centre
        # line of the short direction.
        xi = xi or [nx // 2]
        yj = yj or [ny // 2]
    for i in xi:
        for j in yj:
            for k in range(nz - d):
                cols.add((i, j, k))
    return cols


def layout(cfg):
    n, extent = grid_extent(cfg)
    nx, ny, nz = n
    d = cfg["envelope_depth_modules"]
    doors = entrance_cells(cfg, n)
    cols = column_cells(cfg, n)
    cells = {}
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                edge = min(i, nx - 1 - i, j, ny - 1 - j) < d
                roof = k >= nz - d
                if (i, j, k) in doors:
                    continue
                if roof and edge:
                    role = "roof_edge"
                elif roof:
                    role = "roof"
                elif edge:
                    corner = min(i, nx - 1 - i) < d and min(j, ny - 1 - j) < d
                    role = "corner" if corner else "wall"
                elif (i, j, k) in cols:
                    role = "column"
                else:
                    continue
                cells[(i, j, k)] = role
    return n, extent, cells


def level_thickness(cfg, k, nz):
    t = cfg["sheet_thickness_m"]
    if t.get("levels") == "uniform" or nz == 1:
        return t["base"]
    return t["base"] + (t["top"] - t["base"]) * k / (nz - 1)


# --------------------------------------------------------------------------- #
# Connectors: where the sheet touches the shared face of two modules
# --------------------------------------------------------------------------- #

def face_connectors(cfg, thickness, axis, samples=192):
    """Connector points (dowels / bolts) on a module face normal to `axis`.

    The sheet meets the shared face of two modules along thin curved lands. The
    connectors are spread along the mid-line of those lands at the configured
    spacing. The pattern is identical on every face with the same normal
    (periodicity), so it is computed once per module type and axis.
    """
    M = cfg["module_size_m"]
    spacing = cfg["connectors"]["spacing_m"]
    u = (np.arange(samples) + 0.5) * M / samples
    grids = [u, u, u]
    grids[axis] = np.array([0.0, M / samples])  # second sample only for the gradient
    d = np.take(sheet_distance(cfg["surface"], M, *grids), 0, axis=axis)

    # Mid-line of the lands, kept clear of the face border by the edge distance.
    edge = cfg["connectors"]["edge_distance_m"]
    mid = d <= min(thickness / 4, M / samples)
    uv = np.argwhere(mid) * (M / samples) + M / samples / 2
    uv = uv[np.all((uv > edge) & (uv < M - edge), axis=1)]

    chosen = []
    for q in uv[np.lexsort((uv[:, 1], uv[:, 0]))]:
        if all(np.hypot(*(q - c)) >= spacing for c in chosen):
            chosen.append(q)
    return [[round(float(a), 4), round(float(b), 4)] for a, b in chosen]


# --------------------------------------------------------------------------- #
# Print segmentation
# --------------------------------------------------------------------------- #

def split_counts(cfg):
    M = cfg["module_size_m"]
    return [max(1, math.ceil(round(M / b, 6))) for b in cfg["printer"]["build_volume_m"]]


def print_parts(cfg, thickness, res):
    M = cfg["module_size_m"]
    sx, sy, sz = split_counts(cfg)
    parts = []
    for a in range(sx):
        for b in range(sy):
            for c in range(sz):
                lo = [a * M / sx, b * M / sy, c * M / sz]
                hi = [(a + 1) * M / sx, (b + 1) * M / sy, (c + 1) * M / sz]
                parts.append(sheet_solid(cfg, thickness, lo, hi, res))
    return parts


# --------------------------------------------------------------------------- #
# Build
# --------------------------------------------------------------------------- #

def build(cfg, out_dir, write_meshes=True):
    M = cfg["module_size_m"]
    n, extent, cells = layout(cfg)
    nz = n[2]
    res = cfg["resolution"]
    out = Path(out_dir)
    (out / "modules").mkdir(parents=True, exist_ok=True)
    (out / "parts").mkdir(parents=True, exist_ok=True)

    # One module type per distinct sheet thickness (i.e. per level when graded).
    type_of_level = {}
    types = {}
    for k in sorted({k for (_, _, k) in cells}):
        t = round(level_thickness(cfg, k, nz), 4)
        if t not in type_of_level.values():
            name = f"T{len(types) + 1:02d}"
            types[name] = {"sheet_thickness_m": t, "levels": []}
        name = next(nm for nm, v in types.items() if v["sheet_thickness_m"] == t)
        type_of_level[k] = name
        types[name]["levels"].append(k)

    density = cfg["printer"]["density_kg_m3"]
    splits = split_counts(cfg)
    for name, info in types.items():
        t = info["sheet_thickness_m"]
        print(f"[{name}] sheet {t * 1000:.0f} mm ...", flush=True)
        mesh = sheet_solid(cfg, t, [0, 0, 0], [M, M, M], res["module"])
        info["volume_m3"] = round(float(mesh.volume), 5)
        info["mass_kg"] = round(float(mesh.volume) * density, 1)
        info["watertight"] = bool(mesh.is_watertight)
        info["triangles"] = int(len(mesh.faces))
        info["connectors_per_face"] = {ax: face_connectors(cfg, t, a)
                                       for a, ax in enumerate("xyz")}
        info["print_parts"] = int(np.prod(splits))
        if write_meshes:
            mesh.export(out / "modules" / f"{name}.stl")
            if np.prod(splits) > 1:
                for p, part in enumerate(print_parts(cfg, t, res["module"])):
                    part.export(out / "parts" / f"{name}_p{p + 1:02d}.stl")

    # Placement list and module-to-module joints.
    placements = []
    joints = 0
    bolts = 0
    for (i, j, k), role in sorted(cells.items(), key=lambda c: (c[0][2], c[0][1], c[0][0])):
        tname = type_of_level[k]
        placements.append({"cell": [i, j, k], "origin_m": [i * M, j * M, k * M],
                           "type": tname, "role": role})
        for a, nb in enumerate([(i + 1, j, k), (i, j + 1, k), (i, j, k + 1)]):
            if nb in cells:
                joints += 1
                # The thinner of the two sheets governs the shared contact lands.
                thinner = min(tname, type_of_level[nb[2]],
                              key=lambda nm: types[nm]["sheet_thickness_m"])
                bolts += len(types[thinner]["connectors_per_face"]["xyz"[a]])

    ground = sum(1 for (_, _, k) in cells if k == 0)
    counts = {nm: sum(1 for p in placements if p["type"] == nm) for nm in types}
    total_mass = sum(types[nm]["mass_kg"] * c for nm, c in counts.items())
    manifest = {
        "config": cfg,
        "grid_modules": n,
        "extent_m": extent,
        "module_count": len(placements),
        "roles": {r: sum(1 for p in placements if p["role"] == r)
                  for r in sorted({p["role"] for p in placements})},
        "module_types": types,
        "quantities": counts,
        "print_parts_total": len(placements) * int(np.prod(splits)),
        "module_joints": joints,
        "connectors_total": bolts,
        "ground_footings": ground,
        "total_mass_kg": round(total_mass, 1),
        "placements": placements,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

    if write_meshes:
        # Light assembly for viewing: the mid-surface of every module. (The solid
        # sheets are thinner than a coarse preview grid could resolve.)
        cell = mid_surface(cfg, res["preview"])
        pieces = [cell.copy().apply_translation(p["origin_m"]) for p in placements]
        trimesh.util.concatenate(pieces).export(out / "assembly_preview.stl")
    return manifest


# --------------------------------------------------------------------------- #
# Preview image
# --------------------------------------------------------------------------- #

ROLE_COLOURS = {"corner": "#3d5a80", "wall": "#98c1d9", "roof": "#ee6c4d",
                "roof_edge": "#c4553b", "column": "#293241"}


def render_preview(cfg, manifest, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    M = cfg["module_size_m"]
    nx, ny, nz = manifest["grid_modules"]
    mesh = trimesh.load(Path(out_dir) / "assembly_preview.stl")

    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    tris = mesh.vertices[mesh.faces]
    light = np.array([0.4, -0.6, 0.7]) / np.linalg.norm([0.4, -0.6, 0.7])
    shade = np.clip(np.abs(mesh.face_normals @ light), 0.15, 1.0)
    base = np.array([0.93, 0.91, 0.86])
    pc = Poly3DCollection(tris, facecolors=base * shade[:, None], linewidths=0)
    ax.add_collection3d(pc)
    ext = manifest["extent_m"]
    ax.set_xlim(0, max(ext)); ax.set_ylim(0, max(ext)); ax.set_zlim(0, max(ext))
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-58)
    ax.set_axis_off()
    ax.set_title(f"{cfg['surface']} pavilion  {ext[0]:g} x {ext[1]:g} x {ext[2]:g} m  "
                 f"({manifest['module_count']} modules)")

    # Plan of the ground level and the roof, coloured by role.
    ax2 = fig.add_subplot(1, 2, 2)
    for p in manifest["placements"]:
        i, j, k = p["cell"]
        if k not in (0, nz - 1):
            continue
        dx = 0 if k == 0 else (nx + 1) * M
        ax2.add_patch(plt.Rectangle((i * M + dx, j * M), M, M,
                                    facecolor=ROLE_COLOURS[p["role"]], edgecolor="white"))
        ax2.text(i * M + dx + M / 2, j * M + M / 2, p["type"], ha="center", va="center",
                 fontsize=7, color="white")
    for k, dx, label in [(0, 0, "ground level"), (nz - 1, (nx + 1) * M, "roof level")]:
        ax2.add_patch(plt.Rectangle((dx, 0), nx * M, ny * M, fill=False, lw=1.2))
        ax2.text(dx + nx * M / 2, ny * M + 0.4, label, ha="center")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ROLE_COLOURS.values()]
    ax2.legend(handles, ROLE_COLOURS.keys(), loc="lower center",
               bbox_to_anchor=(0.5, -0.22), ncol=5, frameon=False)
    ax2.set_xlim(-0.5, (2 * nx + 1) * M + 0.5); ax2.set_ylim(-0.5, ny * M + 1.2)
    ax2.set_aspect("equal")
    ax2.set_title("module layout (m, south at the bottom)")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "preview.png", dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #

def merge(base, over):
    out = dict(base)
    for k, v in over.items():
        out[k] = merge(base[k], v) if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", help="JSON config (merged over the defaults)")
    ap.add_argument("--size", nargs=3, type=float, metavar=("X", "Y", "Z"), help="extent in m")
    ap.add_argument("--surface", choices=["gyroid", "schwarz_p", "diamond", "neovius"])
    ap.add_argument("--out", help="output directory")
    ap.add_argument("--manifest-only", action="store_true", help="skip meshes, write manifest")
    args = ap.parse_args()

    cfg = DEFAULT_CONFIG
    if args.config:
        cfg = merge(cfg, json.loads(Path(args.config).read_text()))
    if args.size:
        cfg = merge(cfg, {"size_m": args.size})
    if args.surface:
        cfg = merge(cfg, {"surface": args.surface})
    out = args.out or Path(__file__).parent / cfg["output_dir"]

    m = build(cfg, out, write_meshes=not args.manifest_only)
    if not args.manifest_only:
        render_preview(cfg, m, out)
    print(f"extent {m['extent_m']} m, grid {m['grid_modules']}, {m['module_count']} modules "
          f"({m['quantities']}), {m['print_parts_total']} print parts, "
          f"{m['connectors_total']} connectors, ~{m['total_mass_kg']} kg")


if __name__ == "__main__":
    main()
