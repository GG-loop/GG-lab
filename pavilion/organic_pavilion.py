#!/usr/bin/env python3
"""
Organic pavilion: one minimal-surface cell is one pavilion ("pod").

A triply periodic minimal surface splits space into two labyrinths. Each pod
keeps the central chamber of one cell and trims everything else with a soft
envelope: no plane cuts anything except the ground. For the default Neovius
cell the 12 arms of the chamber that run towards the cell edges become the
openings:

    4 diagonal arms at mid height  -> walk-through portals at the corners
    4 upper arms on the axes        -> clerestory openings above the walls
    the 6 face caps                 -> 4 lobed walls, a dished roof, a buried floor

Mathematics with a twist: an organic layer bends the cells. Each column of
pods twists and tapers about its own axis and the whole thing sways slowly,
all fading to zero at the ground so the footprint stays planted. The sheet is
graded in thickness and thickens along free rims, which work as edge beams.

Scalable and modular: pods sit on a 3D grid (i, j, level). Neighbouring cells
of a periodic surface share their arms, so pods fuse into one form:
  * side by side they are joined by tunnels along the shared cell edges,
    whose floor hump is trimmed so they can be walked through;
  * stacked, the dished roof of one chamber is the floor of the next.

    python organic_pavilion.py                         # one pod, >= 6 x 6 x 6 m
    python organic_pavilion.py --pods 0,0,0 1,1,0      # two pods on the diagonal
    python organic_pavilion.py --pods 0,0,0 0,0,1      # a two-storey tower
    python organic_pavilion.py --surface schwarz_p -c my.json

Outputs (out_organic/ by default): pavilion.stl, parts/*.stl split to the
printer build volume, manifest.json (sizes, clearances, mass), preview.png.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from skimage import measure
import trimesh

from pavilion import tpms, merge

MIN_EXTENT_M = 6.0

# The climbable gyroid sculpture: one twisted gyroid cell whose up-facing skin is
# split into contour ledges (seats, steps, narrow ledges).
ART = {
    "surface": "gyroid",
    "cell_m": [8.0, 8.0, 8.0],
    "chamber_centre_m": 3.0,
    "envelope": {"radius_m": [4.8, 4.8, 4.4], "plan_exponent": 1.6},
    "organic": {"twist_deg": 90.0},
    "terraces": {"enabled": True, "riser_m": 0.18, "max_slope_deg": 58.0, "start_m": 0.0, "rim_keep_m": 0.30},
}

DEFAULT_CONFIG = {
    "surface": "neovius",               # neovius | schwarz_p | gyroid | diamond
    "cell_m": [7.2, 7.2, 10.0],         # one period in x, y, z (stretched upwards)
    "chamber_centre_m": 1.3,            # height of the ground-level chamber centre
    "pods": [[0, 0, 0]],                # grid positions (i, j, level) of the cells
    "envelope": {
        "radius_m": [3.95, 3.95, 5.1],  # semi-axes of the trimming envelope around each chamber
        "plan_exponent": 1.15,          # < 2 pulls the diagonals in: corners open into portals
        "vertical_exponent": 2.6,
        "rim_round_m": 0.35,            # fillet radius where the sheet meets the envelope
        "pod_blend_m": 0.8,             # smooth union between pod envelopes and links
        "link_radius_m": 1.9,           # half-width of the tunnels that join neighbouring pods
        "link_floor_m": 0.65            # floor hump inside ground-level tunnels is trimmed below this
    },
    "organic": {
        "twist_deg": 28.0,              # rotation of the top of a pod relative to the ground
        "taper": 0.10,                  # plan shrink at the top (fraction)
        "sway_m": 0.35,                 # amplitude of the slow lateral drift
        "seed": 7
    },
    "sheet": {
        "base_m": 0.10,                 # at the ground
        "top_m": 0.06,                  # at the crown
        "rim_extra_m": 0.05             # added along free rims (edge beam)
    },
    "printer": {
        "build_volume_m": [2.0, 2.0, 1.5],   # large-format robotic cell
        "layer_m": 0.012,
        "material": "rPETG",
        "density_kg_m3": 1270
    },
    "voxel_m": 0.04,
    "output_dir": "out_organic"
}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def smax(a, b, k):
    """Smooth maximum (rounded intersection) with blend radius k."""
    h = np.clip(0.5 - 0.5 * (a - b) / k, 0, 1)
    return a * (1 - h) + b * h + k * h * (1 - h)


def smin(a, b, k):
    return -smax(-a, -b, k)


def pods(cfg):
    return [tuple(p) + (0,) * (3 - len(p)) for p in cfg["pods"]]


def centre(cfg, pod):
    Lx, Ly, Lz = cfg["cell_m"]
    i, j, k = pod
    return np.array([i * Lx, j * Ly, cfg["chamber_centre_m"] + k * Lz])


def height_m(cfg):
    return max(centre(cfg, p)[2] for p in pods(cfg)) + cfg["envelope"]["radius_m"][2]


# --------------------------------------------------------------------------- #
# Topology: links between neighbouring pods
# --------------------------------------------------------------------------- #

def links(cfg):
    """Segments from each chamber centre to the cell edges it shares with another pod.

    In the Neovius every cell edge carries a tube of the chamber phase that
    connects the four cells around it, so two pods sharing an edge are joined
    by opening the arms from each chamber to that edge. The same rule is used
    for the other surfaces; their chambers simply meet through what is there.
    Returns (pod, start xyz, end xyz, vertical_edge) tuples.
    """
    L = np.array(cfg["cell_m"])
    P = set(pods(cfg))
    half = [-0.5, 0.5]
    edges = {}
    for p in P:
        es = set()
        for ax in range(3):  # edge parallel to axis `ax`: offsets in the other two
            o = [a for a in range(3) if a != ax]
            for s1, s2 in itertools.product(half, half):
                m = [0.0, 0.0, 0.0]
                m[o[0]], m[o[1]] = s1, s2
                es.add((p[0] + m[0], p[1] + m[1], p[2] + m[2], ax))
        edges[p] = es
    out = []
    for p, q in itertools.combinations(sorted(P), 2):
        for e in edges[p] & edges[q]:
            mid = np.array(e[:3]) * L + [0, 0, cfg["chamber_centre_m"]]
            for r in (p, q):
                out.append((r, centre(cfg, r), mid, e[3] == 2))
    return out


# --------------------------------------------------------------------------- #
# Field
# --------------------------------------------------------------------------- #

def warp(cfg, X, Y, Z):
    """Organic deformation, zero at the ground and growing with height.

    Every column of pods twists and tapers about its own axis; the per-column
    deformations are blended with Gaussian weights so clusters stay continuous.
    """
    o = cfg["organic"]
    Lx, Ly, Lz = cfg["cell_m"]
    # the twist is per storey: a tower keeps turning as it rises
    s = np.clip(Z / (Lz * 0.65), 0, None)
    ease = np.where(s < 1, s * s * (3 - 2 * s), s)
    a = -math.radians(o["twist_deg"]) * ease   # the sampled point turns opposite to the form
    ca, sa = np.cos(a), np.sin(a)
    sz = np.clip(Z / height_m(cfg), 0, 1)
    scale = 1.0 / (1.0 - o["taper"] * sz * sz * (3 - 2 * sz))
    sigma2 = (0.6 * max(Lx, Ly)) ** 2

    dx = np.zeros_like(X); dy = np.zeros_like(Y); wsum = np.zeros_like(X)
    for i, j in sorted({(p[0], p[1]) for p in pods(cfg)}):
        x, y = X - i * Lx, Y - j * Ly
        w = np.exp(-(x * x + y * y) / sigma2) + 1e-12
        dx += w * ((x * ca - y * sa) * scale - x)
        dy += w * ((x * sa + y * ca) * scale - y)
        wsum += w
    xr, yr = X + dx / wsum, Y + dy / wsum

    ph = sway_phases(o["seed"])
    k = 2 * math.pi / (1.6 * max(cfg["cell_m"]))
    A = o["sway_m"] * np.clip(Z / (Lz * 0.65), 0, 1)
    xr = xr + A * np.sin(k * Z + ph[0]) * np.cos(0.7 * k * Y + ph[1])
    yr = yr + A * np.sin(0.9 * k * Z + ph[2]) * np.cos(0.6 * k * X + ph[3])
    zr = Z + 0.4 * A * np.sin(0.8 * k * X + ph[4]) * np.sin(0.8 * k * Y + ph[5])
    return xr, yr, zr


def sway_phases(seed):
    """Six phases from the seed (a plain hash, so the web viewer can reproduce it)."""
    return [((math.sin(seed * 12.9898 + i * 78.233) * 43758.5453) % 1.0) * 2 * math.pi for i in range(6)]


def capsule(X, Y, Z, a, b, r, zsquash):
    v = b - a
    t = np.clip(((X - a[0]) * v[0] + (Y - a[1]) * v[1] + (Z - a[2]) * v[2]) / (v @ v), 0, 1)
    px, py, pz = X - a[0] - t * v[0], Y - a[1] - t * v[1], (Z - a[2] - t * v[2]) * zsquash
    return np.sqrt(px * px + py * py + pz * pz) - r


def envelope_sdf(cfg, X, Y, Z):
    e = cfg["envelope"]
    Rx, Ry, Rz = e["radius_m"]
    n, m = e["plan_exponent"], e["vertical_exponent"]
    out = None
    for p in pods(cfg):
        c = centre(cfg, p)
        plan = ((np.abs(X - c[0]) / Rx) ** n + (np.abs(Y - c[1]) / Ry) ** n) ** (1 / n)
        g = (plan ** m + (np.abs(Z - c[2]) / Rz) ** m) ** (1 / m)
        d = (g - 1) * min(Rx, Ry)
        out = d if out is None else smin(out, d, e["pod_blend_m"])
    r = e["link_radius_m"]
    for _, a, b, vertical in links(cfg):
        out = smin(out, capsule(X, Y, Z, a, b, r, 1.4 * r / Rz if vertical else 1.0), e["pod_blend_m"])
    return out


def tunnel_floor_sdf(cfg, X, Y, Z):
    """Region trimmed out of the floor of each ground-level tunnel (negative inside).

    Between two chambers the edge tube has a hump of sheet on its floor. It is
    removed below link_floor_m over the outer part of each link, where the hump
    sits, so the tunnel is walkable without touching the pod walls.
    """
    e = cfg["envelope"]
    out = np.full(X.shape, 1e3, dtype=np.float32)
    for pod, a, b, vertical in links(cfg):
        if not vertical or pod[2] != 0:
            continue
        vx, vy = b[0] - a[0], b[1] - a[1]
        t = np.clip(((X - a[0]) * vx + (Y - a[1]) * vy) / (vx * vx + vy * vy), 0.6, 1.0)
        d_plan = np.sqrt((X - a[0] - t * vx) ** 2 + (Y - a[1] - t * vy) ** 2) - 0.65 * e["link_radius_m"]
        out = np.minimum(out, np.maximum(d_plan, Z - e["link_floor_m"]))
    return out


def bounds(cfg):
    e = cfg["envelope"]
    C = np.array([centre(cfg, p) for p in pods(cfg)])
    pad = 0.6 + cfg["organic"]["sway_m"]
    lo = [C[:, 0].min() - e["radius_m"][0] - pad, C[:, 1].min() - e["radius_m"][1] - pad, 0.0]
    hi = [C[:, 0].max() + e["radius_m"][0] + pad, C[:, 1].max() + e["radius_m"][1] + pad, height_m(cfg) + 0.6]
    # tunnels between pods reach out to the shared cell edges
    r = e["link_radius_m"] + pad
    for _, _, b, _ in links(cfg):
        lo = [min(lo[0], b[0] - r), min(lo[1], b[1] - r), 0.0]
        hi = [max(hi[0], b[0] + r), max(hi[1], b[1] + r), max(hi[2], b[2] + r)]
    return np.array(lo), np.array(hi)


def field(cfg, extras=False):
    """Signed field on a regular grid: negative inside the printed material.

    With extras=True also returns the raw surface value f (for curvature), the
    envelope distance and the local sheet thickness, used by analysis.py."""
    h = cfg["voxel_m"]
    lo, hi = bounds(cfg)
    axes = [np.arange(lo[a] - 1.5 * h, hi[a] + 2 * h, h, dtype=np.float32) for a in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    Xw, Yw, Zw = warp(cfg, X, Y, Z)

    Lx, Ly, Lz = cfg["cell_m"]
    zc = cfg["chamber_centre_m"]
    # chamber centres sit on cell centres: shift by half a period
    f = tpms(cfg["surface"], Xw * 2 * math.pi / Lx + math.pi, Yw * 2 * math.pi / Ly + math.pi,
             (Zw - zc) * 2 * math.pi / Lz + math.pi).astype(np.float32)
    g = np.gradient(f, h)
    dist = np.abs(f) / np.maximum(np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2), 1e-6)
    del g

    env = envelope_sdf(cfg, Xw, Yw, Zw)
    s = cfg["sheet"]
    t = s["base_m"] + (s["top_m"] - s["base_m"]) * np.clip(Z / height_m(cfg), 0, 1)
    t = t + s["rim_extra_m"] * np.exp(-(env / 0.35) ** 2)
    F = smax(dist - t / 2, env, cfg["envelope"]["rim_round_m"])
    F = np.maximum(F, -Z)                                # flat cut at the ground only
    F = np.maximum(F, -tunnel_floor_sdf(cfg, X, Y, Z))   # walk-through floor between pods
    if cfg.get("terraces", {}).get("enabled"):
        F = terrace(F, axes, cfg["terraces"], env)
        F = np.maximum(F, -Z)
    if extras:
        return F, axes, {"f": f, "env": env, "thickness": t}
    return F, axes


def terrace(F, axes, t, env=None):
    """Split the up-facing skin of the sheet into horizontal ledges.

    Every column whose top surface faces upwards (slope below max_slope_deg) is
    built up to the next multiple of the riser height, so the treads follow the
    contour lines of the surface: deep seats where it is shallow, steps where it
    climbs, narrow ledges where it is steep. Steeper parts stay smooth.
    """
    h = float(axes[0][1] - axes[0][0])
    R = t["riser_m"]
    zs = axes[2]
    solid = F <= 0
    top = np.zeros_like(solid)
    top[:, :, :-1] = solid[:, :, :-1] & ~solid[:, :, 1:]
    gx, gy, gz = np.gradient(F, h)
    nz = gz / np.maximum(np.sqrt(gx * gx + gy * gy + gz * gz), 1e-9)
    del gx, gy, gz
    # top voxels whose surface normal (the gradient, pointing out of the material) faces up enough
    sel = top & (nz > math.cos(math.radians(t["max_slope_deg"])))
    if env is not None:
        sel &= env < -t.get("rim_keep_m", 0.12)        # free rims stay smooth lips
    i, j, k = np.nonzero(sel)
    f0, f1 = F[i, j, k], F[i, j, k + 1]
    ztop = zs[k] + f0 / (f0 - f1) * h            # exact height of the skin in this column
    zq = np.ceil((ztop - 0.01) / R) * R           # the ledge level just above it
    zq = np.where(zq < t.get("start_m", 0.0) + 1e-6, ztop, zq)
    F = F.copy()
    for m in range(1, int(math.ceil(R / h)) + 2):
        kk = k + m
        ok = kk < F.shape[2]
        z = zs[np.minimum(kk, F.shape[2] - 1)]
        fill = ok & (z <= zq + h) & (zq > ztop + 1e-4)
        ii, jj, kk2 = i[fill], j[fill], kk[fill]
        F[ii, jj, kk2] = np.minimum(F[ii, jj, kk2], (z[fill] - zq[fill]).astype(F.dtype))
    return F


def to_mesh(F, axes, h):
    verts, faces, _, _ = measure.marching_cubes(F, 0.0, spacing=(h, h, h))
    verts += np.array([axes[0][0], axes[1][0], axes[2][0]])
    mesh = trimesh.Trimesh(verts, faces, process=True)
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    bodies = mesh.split(only_watertight=False)
    if len(bodies) > 1:
        keep = max(b.volume for b in bodies) * 1e-3
        mesh = trimesh.util.concatenate([b for b in bodies if b.volume > keep])
    mesh.fix_normals()
    return mesh


# --------------------------------------------------------------------------- #
# Checks: can people use it?
# --------------------------------------------------------------------------- #

def _column(F, axes, x, y):
    xs, ys, zs = axes
    return F[np.abs(xs - x).argmin(), np.abs(ys - y).argmin(), :] <= 0, zs


def chamber_clearance(F, axes, cfg):
    """Floor level and clear height of every chamber, on its axis."""
    out = {}
    for p in pods(cfg):
        c = centre(cfg, p)
        solid, zs = _column(F, axes, c[0], c[1])
        k = np.abs(zs - c[2]).argmin()
        lo = k
        while lo > 0 and not solid[lo] and zs[lo] > 0:
            lo -= 1
        hi = k
        while hi < len(zs) - 1 and not solid[hi]:
            hi += 1
        floor = max(0.0, float(zs[lo]))
        roof = float(zs[hi]) if solid[hi] else float("inf")
        out[",".join(map(str, p))] = {"floor_m": round(floor, 2),
                                       "clear_height_m": round(roof - floor, 2) if roof < 1e9 else "open"}
    return out


def portal_heights(F, axes, cfg):
    """Clear height of the corner portals of ground-level pods (worst point
    along each diagonal, from the chamber to the outside)."""
    Lx = cfg["cell_m"][0]
    out = {}
    for p in pods(cfg):
        if p[2] != 0:
            continue
        c = centre(cfg, p)
        for name, (sx, sy) in {"NE": (1, 1), "NW": (-1, 1), "SW": (-1, -1), "SE": (1, -1)}.items():
            best = None
            for r in np.arange(0.5, 1.2 * Lx, cfg["voxel_m"]):
                x, y = c[0] + sx * r / math.sqrt(2), c[1] + sy * r / math.sqrt(2)
                xs, ys, zs = axes
                if not (xs[0] < x < xs[-1] and ys[0] < y < ys[-1]):
                    break
                solid, zs = _column(F, axes, x, y)
                if not solid.any():
                    break                                  # reached the outside
                if solid[np.abs(zs - 0.05).argmin()]:
                    best = 0.0                             # a wall stands here
                    break
                above = np.where(solid & (zs > 0.1))[0]
                clear = float(zs[above[0]])
                best = clear if best is None else min(best, clear)
            out[f"{','.join(map(str, p))}-{name}"] = round(best or 0.0, 2)
    return out


# --------------------------------------------------------------------------- #
# Print segmentation
# --------------------------------------------------------------------------- #

def segments(F, axes, cfg):
    h = cfg["voxel_m"]
    xs, ys, zs = axes
    lo = np.array([xs[0], ys[0], 0.0])
    hi = np.array([xs[-1], ys[-1], zs[-1]])
    b = np.array(cfg["printer"]["build_volume_m"])
    n = np.ceil((hi - lo) / b).astype(int)
    parts = []
    for a, c, d in itertools.product(range(n[0]), range(n[1]), range(n[2])):
        plo = lo + b * [a, c, d]
        phi = np.minimum(plo + b, hi)
        sl = [slice(max(0, int((plo[k] - axes[k][0]) / h) - 2), int((phi[k] - axes[k][0]) / h) + 3)
              for k in range(3)]
        sub = F[sl[0], sl[1], sl[2]]
        if sub.size == 0 or sub.min() > 0:
            continue
        sx, sy, sz = xs[sl[0]], ys[sl[1]], zs[sl[2]]
        Xs, Ys, Zs = np.meshgrid(sx, sy, sz, indexing="ij")
        q = np.stack([np.maximum(plo[0] - Xs, Xs - phi[0]), np.maximum(plo[1] - Ys, Ys - phi[1]),
                      np.maximum(plo[2] - Zs, Zs - phi[2])]).max(axis=0)
        G = np.maximum(sub, q)
        if G.min() > 0:
            continue
        m = to_mesh(G, (sx, sy, sz), h)
        if m.volume < 1e-4:
            continue
        parts.append((f"P{a:02d}{c:02d}{d:02d}", m))
    return parts, n.tolist()


# --------------------------------------------------------------------------- #

def build(cfg, out_dir, write_parts=True):
    out = Path(out_dir)
    (out / "parts").mkdir(parents=True, exist_ok=True)
    F, axes = field(cfg)
    mesh = to_mesh(F, axes, cfg["voxel_m"])
    ext = mesh.bounds[1] - mesh.bounds[0]
    if ext.min() < MIN_EXTENT_M:
        raise ValueError(f"built pavilion is {ext.round(2).tolist()} m; every side must reach "
                         f"{MIN_EXTENT_M} m - enlarge cell_m / envelope radius_m")
    mesh.export(out / "pavilion.stl")
    density = cfg["printer"]["density_kg_m3"]

    part_info, grid = [], None
    for p in out.glob("parts/*.stl"):
        p.unlink()
    if write_parts:
        parts, grid = segments(F, axes, cfg)
        for name, m in parts:
            m.export(out / "parts" / f"{name}.stl")
            part_info.append({"name": name, "extent_m": (m.bounds[1] - m.bounds[0]).round(2).tolist(),
                              "mass_kg": round(m.volume * density, 1)})

    b = mesh.bounds
    manifest = {
        "config": cfg,
        "bounding_box_m": ext.round(2).tolist(),
        "footprint_m": [round(b[0][0], 2), round(b[0][1], 2), round(b[1][0], 2), round(b[1][1], 2)],
        "chambers": chamber_clearance(F, axes, cfg),
        "portal_clear_height_m": portal_heights(F, axes, cfg),
        "watertight": bool(mesh.is_watertight),
        "volume_m3": round(float(mesh.volume), 3),
        "surface_area_m2": round(float(mesh.area) / 2, 1),
        "mass_kg": round(float(mesh.volume) * density, 0),
        "print_grid": grid,
        "print_parts": part_info,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return mesh, manifest


def render_preview(mesh, manifest, cfg, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    try:
        m = mesh.simplify_quadric_decimation(face_count=min(len(mesh.faces), 120000))
    except ImportError:  # fast-simplification not installed: render the full mesh
        m = mesh
    fig = plt.figure(figsize=(15, 6.5))
    light = np.array([0.35, -0.55, 0.75]); light /= np.linalg.norm(light)
    b = mesh.bounds
    mid = (b[0] + b[1]) / 2
    r = (b[1] - b[0]).max() / 2
    for n, (elev, azim, title) in enumerate([(12, -55, "view from south-east"), (38, 30, "view from north-west"),
                                             (88, -90, "plan")]):
        ax = fig.add_subplot(1, 3, n + 1, projection="3d")
        tris = m.vertices[m.faces]
        sh = np.clip(m.face_normals @ light, 0, 1) * 0.75 + 0.25
        ax.add_collection3d(Poly3DCollection(tris, facecolors=np.array([0.93, 0.92, 0.89]) * sh[:, None],
                                             linewidths=0))
        ax.set_xlim(mid[0] - r, mid[0] + r); ax.set_ylim(mid[1] - r, mid[1] + r); ax.set_zlim(0, 2 * r)
        ax.set_box_aspect((1, 1, 1)); ax.view_init(elev, azim); ax.set_axis_off()
        ax.set_title(title, fontsize=10)
    bb = manifest["bounding_box_m"]
    fig.suptitle(f"{cfg['surface']} pavilion, {len(cfg['pods'])} pod(s), "
                 f"{bb[0]:.1f} x {bb[1]:.1f} x {bb[2]:.1f} m, ~{manifest['mass_kg'] / 1000:.1f} t", fontsize=12)
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "preview.png", dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config")
    ap.add_argument("--pods", nargs="+", help="cells as i,j or i,j,level (e.g. 0,0 1,0 0,0,1)")
    ap.add_argument("--surface", choices=["neovius", "schwarz_p", "gyroid", "diamond"])
    ap.add_argument("--no-parts", action="store_true", help="skip print segmentation")
    ap.add_argument("--art", action="store_true", help="the climbable terraced gyroid sculpture")
    ap.add_argument("--out")
    args = ap.parse_args()
    cfg = merge(DEFAULT_CONFIG, ART) if args.art else DEFAULT_CONFIG
    if args.config:
        cfg = merge(cfg, json.loads(Path(args.config).read_text()))
    if args.pods:
        cfg = merge(cfg, {"pods": [[int(v) for v in p.split(",")] for p in args.pods]})
    if args.surface:
        cfg = merge(cfg, {"surface": args.surface})
    out = args.out or Path(__file__).parent / cfg["output_dir"]
    mesh, m = build(cfg, out, write_parts=not args.no_parts)
    render_preview(mesh, m, cfg, out)
    print(json.dumps({k: v for k, v in m.items() if k not in ("config", "print_parts")}, indent=1))
    print(f"{len(m['print_parts'])} print parts")


if __name__ == "__main__":
    main()
