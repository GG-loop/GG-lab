#!/usr/bin/env python3
"""
Closed rooms in the open sculpture, made with tensioned textile membranes.

The gyroid splits space into two labyrinths, each a chain of chambers joined by
narrow necks. A room is one or more chambers; it is closed by stretching a
textile membrane across every opening it has: the necks to the next chamber
("diaphragms") and the gaps between the printed rims where it meets the
outside. One opening at ground level stays open as the door.

Each membrane is form-found as a minimal surface (a soap film) on its boundary:
the printed rim of the sheet and, where it reaches down, the ground. That is
the shape a uniformly prestressed fabric takes, so the textile and the printed
shell are the same mathematics at two scales.

Steps:
  1. chambers: watershed of the air on its distance to the sheet; the
     watershed lines fall on the necks, the narrowest sections
  2. room: the chamber (or adjacent group) with the most floor per m2 of
     membrane, and preferably off the climbing line
  3. openings: the room's boundary away from the sheet and the ground, split
     into patches; the tallest one touching the ground is the door
  4. form finding: each patch relaxed to a minimal surface with its edge
     snapped onto the sheet (or the ground)
  5. checks: area, double curvature (anticlastic share, flat zones that
     flutter), ponding, textile force under wind and snow, edge pull on the rims

    python membranes.py            # writes out_sculpture/membranes.stl and membranes.json
"""

import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from skimage import measure
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
import trimesh

import organic_pavilion as op
from pavilion import merge

HERE = Path(__file__).parent

TEXTILE = {
    "name": "PVC-coated polyester, type II",
    "strength_kN_m": 84.0,       # warp, characteristic
    "safety_factor": 4.0,        # membrane practice for wind / snow (long-term)
    "prestress_kN_m": 1.5,       # both directions
    "weight_kg_m2": 0.9,
}


def chambers(F, axes, env, min_peak_m=1.0, min_depth_m=0.8):
    """Watershed the air (outside the sheet, inside the envelope, above ground) into chambers."""
    h = float(axes[0][1] - axes[0][0])
    Z = axes[2][None, None, :]
    air = (F > 0) & (env < 0) & (Z > 0)
    dist = ndimage.distance_transform_edt(air) * h
    peaks = peak_local_max(dist, min_distance=max(1, int(min_peak_m / h)), threshold_abs=min_depth_m,
                           labels=air.astype(int))
    markers = np.zeros(F.shape, np.int32)
    for i, p in enumerate(peaks):
        markers[tuple(p)] = i + 1
    return watershed(-dist, markers, mask=air), air


def choose_room(lab, F, axes, env, route_pts, max_group=3, min_floor=6.0):
    """The chamber group with the most floor (2.1 m headroom) per m2 of opening."""
    h = float(axes[0][1] - axes[0][0])
    zs = axes[2]
    Z = zs[None, None, :]
    outside = (env >= 0) & (Z > 0)
    k0, k1 = np.searchsorted(zs, 0.1), np.searchsorted(zs, 2.1)
    n = int(lab.max())
    adj = set()
    for c in range(1, n + 1):
        d = ndimage.binary_dilation(lab == c) & (lab > 0) & (lab != c)
        for o in np.unique(lab[d]):
            adj.add((min(c, int(o)), max(c, int(o))))
    ix = np.clip(np.round((route_pts - [a[0] for a in axes]) / h).astype(int), 0, np.array(F.shape) - 1)
    ix[:, 2] = np.minimum(ix[:, 2] + int(0.3 / h), F.shape[2] - 1)
    route_lab = lab[ix[:, 0], ix[:, 1], ix[:, 2]]
    best = None
    for size in range(1, max_group + 1):
        for combo in itertools.combinations(range(1, n + 1), size):
            if size > 1 and not all(any((min(a, b), max(a, b)) in adj for b in combo if b != a) for a in combo):
                continue
            m = np.isin(lab, combo)
            floor = m[:, :, k0:k1 + 1].all(axis=2).sum() * h * h
            if floor < min_floor:
                continue
            sh = ndimage.binary_dilation(m) & ~m
            opening = ((sh & (lab > 0)).sum() + (sh & outside).sum()) * h * h
            score = floor / max(opening, 1e-6)
            if best is None or score > best[0]:
                best = (score, combo, floor, opening, float(np.isin(route_lab, combo).mean()))
    if best is None:
        return None
    score, combo, floor, opening, on_route = best
    return {"chambers": list(combo), "floor_m2": round(floor, 1), "opening_m2": round(opening, 1),
            "floor_per_membrane": round(score, 3), "route_share_inside": round(on_route, 2),
            "mask": np.isin(lab, combo)}


def opening_patches(room, F, axes, min_membrane=0.5):
    """Room boundary away from the sheet and the ground. Large pieces become
    membranes; the thin slivers that run a few centimetres off the sheet are the
    membrane hems, clamped onto the printed rims."""
    h = float(axes[0][1] - axes[0][0])
    o = np.array([a[0] for a in axes])
    ind = ndimage.gaussian_filter(room.astype(np.float32), 1.2)
    verts, faces, _, _ = measure.marching_cubes(ind, 0.5, spacing=(h, h, h))
    blob = trimesh.Trimesh(verts + o, faces, process=True)
    c = blob.triangles_center
    ix = np.clip(np.round((c - o) / h).astype(int), 0, np.array(F.shape) - 1)
    keep = (F[ix[:, 0], ix[:, 1], ix[:, 2]] > 0.01) & (c[:, 2] > 0.06)   # not inside the sheet, not the floor
    sub = blob.submesh([np.where(keep)[0]], append=True)
    parts = sub.split(only_watertight=False)
    big = [p for p in parts if p.area >= min_membrane]
    hems = [p for p in parts if p.area < min_membrane]
    return big, hems, blob


def sealed(room, F, axes, env, walls, core_m=0.6):
    """True when no air path leads from the room's core (cells more than core_m
    from its boundary) to the outside, with the given membranes as walls."""
    h = float(axes[0][1] - axes[0][0])
    o = np.array([a[0] for a in axes])
    shape = np.array(F.shape)
    wall = np.zeros(F.shape, bool)
    for m in walls:
        pts, _ = trimesh.sample.sample_surface_even(m, int(m.area * 3000) + 100)
        ix = np.clip(np.round((pts - o) / h).astype(int), 0, shape - 1)
        wall[tuple(ix.T)] = True
    wall = ndimage.binary_dilation(wall)
    solid = ndimage.binary_dilation(F <= 0)                 # no digital tunnels through the thin sheet
    Z = axes[2][None, None, :]
    free = ~solid & ~wall & (Z > 0.1)                       # membranes are anchored into a ground beam
    L, _ = ndimage.label(free)
    out_ids = np.unique(L[free & (env > 0.5)])
    out_ids = out_ids[out_ids > 0]
    core = ndimage.binary_erosion(room, iterations=max(1, int(core_m / h)))
    return not bool((np.isin(L, out_ids) & core).any())


def form_find(patch, sheet, ground_tol=0.3, sheet_snap=0.35, high_point=None):
    """Minimal-surface-like membrane on the patch boundary: boundary vertices are
    snapped to the sheet (or the ground) and fixed; interior vertices solve the
    discrete Laplace equation (a uniform-weight harmonic surface), then a few
    cotangent-weighted smoothing passes approach the minimal surface."""
    m = patch.copy()
    m.merge_vertices()
    v = m.vertices.copy()
    edges = m.edges_sorted
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    bnd = np.unique(uniq[counts == 1].ravel())
    tree = cKDTree(sheet.vertices)
    on_ground = bnd[v[bnd, 2] < ground_tol]
    on_sheet = np.setdiff1d(bnd, on_ground)
    d, k = tree.query(v[on_sheet])
    snap = d < sheet_snap
    # clamp the edge into the printed rim: 3 cm past the sheet surface (a keder or clamping detail)
    target = sheet.vertices[k[snap]]
    inward = -sheet.vertex_normals[k[snap]]
    v[on_sheet[snap]] = target + inward * 0.03
    v[on_ground, 2] = 0.0
    fixed = np.zeros(len(v), bool)
    fixed[bnd] = True
    if high_point is not None:
        # a strut or ring pushes the middle of the membrane out: a cone-like, doubly curved shape
        inner = np.setdiff1d(np.arange(len(v)), bnd)
        c = v[inner].mean(0)
        k = inner[np.argmin(np.linalg.norm(v[inner] - c, axis=1))]
        ring = inner[np.linalg.norm(v[inner] - v[k], axis=1) < 0.15]
        v[ring] += high_point
        fixed[ring] = True
    # graph Laplacian solve for the interior
    e = m.edges_unique
    n = len(v)
    W = sparse.coo_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n)).tocsr()
    L = sparse.diags(np.asarray(W.sum(1)).ravel()) - W
    free = np.where(~fixed)[0]
    if len(free):
        A = L[free][:, free].tocsc()
        B = L[free][:, np.where(fixed)[0]]
        for ax in range(3):
            v[free, ax] = spsolve(A, -B @ v[fixed, ax])
    m.vertices = v
    # mean-curvature smoothing (cotangent weights) with the boundary held
    for _ in range(30):
        lap = cot_laplacian_step(m)
        v = m.vertices.copy()
        v[free] += 0.5 * lap[free]
        m.vertices = v
    m.fix_normals()
    return m, int(len(bnd)), float(np.ptp(np.asarray(m.vertices)[bnd], axis=0).max()) if len(bnd) else 0.0


def cot_laplacian_step(m):
    """Cotangent Laplacian displacement (towards the mean-curvature flow)."""
    v, f = m.vertices, m.faces
    out = np.zeros_like(v)
    wsum = np.zeros(len(v))
    for i in range(3):
        a, b, c = f[:, i], f[:, (i + 1) % 3], f[:, (i + 2) % 3]
        u, w = v[b] - v[a], v[c] - v[a]
        cot = (u * w).sum(1) / (np.linalg.norm(np.cross(u, w), axis=1) + 1e-12)
        cot = np.clip(cot, 0, 10)
        for p, q in ((b, c), (c, b)):
            np.add.at(out, p, cot[:, None] * (v[q] - v[p]))
            np.add.at(wsum, p, cot)
    return out / np.maximum(wsum, 1e-9)[:, None]


def curvatures(m, smooth=3):
    """Per-vertex Gaussian curvature (angle defect / barycentric area) and mean
    curvature (half the cotangent Laplacian), lightly averaged over neighbours."""
    v, f = m.vertices, m.faces
    area = np.zeros(len(v))
    np.add.at(area, f.ravel(), np.repeat(m.area_faces / 3, 3))
    K = trimesh.curvature.vertex_defects(m) / np.maximum(area, 1e-9)
    lap = np.zeros_like(v)
    for i in range(3):
        a, b, c = f[:, i], f[:, (i + 1) % 3], f[:, (i + 2) % 3]
        u, w = v[b] - v[a], v[c] - v[a]
        cot = (u * w).sum(1) / (np.linalg.norm(np.cross(u, w), axis=1) + 1e-12)
        np.add.at(lap, b, cot[:, None] * (v[c] - v[b]))
        np.add.at(lap, c, cot[:, None] * (v[b] - v[c]))
    H = np.linalg.norm(lap, axis=1) / (4 * np.maximum(area, 1e-9))
    nb = m.vertex_neighbors
    for _ in range(smooth):
        K = np.array([np.mean(np.r_[K[i], K[nb[i]]]) for i in range(len(v))])
        H = np.array([np.mean(np.r_[H[i], H[nb[i]]]) for i in range(len(v))])
    return K, H


def membrane_checks(m, textile=TEXTILE, wind=0.8, snow=1.2):
    """Curvature, ponding and textile force (screening)."""
    pts = m.vertices
    K, H = curvatures(m)
    edges = m.edges_sorted
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    bnd = np.unique(uniq[counts == 1].ravel())
    inner = np.setdiff1d(np.arange(len(pts)), bnd)
    Ki = K[inner] if len(inner) else K
    R = 1 / np.sqrt(np.abs(Ki) + 1e-9)
    anticlastic = float((Ki < -1e-3).mean())
    flat = float((np.abs(Ki) < 1 / 12 ** 2).mean())
    # ponding: interior points facing up, nearly level and lower than their neighbours
    n = m.vertex_normals
    level = np.abs(n[:, 2]) > 0.95
    nbr = m.vertex_neighbors
    lowest = np.array([all(pts[i, 2] <= pts[j, 2] + 1e-4 for j in nbr[i]) if nbr[i] else False for i in range(len(pts))])
    pond = int((level & lowest)[inner].sum()) if len(inner) else 0
    # textile force: p * R (single-curvature bound, conservative) at the 90th-percentile radius
    r90 = float(np.percentile(np.minimum(R, 50), 90))
    up = np.abs(n[:, 2]).mean()
    p_design = 1.5 * max(wind, snow * up)       # kN/m2
    T = textile["prestress_kN_m"] + p_design * r90
    allow = textile["strength_kN_m"] / textile["safety_factor"]
    return {
        "area_m2": round(float(m.area), 2),
        "span_m": round(float(np.linalg.norm(m.bounds[1] - m.bounds[0])), 2),
        "anticlastic_share": round(anticlastic, 2),
        "flat_share": round(flat, 2),
        "mean_curvature_abs_median": round(float(np.median(np.abs(H[inner] if len(inner) else H))), 3),
        "radius_p90_m": round(r90, 2),
        "ponding_points": pond,
        "design_pressure_kN_m2": round(p_design, 2),
        "textile_force_kN_m": round(T, 2),
        "textile_allowable_kN_m": round(allow, 2),
        "utilisation": round(T / allow, 2),
    }


def cut_door(m, width=1.2, height=2.2):
    """Open a zipped door in a membrane: the faces within a 1.2 x 2.2 m frame at the
    middle of its ground edge are removed (the frame is a printed or steel door ring)."""
    edges = m.edges_sorted
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    bnd = np.unique(uniq[counts == 1].ravel())
    low = bnd[m.vertices[bnd, 2] < 0.1]
    if len(low) == 0:
        return m, None
    g = m.vertices[low]
    centre = g[np.argmin(np.linalg.norm(g - g.mean(0), axis=1))]
    c = m.triangles_center
    keep = ~((np.hypot(c[:, 0] - centre[0], c[:, 1] - centre[1]) < width / 2) & (c[:, 2] < height))
    out = m.submesh([np.where(keep)[0]], append=True)
    return out, {"width_m": width, "clear_height_m": height, "at_m": [round(float(x), 2) for x in centre],
                 "type": "zipped door in the membrane, with a frame"}


def cut_hatches(membranes, route, radius=0.6, near=0.3):
    """Where the climbing line passes through a membrane, open a round hatch
    (a zipped cover closes it) so climbers can leave the room upwards."""
    hatches = []
    out = []
    for k, m in enumerate(membranes):
        d, idx = cKDTree(m.vertices).query(route)
        hit = np.where(d < near)[0]
        if len(hit) == 0:
            out.append(m)
            continue
        c = route[hit[len(hit) // 2]]
        keep = np.linalg.norm(m.triangles_center - c, axis=1) > radius
        out.append(m.submesh([np.where(keep)[0]], append=True))
        same = [hh for hh in hatches if np.linalg.norm(np.array(hh["at_m"]) - c) < 1.0]
        if same:                     # one passage through two neighbouring membranes
            same[0]["membranes"].append(k)
        else:
            hatches.append({"membranes": [k], "at_m": [round(float(x), 2) for x in c], "diameter_m": 2 * radius,
                            "type": "round climbing hatch with a zipped cover"})
    return out, hatches


def boundary_forces(membranes, textile=TEXTILE, wind=0.8):
    """Line loads the membranes put on the printed rims (for the structural model):
    prestress pulling along the membrane into its span at every boundary vertex,
    plus each membrane's wind resultant shared over its boundary."""
    forces = []
    for m in membranes:
        edges = m.edges_sorted
        uniq, counts = np.unique(edges, axis=0, return_counts=True)
        be = uniq[counts == 1]
        if len(be) == 0:
            continue
        length = np.linalg.norm(m.vertices[be[:, 1]] - m.vertices[be[:, 0]], axis=1)
        share = np.zeros(len(m.vertices))
        np.add.at(share, be[:, 0], length / 2)
        np.add.at(share, be[:, 1], length / 2)
        centre = m.vertices.mean(0)
        nrm = m.vertex_normals
        wind_res = 1.5 * wind * 1e3 * float((m.area_faces * np.abs(m.face_normals[:, 0])).sum())  # N along x
        total_len = share.sum()
        for i in np.where(share > 0)[0]:
            p = m.vertices[i]
            d = centre - p
            d -= nrm[i] * (d @ nrm[i])                       # in the membrane's tangent plane
            d /= np.linalg.norm(d) + 1e-9
            f = d * textile["prestress_kN_m"] * 1e3 * share[i]
            f[0] += wind_res * share[i] / total_len
            forces.append((p, f))
    return forces


def build(out=HERE / "out_sculpture", voxel=0.05, log=print):
    import analysis as an
    cfg = merge(op.DEFAULT_CONFIG, op.ART)
    report = json.loads((out / "report.json").read_text())
    fine = merge(cfg, {"voxel_m": voxel})
    log("field ...")
    F, axes, ex = op.field(fine, extras=True)
    sheet = op.to_mesh(F, axes, voxel)
    lab, air = chambers(F, axes, ex["env"])
    log(f"chambers: {lab.max()}")
    route = np.array(report["climb"]["route"])
    room = choose_room(lab, F, axes, ex["env"], route)
    log(f"room: chambers {room['chambers']}, floor {room['floor_m2']} m2, openings {room['opening_m2']} m2")
    patches, hems, blob = opening_patches(room["mask"], F, axes)
    # the door goes into the smallest opening that reaches the ground and is tall enough
    ground = [i for i, p in enumerate(patches) if p.bounds[0][2] < 0.25 and p.bounds[1][2] > 2.3]
    door = min(ground, key=lambda i: patches[i].area) if ground else None
    room_c = np.argwhere(room["mask"]).mean(0) * voxel + np.array([a[0] for a in axes])
    membranes, info, door_info = [], [], None
    for i, p in enumerate(patches):
        m, nb, span = form_find(p, sheet)
        chk = membrane_checks(m)
        high = None
        if chk["flat_share"] > 0.25 or chk["utilisation"] > 1:
            # too flat to carry wind: push its middle out of the room by span / 8
            n = m.vertex_normals.mean(0)
            n /= np.linalg.norm(n) + 1e-9
            if n @ (m.vertices.mean(0) - room_c) < 0:
                n = -n
            rise = chk["span_m"] / 8
            m, nb, span = form_find(p, sheet, high_point=n * rise)
            chk = membrane_checks(m)
            high = {"rise_m": round(rise, 2), "at_m": [round(float(c), 2) for c in (p.vertices.mean(0) + n * rise)]}
        membranes.append(m)
        info.append({"id": f"M{len(membranes)}", "centre_m": [round(float(c), 2) for c in m.vertices.mean(0)],
                     "high_point": high, "has_door": i == door, **chk})
    is_sealed = sealed(room["mask"], F, axes, ex["env"], membranes + hems)
    log(f"sealed with the door closed: {is_sealed}")
    if door is not None:
        membranes[door], door_info = cut_door(membranes[door])
    # the climbing line (computed on the untreaded sheet, as in analysis.py)
    plain = merge(cfg, {"terraces": {"enabled": False}})
    F0, ax0 = op.field(plain)
    clim = an.climb(plain, F0, ax0, op.to_mesh(F0, ax0, plain["voxel_m"]), cfg["terraces"]["riser_m"], log=log)
    route = np.array(clim["route"])
    membranes, hatches = cut_hatches(membranes, route)
    for hch in hatches:
        for k in hch["membranes"]:
            info[k]["hatch"] = True
    log(f"hatches: {hatches}")
    log(f"membranes: {len(membranes)}, total {sum(x['area_m2'] for x in info):.1f} m2; door {door_info}")
    trimesh.util.concatenate(membranes + hems).export(out / "membranes.stl")
    h = voxel
    zs = axes[2]
    k0, k1 = np.searchsorted(zs, 0.1), np.searchsorted(zs, 2.1)
    roomm = room["mask"]
    zz = np.nonzero(roomm)[2]
    result = {
        "textile": TEXTILE,
        "room": {"chambers": room["chambers"], "volume_m3": round(float(roomm.sum() * h ** 3), 1),
                 "floor_m2": room["floor_m2"], "max_height_m": round(float(zs[zz].max()), 2),
                 "floor_per_membrane_m2": room["floor_per_membrane"]},
        "door": door_info,
        "sealed_with_door_closed": is_sealed,
        "hatches": hatches,
        "route_summit_m": clim["summit_m"],
        "hem_area_m2": round(float(sum(p.area for p in hems)), 1),
        "membranes": info,
        "total_area_m2": round(sum(x["area_m2"] for x in info), 1),
        "voxel_m": voxel,
    }
    np.savez_compressed(out / "room.npz", mask=roomm, origin=np.array([a[0] for a in axes]), h=voxel)
    forces = boundary_forces(membranes)
    np.savez_compressed(out / "membrane_forces.npz", p=np.array([f[0] for f in forces]), f=np.array([f[1] for f in forces]))
    (out / "membranes.json").write_text(json.dumps(result, indent=1))
    return result, membranes, roomm, axes


if __name__ == "__main__":
    r, _, _, _ = build()
    print(json.dumps({k: v for k, v in r.items() if k != "textile"}, indent=1))
