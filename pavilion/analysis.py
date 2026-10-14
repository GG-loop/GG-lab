#!/usr/bin/env python3
"""
Stairs and structural feasibility for the organic pavilion.

1. Stairs ("scale"): finds the best line on the walkable top of the sheet from
   the ground to an upper floor and fits a regular, code-compliant stair to it.
   The treads are fused into the printed sheet so stair and shell are one part.

   Rules (Italian DM 236/1989 art. 8.1.10 and common European practice; check
   the local code before building):
     riser 15-17 cm, constant within a flight      tread >= 30 cm
     2 risers + tread = 62-64 cm (Blondel)          width >= 1.20 m (public)
     <= 15 risers per flight, landings >= 1.20 m   headroom >= 2.10 m
     open sides: parapet >= 1.00 m, handrail 0.90-1.00 m

   The route search runs over the mesh: every edge on the up-facing skin of the
   sheet costs its length, weighted by how far its pitch is from the ideal
   stair pitch (atan(riser / tread) ~ 27 deg), its cross slope, the distance to
   free rims and the headroom above it.

2. Structure: a voxel finite-element model of the whole printed solid
   (8-node hexahedra, 4 cm, solved with algebraic multigrid). Load cases follow
   Eurocode / NTC 2018 in simplified form: self weight, crowd on stairs and
   floors, snow on up-facing areas, wind as horizontal pressure. Checked:
     * stress in the printed polymer, with the weaker strength across layers
     * local shell buckling from the surface curvature (minimal surface:
       k1 = -k2, so the radius follows from the Gaussian curvature)
     * every assembly joint between printed parts (tension across the cut)
     * deflection and overturning
   This is a screening model, not a design calculation: linear, isotropic
   stiffness, 2-3 voxels through the sheet (bending is under-estimated where
   the shell does not work as a membrane). A shell FE model by an engineer is
   the next step.

    python analysis.py                     # two-storey tower, stairs + structure
    python analysis.py -c my.json --no-fe  # stairs only
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra
import trimesh

import organic_pavilion as op
from pavilion import merge

STAIR_RULES = {
    "riser_min_m": 0.15,
    "riser_max_m": 0.17,
    "blondel_min_m": 0.62,
    "blondel_max_m": 0.64,
    "tread_min_m": 0.30,
    "width_min_m": 1.20,
    "max_risers_per_flight": 15,
    "landing_min_m": 1.20,
    "headroom_min_m": 2.10,
    "parapet_min_m": 1.00,
    "fill_max_m": 0.45,        # deepest printed fill under a tread before it needs its own support
}

# Large-format printed recycled PETG with glass fibre (typical published values).
MATERIAL = {
    "name": "rPETG-GF20, large-format printed",
    "E_MPa": 4500.0,
    "nu": 0.38,
    "density_kg_m3": 1350.0,
    "f_t_inplane_MPa": 45.0,     # along the print layers
    "f_t_interlayer_MPa": 18.0,  # across layers (weak direction, z when printed upright)
    "f_c_MPa": 55.0,
    "gamma_M": 1.5,
    "k_mod": 0.5,                # long-term load, outdoor temperature, creep (polymer)
    "joint_efficiency": 0.35,    # bolted/bonded joint vs. continuous sheet, tension across the cut
}

LOADS = {
    "live_kN_m2": 4.0,       # crowd on stairs and floors (EN 1991-1-1 / NTC 2018 cat. C2-C3)
    "snow_kN_m2": 1.2,       # characteristic on up-facing roof areas (low-altitude Italy, mu ~ 0.8)
    "wind_kN_m2": 0.8,       # peak pressure, simplified, horizontal on windward faces
}

TOWER = {
    "cell_m": [7.2, 7.2, 7.2],
    "envelope": {"radius_m": [3.95, 3.95, 3.67]},
    "pods": [[0, 0, 0], [0, 0, 1]],
}


# =========================================================================== #
# Stairs
# =========================================================================== #

def grid_index(axes, x, y, z):
    h = float(axes[0][1] - axes[0][0])
    ix = np.clip(np.round((x - axes[0][0]) / h).astype(int), 0, len(axes[0]) - 1)
    iy = np.clip(np.round((y - axes[1][0]) / h).astype(int), 0, len(axes[1]) - 1)
    iz = np.clip(np.round((z - axes[2][0]) / h).astype(int), 0, len(axes[2]) - 1)
    return ix, iy, iz


def headroom(F, axes, pts, start=0.12, stop=2.6):
    """Free height above each point (m), capped at `stop`."""
    h = float(axes[0][1] - axes[0][0])
    ix, iy, iz = grid_index(axes, pts[:, 0], pts[:, 1], pts[:, 2])
    free = np.full(len(pts), stop)
    open_ = np.ones(len(pts), bool)
    nz = F.shape[2]
    for k in range(int(start / h), int(stop / h) + 1):
        kk = np.clip(iz + k, 0, nz - 1)
        hit = open_ & (F[ix, iy, kk] <= 0) & (iz + k < nz)
        free[hit] = k * h
        open_ &= ~hit
    return free


def surface_top(F, axes, x, y, z_guess, window=0.6):
    """Top of the material directly under/around z_guess at (x, y), or None."""
    h = float(axes[0][1] - axes[0][0])
    ix, iy, _ = grid_index(axes, np.array([x]), np.array([y]), np.array([0.0]))
    col = F[ix[0], iy[0], :] <= 0
    zs = axes[2]
    sel = np.where(col & (zs > z_guess - window) & (zs < z_guess + window))[0]
    if len(sel) == 0:
        return None
    # tops of material runs = solid voxels with empty voxel above
    tops = [k for k in sel if k + 1 >= len(zs) or not col[k + 1]]
    if not tops:
        return None
    k = min(tops, key=lambda k: abs(zs[k] - z_guess))
    return float(zs[k] + h / 2)


FLIGHT = {
    "max_span_m": 4.0,         # longest printed flight between two landings on the shell
    "max_span_from_ground_m": 7.0,  # a first flight whose foot stands on the ground
    "pitch_min_deg": 8.0,
    "pitch_max_deg": 27.0,     # the code pitch atan(17/30); landings make a chain of flights flatter
    "stringer_m": 0.20,        # printed slab under the treads of a flight
    "cost_factor": 2.5,        # flights cost more than steps fused onto the shell
    "carve_cost": 0.6,         # per blocked sample: the sheet has to be opened for the stair
    "max_blocked": 70,
}


def _blocked(F, axes, a, b, flight=FLIGHT, rules=STAIR_RULES):
    """Number of sample points (out of 210) of the stair volume a->b that hit the
    sheet: headroom above the stair line over the full width, and the stringer
    below it. Zero means the flight fits as is; otherwise the sheet must be opened."""
    count = np.zeros(len(a), int)
    d = b - a
    hn = np.hypot(d[:, 0], d[:, 1]) + 1e-9
    px, py = -d[:, 1] / hn, d[:, 0] / hn                       # horizontal, across the flight
    half = rules["width_min_m"] / 2
    for t in np.linspace(0.1, 0.9, 14):
        p = a + d * t
        for off in (-half, 0.0, half):
            for dz in (0.12, 0.8, 1.5, rules["headroom_min_m"] - 0.05, -flight["stringer_m"]):
                ix, iy, iz = grid_index(axes, p[:, 0] + off * px, p[:, 1] + off * py, p[:, 2] + dz)
                count += F[ix, iy, iz] <= 0
    return count


def find_route(cfg, F, axes, mesh, target_pod, rules=STAIR_RULES, flight=FLIGHT, log=print):
    """Cheapest stair route from the ground to the floor of `target_pod`.

    Nodes: the walkable skin of the sheet (mesh vertices) plus a grid of points on
    the ground. Edges:
      * along the skin: steps fused onto the sheet, weighted by pitch, cross slope,
        distance to the free rim;
      * on the ground: plain walking;
      * flights: straight printed stair runs (pitch 8-31 deg, <= 4 m) between any
        two nodes, where headroom and stringer depth are free of the sheet.
    """
    from scipy.spatial import cKDTree

    v = mesh.vertices
    n = mesh.vertex_normals
    up = n[:, 2] > 0.5                                     # surface slope < 60 deg, facing up
    room = np.zeros(len(v))
    room[up] = headroom(F, axes, v[up] + n[up] * 0.02)
    xw, yw, zw = op.warp(cfg, v[:, 0], v[:, 1], v[:, 2])
    env = op.envelope_sdf(cfg, xw, yw, zw)
    walk = up & (room >= rules["headroom_min_m"]) & (env < -0.35)

    floor = op.chamber_clearance(F, axes, cfg)[",".join(map(str, target_pod))]["floor_m"]
    c = op.centre(cfg, target_pod)
    target = walk & (np.hypot(v[:, 0] - c[0], v[:, 1] - c[1]) < 1.6) & (np.abs(v[:, 2] - floor) < 0.35)

    # ---- skin edges
    e = mesh.edges_unique
    e = e[walk[e[:, 0]] & walk[e[:, 1]]]
    a, b = e[:, 0], e[:, 1]
    d = v[b] - v[a]
    horiz = np.hypot(d[:, 0], d[:, 1]) + 1e-9
    length = np.linalg.norm(d, axis=1)
    target_pitch = math.degrees(math.atan(0.16 / 0.31))
    nn = (n[a] + n[b]) / 2
    perp = np.stack([-d[:, 1] / horiz, d[:, 0] / horiz, np.zeros(len(d))], axis=1)
    cross = np.degrees(np.arcsin(np.clip(np.abs((nn * perp).sum(1)) / np.linalg.norm(nn, axis=1), 0, 1)))
    rim = (np.maximum(env[a], env[b]) > -0.7) * 2.0

    def cost(dz):
        p = np.degrees(np.arctan2(dz, horiz))
        c_ = 1 + 1.5 * ((p - target_pitch) / target_pitch) ** 2
        c_ = np.where(p > 31, c_ + 3 + 20 * ((p - 31) / 4) ** 2, c_)   # steeper than a stair: notched steps
        c_ = np.where(p > 38, np.inf, c_)
        c_ = np.where(p < -2, c_ + 12, c_)                               # going down on the way up
        c_ += 3 * (np.maximum(cross - 10, 0) / 10) ** 2
        c_ = np.where(cross > 35, np.inf, c_)
        return length * (c_ + rim)

    rows, cols, data, kinds = [], [], [], []
    carve_of = {}

    def add(r, c_, w, kind):
        ok = np.isfinite(w)
        rows.append(r[ok]); cols.append(c_[ok]); data.append(w[ok]); kinds.append(np.full(ok.sum(), kind))

    add(a, b, cost(d[:, 2]), 0)
    add(b, a, cost(-d[:, 2]), 0)

    # ---- ground nodes
    nv = len(v)
    gs = 0.3
    gx = np.arange(axes[0][0] + 0.2, axes[0][-1] - 0.2, gs)
    gy = np.arange(axes[1][0] + 0.2, axes[1][-1] - 0.2, gs)
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    gp = np.stack([GX.ravel(), GY.ravel(), np.full(GX.size, 0.0)], 1)
    gfree = headroom(F, axes, gp + [0, 0, 0.0], start=0.04, stop=rules["headroom_min_m"]) >= rules["headroom_min_m"]
    gp = gp[gfree]
    ng = len(gp)
    gid = nv + np.arange(ng)
    tree_g = cKDTree(gp[:, :2])
    pa, pb = np.array(list(tree_g.query_pairs(gs * 1.5))).T
    dl = np.linalg.norm(gp[pa] - gp[pb], axis=1)
    add(gid[pa], gid[pb], dl * 1.0, 1)
    add(gid[pb], gid[pa], dl * 1.0, 1)

    # ---- flights between sampled skin nodes and ground nodes
    wid = np.where(walk)[0]
    key = np.floor(v[wid] / 0.3).astype(int)
    _, first = np.unique(key, axis=0, return_index=True)
    sk = wid[first]
    fpts = np.concatenate([v[sk] + n[sk] * 0.01, gp])
    fid = np.concatenate([sk, gid])
    tree = cKDTree(fpts)
    pairs = np.array(list(tree.query_pairs(flight["max_span_m"])))
    # first flights may be longer: their foot stands on the ground
    t_sk = cKDTree(fpts[:len(sk)])
    gpairs = [(len(sk) + gi, si) for gi, lst in enumerate(t_sk.query_ball_point(gp, flight["max_span_from_ground_m"]))
              for si in lst]
    if gpairs:
        pairs = np.concatenate([pairs.reshape(-1, 2), np.array(gpairs)])
        pairs = np.unique(np.sort(pairs, axis=1), axis=0)
    log(f"route: {walk.sum()} walkable skin vertices, {ng} ground points, {len(pairs)} candidate flights")
    if len(pairs):
        p, q = pairs[:, 0], pairs[:, 1]
        lowfirst = fpts[p, 2] <= fpts[q, 2]
        lo_ = np.where(lowfirst, p, q); hi_ = np.where(lowfirst, q, p)
        dd = fpts[hi_] - fpts[lo_]
        hz = np.hypot(dd[:, 0], dd[:, 1])
        pitch = np.degrees(np.arctan2(dd[:, 2], hz + 1e-9))
        sel = (pitch >= flight["pitch_min_deg"]) & (pitch <= flight["pitch_max_deg"]) & (hz > 0.6)
        lo_, hi_, dd = lo_[sel], hi_[sel], dd[sel]
        blk = np.zeros(len(lo_), int)
        for c0 in range(0, len(lo_), 100000):
            blk[c0:c0 + 100000] = _blocked(F, axes, fpts[lo_[c0:c0 + 100000]], fpts[hi_[c0:c0 + 100000]])
        ok = blk <= flight["max_blocked"]
        lo_, hi_, dd, blk = lo_[ok], hi_[ok], dd[ok], blk[ok]
        L = np.linalg.norm(dd, axis=1)
        w = (L * flight["cost_factor"] + 1.0 + 1.2 * np.floor(dd[:, 2] / 2.5)   # a landing per 15 risers
             + flight["carve_cost"] * blk)                                    # opening the sheet
        add(fid[lo_], fid[hi_], w, 2)
        for a_, b_, k_ in zip(fid[lo_], fid[hi_], blk):
            carve_of[(int(a_), int(b_))] = int(k_)
        log(f"route: {(blk == 0).sum()} flights fit as they are, {(blk > 0).sum()} need an opening in the sheet")

    S = nv + ng
    R_ = np.concatenate(rows + [np.full(ng, S)])
    C_ = np.concatenate(cols + [gid])
    D_ = np.concatenate(data + [np.full(ng, 1e-6)])
    K_ = np.concatenate(kinds + [np.full(ng, 1)])
    order = np.lexsort((D_, C_, R_))
    R_, C_, D_, K_ = R_[order], C_[order], D_[order], K_[order]
    first_edge = np.ones(len(R_), bool)
    first_edge[1:] = (R_[1:] != R_[:-1]) | (C_[1:] != C_[:-1])   # keep the cheapest of duplicate edges
    R_, C_, D_, K_ = R_[first_edge], C_[first_edge], D_[first_edge], K_[first_edge]
    G = sparse.csr_matrix((D_, (R_, C_)), shape=(S + 1, S + 1))
    kind_of = {(int(r), int(c_)): int(k) for r, c_, k in zip(R_, C_, K_) if k == 2}
    dist, pred = dijkstra(G, directed=True, indices=S, return_predecessors=True)
    tgt = np.where(target)[0]
    allp = np.concatenate([v, gp])
    reach = np.isfinite(dist[:S])
    log(f"route: {reach.sum()} nodes reachable, highest at z = {allp[reach, 2].max():.2f} m; "
        f"{len(tgt)} target nodes on the floor at {floor:.2f} m")
    if len(tgt) == 0 or not np.isfinite(dist[tgt]).any():
        return None
    end = tgt[np.argmin(dist[tgt])]
    path = []
    k = end
    while k != S and k >= 0:
        path.append(int(k))
        k = pred[k]
    path = path[::-1]
    pts = allp[path]
    # segment type for each step of the path: 0 skin, 1 ground, 2 flight
    seg, carve = [], []
    for i in range(len(path) - 1):
        a_, b_ = path[i], path[i + 1]
        seg.append(2 if (a_, b_) in kind_of else (1 if a_ >= nv and b_ >= nv else 0))
        carve.append(carve_of.get((a_, b_), 0) if seg[-1] == 2 else 0)
    # drop the walk on the ground before the first rise
    k0 = 0
    while k0 < len(seg) and seg[k0] == 1:
        k0 += 1
    pts, seg, carve = pts[k0:], seg[k0:], carve[k0:]
    return {"points": pts, "segments": seg, "carve": carve, "floor_m": floor, "cost": float(dist[end])}


def resample_route(F, axes, pts, seg, step=0.05, sigma=0.3):
    """Resample by horizontal length. Skin stretches are smoothed in plan and dropped
    back onto the surface; flights stay straight lines."""
    out, kind = [], []
    i = 0
    while i < len(seg):
        j = i
        while j < len(seg) and seg[j] == seg[i]:
            j += 1
        chunk = pts[i:j + 1]
        hs = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(chunk[:, :2], axis=0).T))])
        S = np.arange(0, hs[-1], step) if hs[-1] > step else np.array([0.0])
        x = np.interp(S, hs, chunk[:, 0]); y = np.interp(S, hs, chunk[:, 1]); z = np.interp(S, hs, chunk[:, 2])
        if seg[i] == 0:
            k = int(sigma / step) * 3
            if len(S) > 2 * k + 2:
                w = np.exp(-0.5 * (np.arange(-k, k + 1) * step / sigma) ** 2); w /= w.sum()
                x = np.convolve(np.pad(x, k, mode="edge"), w, mode="valid")
                y = np.convolve(np.pad(y, k, mode="edge"), w, mode="valid")
            zz = []
            for xi, yi, zi in zip(x, y, z):
                zt = surface_top(F, axes, xi, yi, zi, window=0.3)
                zz.append(zi if zt is None else zt)
            z = np.array(zz)
        for xi, yi, zi in zip(x, y, z):
            out.append((xi, yi, zi)); kind.append(seg[i])
        i = j
    out.append(tuple(pts[-1])); kind.append(seg[-1] if seg else 0)
    out = np.array(out)
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(out[:, :2], axis=0).T))])
    return out, s, np.array(kind)


def fit_stair(route, s, kind, floor, rules=STAIR_RULES):
    """Regular stair laid along the route: constant riser, Blondel tread, a landing
    at least every 15 risers. On the sheet the treads may not cut into it (short
    tread) nor float more than fill_max above it (landing until it rises). In a
    printed flight the treads are free."""
    shell = kind == 0
    zr = route[:, 2]
    H = floor
    n = math.ceil(H / rules["riser_max_m"] - 1e-9)
    R = H / n
    G = min(max(0.63 - 2 * R, rules["tread_min_m"]), 0.34)
    fill_max = rules["fill_max_m"]
    end = s[-1]

    def at(sq):
        k = int(np.clip(np.searchsorted(s, sq), 0, len(s) - 1))
        return k

    steps = []
    s_r = 0.0
    flight, in_flight = 1, 0
    for i in range(1, n + 1):
        z_top = i * R
        if i > 1:
            prev_top = (i - 1) * R
            nxt = s_r + G
            landing = in_flight >= rules["max_risers_per_flight"]
            k = at(nxt)
            if nxt < end and shell[k] and z_top - zr[k] > fill_max + R:
                wait = np.where((s >= nxt) & (~shell | (zr >= z_top - fill_max - R)))[0]
                nxt = s[wait[0]] if len(wait) else end
                landing = True
            if landing:
                nxt = max(nxt, s_r + rules["landing_min_m"])
            hit = np.where((s > s_r + 1e-6) & shell & (zr > prev_top + 0.025))[0]
            if len(hit):
                nxt = min(nxt, s[hit[0]])
            tread = nxt - s_r
            steps[-1]["tread_m"] = round(tread, 3)
            steps[-1]["kind"] = "landing" if tread >= rules["landing_min_m"] - 1e-6 else "tread"
            if steps[-1]["kind"] == "landing":
                flight += 1
                in_flight = 0
        s_r = nxt if i > 1 else 0.0
        in_flight += 1
        steps.append({"n": i, "s_m": round(s_r, 3), "z_top_m": round(z_top, 3), "riser_m": round(R, 4),
                      "tread_m": None, "kind": "tread", "flight": flight})
    steps[-1]["tread_m"] = round(max(end - s_r, G), 3)
    steps[-1]["kind"] = "arrival"
    overshoot = max(0.0, steps[-1]["s_m"] - end)

    # place each tread on the route (extrapolating past its end if the stair overshoots)
    d_end = route[-1, :2] - route[max(len(route) - 6, 0), :2]
    d_end = d_end / (np.hypot(*d_end) + 1e-9)
    for st in steps:
        s0 = st["s_m"]; s1 = s0 + st["tread_m"]
        sm = (s0 + s1) / 2
        if sm <= end:
            k = at(sm)
            xy = route[k, :2]
            t = route[min(k + 3, len(route) - 1), :2] - route[max(k - 3, 0), :2]
        else:
            xy = route[-1, :2] + d_end * (sm - end)
            t = d_end
        t = t / (np.hypot(*t) + 1e-9)
        mask = (s >= s0) & (s <= s1)
        ks = kind[mask] if mask.any() else np.array([kind[at(min(sm, end))]])
        st["support"] = {0: "shell", 1: "ground", 2: "flight"}[int(np.bincount(ks, minlength=3).argmax())]
        sup = mask & shell
        st["fill_max_m"] = round(float(st["z_top_m"] - zr[sup].min()), 3) if sup.any() else 0.0
        st["clash_m"] = round(float(max(0.0, zr[sup].max() - st["z_top_m"])), 3) if sup.any() else 0.0
        st["xyz"] = [round(float(xy[0]), 3), round(float(xy[1]), 3), st["z_top_m"]]
        st["dir"] = [round(float(t[0]), 4), round(float(t[1]), 4)]
    return {"risers": n, "riser_m": round(R, 4), "tread_m": round(G, 3), "blondel_m": round(2 * R + G, 3),
            "pitch_deg": round(math.degrees(math.atan(R / G)), 1), "overshoot_m": round(overshoot, 2),
            "steps": steps}


def lateral_width(F, axes, st, rules=STAIR_RULES, max_half=1.0):
    """Usable width at a tread: how far headroom (and, on the shell, support) reach to each side."""
    x, y, _ = st["xyz"]; zt = st["z_top_m"]
    px, py = -st["dir"][1], st["dir"][0]
    sides = {}
    for sgn, name in [(1, "left"), (-1, "right")]:
        reach, why = max_half, "free"
        for dd in np.arange(0.05, max_half + 1e-6, 0.05):
            qx, qy = x + sgn * dd * px, y + sgn * dd * py
            room = headroom(F, axes, np.array([[qx, qy, zt]]), start=0.06, stop=rules["headroom_min_m"])[0]
            if room < rules["headroom_min_m"] - 1e-6:
                reach, why = dd - 0.05, "wall"
                break
            if st["support"] == "shell":
                top = surface_top(F, axes, qx, qy, zt - 0.2, window=rules["fill_max_m"] + 0.3)
                if top is None:
                    reach, why = dd - 0.05, "open edge"
                    break
        if why == "free" and st["support"] != "ground" and zt > 0.5:
            why = "open edge"   # the flight continues past this width with a drop beside it
        sides[name] = (round(reach, 2), why)
    return sides


def check_stair(fit, widths, rules=STAIR_RULES):
    issues = []
    R, G = fit["riser_m"], fit["tread_m"]
    if not (rules["riser_min_m"] <= R <= rules["riser_max_m"]):
        issues.append(f"riser {R * 100:.1f} cm outside {rules['riser_min_m'] * 100:.0f}-{rules['riser_max_m'] * 100:.0f} cm")
    if not (rules["blondel_min_m"] - 1e-3 <= 2 * R + G <= rules["blondel_max_m"] + 1e-3):
        issues.append(f"2R+G = {(2 * R + G) * 100:.1f} cm outside 62-64 cm")
    rows = []
    for st, sides in zip(fit["steps"], widths):
        left, right = sides["left"], sides["right"]
        w = left[0] + right[0]
        flags = []
        if st["kind"] == "tread" and st["tread_m"] < rules["tread_min_m"] - 1e-3:
            flags.append("short tread: support steeper than the stair")
        if st["kind"] == "tread" and rules["tread_min_m"] + 0.06 < st["tread_m"] < rules["landing_min_m"] - 1e-3:
            flags.append("irregular tread")
        if st["clash_m"] > 0.03:
            flags.append(f"sheet {st['clash_m'] * 100:.0f} cm above tread")
        if st["support"] == "shell" and st["fill_max_m"] > rules["fill_max_m"]:
            flags.append(f"fill {st['fill_max_m']:.2f} m")
        if w < rules["width_min_m"] - 1e-3:
            flags.append(f"width {w:.2f} m")
        if "open edge" in (left[1], right[1]):
            flags.append("parapet")
        rows.append({**st, "width_m": round(w, 2), "left": left, "right": right, "flags": flags})
    flights, count = [], 0
    for r in rows:
        count += 1
        if r["kind"] in ("landing", "arrival"):
            flights.append(count); count = 0
    if count:
        flights.append(count)
    if max(flights) > rules["max_risers_per_flight"]:
        issues.append(f"a flight has {max(flights)} risers (max {rules['max_risers_per_flight']})")
    narrow = [r["n"] for r in rows if r["width_m"] < rules["width_min_m"] - 1e-3]
    if narrow:
        issues.append(f"{len(narrow)} treads narrower than {rules['width_min_m']:.2f} m "
                      f"(min {min(r['width_m'] for r in rows):.2f} m)")
    for flag, label in [("short tread", "treads where the support is steeper than the stair"),
                        ("irregular tread", "irregular treads"), ("sheet", "treads clashing with the sheet")]:
        k = [r["n"] for r in rows if any(f.startswith(flag) for f in r["flags"])]
        if k:
            issues.append(f"{len(k)} {label} (steps {k[0]}-{k[-1]})")
    return rows, flights, issues


def carve_openings(F, axes, rows, rules=STAIR_RULES, width=1.2, edge=0.15):
    """Open the sheet where the stair needs passage: the volume above each tread
    (full width, 2.10 m headroom). Returns the new field and the openings found,
    grouped by consecutive treads, with their size and position."""
    h = float(axes[0][1] - axes[0][0])
    F0 = F
    F = F.copy()
    removed = []
    for r in rows:
        x, y, zt = r["xyz"]
        dx, dy = r["dir"]
        L = (r["tread_m"] or 0.3) + 0.02
        hw = width / 2 + edge
        rad = math.hypot(L / 2, hw) + 0.1
        i0, j0, _ = grid_index(axes, np.array([x - rad]), np.array([y - rad]), np.array([0.0]))
        i1, j1, _ = grid_index(axes, np.array([x + rad]), np.array([y + rad]), np.array([0.0]))
        kz = np.where((axes[2] > zt + 0.005) & (axes[2] < zt + rules["headroom_min_m"] + 0.08))[0]
        zs = axes[2][kz]
        count = 0
        for i in range(i0[0], i1[0] + 1):
            for j in range(j0[0], j1[0] + 1):
                px, py = axes[0][i] - x, axes[1][j] - y
                u = px * dx + py * dy
                vv = -px * dy + py * dx
                d_plan = max(abs(u) - L / 2, abs(vv) - hw)
                if d_plan > h:
                    continue
                box = np.maximum.reduce([np.full(len(kz), d_plan), zt + 0.005 - zs,
                                         zs - (zt + rules["headroom_min_m"] + 0.05)])
                before = F[i, j, kz] <= 0
                F[i, j, kz] = np.maximum(F[i, j, kz], -box)
                count += int((before & (F[i, j, kz] > 0)).sum())
        removed.append(count * h ** 3)
    openings = []
    cur = None
    for r, vol in zip(rows, removed):
        if vol > 2e-3:
            if cur is None:
                cur = {"steps": [r["n"], r["n"]], "volume_m3": 0.0, "pts": []}
            cur["steps"][1] = r["n"]; cur["volume_m3"] += vol; cur["pts"].append(r["xyz"])
        elif cur is not None:
            openings.append(cur); cur = None
    if cur is not None:
        openings.append(cur)
    for o in openings:
        p = np.array(o.pop("pts"))
        o["centre_m"] = [round(float(c), 2) for c in p.mean(0)]
        o["volume_m3"] = round(o["volume_m3"], 3)
        o["area_m2"] = round(o["volume_m3"] / 0.08, 1)   # sheet removed / mean sheet thickness
    return F, openings


def stair_solids(F, axes, rows, rules=STAIR_RULES, flight=FLIGHT, width=1.2, tie=0.06):
    """Fuse the treads into the field. On the shell each tread is a block from the
    sheet up to the tread; in a flight it sits on a printed stringer slab; from the
    ground it is a solid plinth."""
    h = float(axes[0][1] - axes[0][0])
    F = F.copy()
    for r in rows:
        if r["kind"] == "arrival":
            continue
        x, y, _ = r["xyz"]
        if r["support"] == "flight" and r["kind"] == "landing":
            pass
        dx, dy = r["dir"]
        L = r["tread_m"] + 0.02
        left, right = r["left"][0], r["right"][0]
        w_l = min(left, width / 2 + max(0, width / 2 - right))
        w_r = min(right, width / 2 + max(0, width / 2 - left))
        off = (w_l - w_r) / 2
        hw = (w_l + w_r) / 2
        if hw < 0.2:
            continue
        zt = r["z_top_m"]
        if r["support"] == "flight":
            depth = r["riser_m"] + flight["stringer_m"]
        elif r["support"] == "ground" or zt <= rules["fill_max_m"] + r["riser_m"]:
            depth = zt + 0.05
        else:
            depth = r["fill_max_m"] + 0.25
        rad = math.hypot(L / 2, hw) + abs(off) + 0.1
        i0, j0, _ = grid_index(axes, np.array([x - rad]), np.array([y - rad]), np.array([0.0]))
        i1, j1, _ = grid_index(axes, np.array([x + rad]), np.array([y + rad]), np.array([0.0]))
        kz = np.where((axes[2] > zt - depth - 0.02) & (axes[2] < zt + 0.03))[0]
        if len(kz) == 0:
            continue
        zs = axes[2][kz]
        for i in range(i0[0], i1[0] + 1):
            for j in range(j0[0], j1[0] + 1):
                px, py = axes[0][i] - x, axes[1][j] - y
                u = px * dx + py * dy
                vv = -px * dy + py * dx - off
                d_plan = max(abs(u) - L / 2, abs(vv) - hw)
                if d_plan > h:
                    continue
                if r["support"] == "shell":
                    col = F[i, j, kz] <= 0
                    below = np.where(col & (zs < zt - 0.02))[0]
                    base = zs[below[-1]] - tie if len(below) else zt - depth
                else:
                    base = max(zt - depth, -0.01)
                sdf = np.maximum.reduce([np.full(len(kz), d_plan), zs - zt, base - zs])
                F[i, j, kz] = np.minimum(F[i, j, kz], sdf)
    return F



# =========================================================================== #
# Climbing the terraced sculpture
# =========================================================================== #

CLIMB_CLASSES = [
    # name, max surface slope (deg), what it is for
    ("seat", 12.0, "seats and landings: ledges deeper than 0.85 m"),
    ("steps", 32.0, "comfortable steps: treads 0.29-0.85 m"),
    ("steep steps", 45.0, "steep steps: treads 0.18-0.29 m, a hand on the next ledge"),
    ("ledges", 58.0, "climbing ledges: treads 0.11-0.18 m, hands and feet"),
]


def slope_class(slope_deg):
    c = np.full(slope_deg.shape, len(CLIMB_CLASSES), int)   # last = too steep / smooth
    for k in range(len(CLIMB_CLASSES) - 1, -1, -1):
        c[slope_deg < CLIMB_CLASSES[k][1]] = k
    return c


def ledge_classes(mesh, mesh0, radius=0.3):
    """Ledge class for every vertex of the terraced mesh, from the slope of the
    untreaded up-facing skin under it (averaged over a small neighbourhood so the
    classes form clean bands). Undersides and anything far from that skin are smooth."""
    from scipy.spatial import cKDTree
    n0 = mesh0.vertex_normals
    up = np.where(n0[:, 2] > 0.05)[0]
    slope0 = np.degrees(np.arccos(np.clip(n0[up, 2], -1, 1)))
    tree = cKDTree(mesh0.vertices[up])
    d, idx = tree.query(mesh.vertices, k=24, distance_upper_bound=radius)
    valid = np.isfinite(d)
    sl = np.where(valid, slope0[np.minimum(idx, len(up) - 1)], 0.0)
    cnt = valid.sum(1)
    mean = np.where(cnt > 0, sl.sum(1) / np.maximum(cnt, 1), 90.0)
    cls = slope_class(mean)
    # underside: nearer to the plain sheet's down-facing skin than to its up-facing one
    # (a vertex's own normal is no guide on the ledges: risers and noses face sideways)
    down = np.where(n0[:, 2] < -0.05)[0]
    d_dn, _ = cKDTree(mesh0.vertices[down]).query(mesh.vertices)
    d_up, _ = cKDTree(mesh0.vertices[up]).query(mesh.vertices)
    cls[(cnt == 0) | (d_dn < d_up)] = len(CLIMB_CLASSES)
    # majority vote over the 1-ring, twice, for clean bands
    nb = mesh.vertex_neighbors
    k = len(CLIMB_CLASSES) + 1
    for _ in range(2):
        cls = np.array([np.bincount(np.r_[cls[i], cls[nb[i]]], minlength=k).argmax() for i in range(len(cls))])
    return cls.astype(np.uint8)


def climb(cfg, F, axes, mesh, riser, headroom_min=1.9, rim_keep=0.3, log=print, avoid=None):
    """Best climbing line over the terraced skin, from the ground to the summit.

    Works on the untreaded sheet: wherever its up-facing skin is below 58 deg the
    terraces turn it into ledges, and the tread depth follows from the slope
    (tread = riser / tan(slope))."""
    v = mesh.vertices
    n = mesh.vertex_normals
    slope = np.degrees(np.arccos(np.clip(n[:, 2], -1, 1)))
    cls = slope_class(slope)
    up = cls < len(CLIMB_CLASSES)
    room = np.zeros(len(v))
    room[up] = headroom(F, axes, v[up] + n[up] * 0.02)
    xw, yw, zw = op.warp(cfg, v[:, 0], v[:, 1], v[:, 2])
    env = op.envelope_sdf(cfg, xw, yw, zw)
    ok = up & (room >= headroom_min) & (env < -rim_keep)
    if avoid is not None:                       # e.g. a room closed by membranes
        ok &= ~avoid(v + n * 0.3)

    # areas per class (whole sculpture, up-facing skin)
    fcls = slope_class(np.degrees(np.arccos(np.clip(mesh.face_normals[:, 2], -1, 1))))
    areas = {CLIMB_CLASSES[k][0]: round(float(mesh.area_faces[fcls == k].sum()), 1) for k in range(len(CLIMB_CLASSES))}

    e = mesh.edges_unique
    e = e[ok[e[:, 0]] & ok[e[:, 1]]]
    L = np.linalg.norm(v[e[:, 1]] - v[e[:, 0]], axis=1)
    weight = np.array([1.0, 1.0, 1.8, 3.5])[np.maximum(cls[e[:, 0]], cls[e[:, 1]])]
    G = sparse.coo_matrix((np.concatenate([L * weight] * 2),
                           (np.concatenate([e[:, 0], e[:, 1]]), np.concatenate([e[:, 1], e[:, 0]]))),
                          shape=(len(v) + 1, len(v) + 1)).tocsr()
    src = np.where(ok & (v[:, 2] < riser + 0.1))[0]
    S = len(v)
    G = G.tolil(); G[S, src] = 1e-6; G = G.tocsr()
    dist, pred = dijkstra(G, directed=True, indices=S, return_predecessors=True)
    reach = np.isfinite(dist[:S]) & ok
    if not reach.any():
        return None
    summit = np.where(reach)[0][np.argmax(v[reach, 2])]
    path = []
    k = summit
    while k != S and k >= 0:
        path.append(int(k)); k = pred[k]
    path = path[::-1]
    pts = v[path]
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    pc = cls[path]
    comp = {CLIMB_CLASSES[k][0]: round(float(seg[np.minimum(pc[:-1], pc[1:]) == k].sum()), 2)
            for k in range(len(CLIMB_CLASSES))}
    # fall height beside the line: drop to the next surface (or the ground) 0.6 m to either side
    drops = []
    for i in range(0, len(pts), max(1, len(pts) // 120)):
        p = pts[i]
        t = pts[min(i + 3, len(pts) - 1)] - pts[max(i - 3, 0)]
        hn = math.hypot(t[0], t[1]) + 1e-9
        px, py = -t[1] / hn, t[0] / hn
        worst = 0.0
        for sg in (1, -1):
            qx, qy = p[0] + sg * 0.6 * px, p[1] + sg * 0.6 * py
            ix, iy, iz = grid_index(axes, np.array([qx]), np.array([qy]), np.array([p[2]]))
            col = F[ix[0], iy[0], :iz[0] + 1] <= 0
            below = np.where(col)[0]
            z_below = axes[2][below[-1]] if len(below) else 0.0
            worst = max(worst, p[2] - max(z_below, 0.0))
        drops.append((float(s[i]), float(p[2]), worst))
    drops = np.array(drops)
    guard = drops[:, 2] > 1.0
    log(f"climb: summit {v[summit, 2]:.2f} m reached over {s[-1]:.1f} m, {int(round(v[summit, 2] / riser))} ledges; "
        f"{guard.mean() * 100:.0f}% of the line has a drop over 1 m beside it")
    return {
        "riser_m": riser,
        "summit_m": round(float(v[summit, 2]), 2),
        "top_of_sculpture_m": round(float(v[:, 2].max()), 2),
        "length_m": round(float(s[-1]), 2),
        "ledges_climbed": int(round(v[summit, 2] / riser)),
        "composition_m": comp,
        "areas_m2": areas,
        "reachable_share": round(float(reach.sum() / max(ok.sum(), 1)), 3),
        "max_drop_m": round(float(drops[:, 2].max()), 2),
        "guarded_share": round(float(guard.mean()), 3),
        "classes": [{"name": c[0], "max_slope_deg": c[1], "use": c[2]} for c in CLIMB_CLASSES],
        "route": np.round(pts, 3).tolist(),
        "route_class": pc.tolist(),
        "drops": np.round(drops, 2).tolist(),
    }

# =========================================================================== #
# Structure: voxel finite elements
# =========================================================================== #

def hex_stiffness(E, nu, h):
    """24x24 stiffness of a cubic 8-node hexahedron of side h (2x2x2 Gauss)."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam; D[np.arange(3), np.arange(3)] += 2 * mu; D[3:, 3:] = np.eye(3) * mu
    corners = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]) * 2 - 1
    g = 1 / math.sqrt(3)
    K = np.zeros((24, 24))
    for gp in np.array(np.meshgrid([-g, g], [-g, g], [-g, g])).T.reshape(-1, 3):
        B = b_matrix(corners, gp, h)
        K += B.T @ D @ B * (h / 2) ** 3
    return K, D


def b_matrix(corners, xi, h):
    dN = np.zeros((8, 3))
    for a, c in enumerate(corners):
        dN[a, 0] = c[0] * (1 + c[1] * xi[1]) * (1 + c[2] * xi[2]) / 8
        dN[a, 1] = c[1] * (1 + c[0] * xi[0]) * (1 + c[2] * xi[2]) / 8
        dN[a, 2] = c[2] * (1 + c[0] * xi[0]) * (1 + c[1] * xi[1]) / 8
    dN *= 2 / h
    B = np.zeros((6, 24))
    for a in range(8):
        B[0, 3 * a] = dN[a, 0]; B[1, 3 * a + 1] = dN[a, 1]; B[2, 3 * a + 2] = dN[a, 2]
        B[3, 3 * a] = dN[a, 1]; B[3, 3 * a + 1] = dN[a, 0]
        B[4, 3 * a + 1] = dN[a, 2]; B[4, 3 * a + 2] = dN[a, 1]
        B[5, 3 * a] = dN[a, 2]; B[5, 3 * a + 2] = dN[a, 0]
    return B


def gaussian_curvature_radius(f, h):
    """Radius 1/sqrt|K| of the level surfaces of f (implicit-surface formula)."""
    fx, fy, fz = np.gradient(f, h)
    fxx, fxy, fxz = np.gradient(fx, h)
    _, fyy, fyz = np.gradient(fy, h)
    _, _, fzz = np.gradient(fz, h)
    g2 = fx * fx + fy * fy + fz * fz + 1e-12
    num = (fx * fx * (fyy * fzz - fyz * fyz) + fy * fy * (fxx * fzz - fxz * fxz) + fz * fz * (fxx * fyy - fxy * fxy)
           + 2 * (fx * fy * (fxz * fyz - fxy * fzz) + fy * fz * (fxy * fxz - fyz * fxx) + fx * fz * (fxy * fyz - fxz * fyy)))
    K = num / (g2 * g2)
    return 1 / np.sqrt(np.abs(K) + 1e-6)


def structure(cfg, F, axes, extras, rows, material=MATERIAL, loads=LOADS, log=print, crowd_everywhere=False,
              extra_forces=None):
    from pyamg import smoothed_aggregation_solver

    h = float(axes[0][1] - axes[0][0])
    nx, ny, nz = F.shape
    solid = F <= 0
    # drop voxels below ground level (the marching-cubes padding)
    solid &= (axes[2] >= -1e-6)[None, None, :]
    eid = -np.ones(F.shape, dtype=np.int64)
    eidx = np.argwhere(solid)
    ne = len(eidx)
    eid[tuple(eidx.T)] = np.arange(ne)
    log(f"FE: {ne} elements of {h * 1000:.0f} mm")

    # nodes: corners of solid voxels
    offs = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    corner = (eidx[:, None, :] + offs[None, :, :]).reshape(-1, 3)
    key = (corner[:, 0] * (ny + 1) + corner[:, 1]) * (nz + 1) + corner[:, 2]
    ukeys, inv = np.unique(key, return_inverse=True)
    nn = len(ukeys)
    enodes = inv.reshape(ne, 8)
    edof = (enodes[:, :, None] * 3 + np.arange(3)[None, None, :]).reshape(ne, 24)
    ndof = 3 * nn
    log(f"FE: {nn} nodes, {ndof} unknowns")

    E = material["E_MPa"] * 1e6
    KE, D = hex_stiffness(E, material["nu"], h)

    # assemble in chunks
    K = None
    ii = np.repeat(np.arange(24), 24); jj = np.tile(np.arange(24), 24)
    kv = KE.ravel()
    for c0 in range(0, ne, 40000):
        ed = edof[c0:c0 + 40000]
        r = ed[:, ii].ravel(); c = ed[:, jj].ravel()
        part = sparse.coo_matrix((np.tile(kv, len(ed)), (r, c)), shape=(ndof, ndof)).tocsr()
        K = part if K is None else K + part
    # supports: nodes at ground level fixed
    node_ijk = np.stack([ukeys // ((ny + 1) * (nz + 1)), (ukeys // (nz + 1)) % (ny + 1), ukeys % (nz + 1)], 1)
    z_node = axes[2][0] - h / 2 + node_ijk[:, 2] * h
    fixed_nodes = np.where(z_node <= h * 0.6)[0]
    fixed = np.zeros(ndof, bool)
    for d_ in range(3):
        fixed[fixed_nodes * 3 + d_] = True
    free = np.where(~fixed)[0]
    Kf = K[free][:, free].tocsr()
    log(f"FE: {len(fixed_nodes)} supported nodes, building multigrid")

    # rigid body modes for the multigrid near-nullspace
    xyz = np.stack([axes[0][0] - h / 2 + node_ijk[:, 0] * h, axes[1][0] - h / 2 + node_ijk[:, 1] * h, z_node], 1)
    xyz = xyz - xyz.mean(0)
    Bm = np.zeros((ndof, 6))
    Bm[0::3, 0] = 1; Bm[1::3, 1] = 1; Bm[2::3, 2] = 1
    Bm[0::3, 3] = -xyz[:, 1]; Bm[1::3, 3] = xyz[:, 0]
    Bm[1::3, 4] = -xyz[:, 2]; Bm[2::3, 4] = xyz[:, 1]
    Bm[0::3, 5] = xyz[:, 2]; Bm[2::3, 5] = -xyz[:, 0]
    ml = smoothed_aggregation_solver(Kf, B=Bm[free], strength=("symmetric", {"theta": 0.0}),
                                     smooth="jacobi", max_coarse=2000)

    # ---- loads
    vol = h ** 3
    g = 9.81
    w_self = material["density_kg_m3"] * g * vol           # N per element
    # up-facing skin elements (solid with empty voxel above) carry area loads
    above_empty = np.zeros_like(solid)
    above_empty[:, :, :-1] = solid[:, :, :-1] & ~solid[:, :, 1:]
    skin = above_empty[tuple(eidx.T)]
    # walkable areas: stair treads and chamber floors (skin within 0.4 m of a tread or a floor)
    zc = axes[2][eidx[:, 2]]
    xe, ye = axes[0][eidx[:, 0]], axes[1][eidx[:, 1]]
    walk = np.zeros(ne, bool)
    for r in rows:
        x, y, _ = r["xyz"]
        walk |= skin & (np.hypot(xe - x, ye - y) < 0.7) & (np.abs(zc - r["z_top_m"]) < 0.08)
    for p in op.pods(cfg):
        c = op.centre(cfg, p)
        fl = op.chamber_clearance(F, axes, cfg)[",".join(map(str, p))]["floor_m"]
        if fl > 0.05:
            walk |= skin & (np.hypot(xe - c[0], ye - c[1]) < 2.2) & (np.abs(zc - fl) < 0.6)
    if crowd_everywhere:              # a climbable sculpture: people anywhere on its up-facing skin
        walk = skin.copy()
    roof = skin & ~walk
    # wind: windward faces, +x direction, solid with empty voxel on the -x side
    wind_face = np.zeros_like(solid)
    wind_face[1:, :, :] = solid[1:, :, :] & ~solid[:-1, :, :]
    windward = wind_face[tuple(eidx.T)]

    def vector(per_elem_xyz):
        f = np.zeros(ndof)
        np.add.at(f, edof[:, 0::3].ravel(), np.repeat(per_elem_xyz[:, 0] / 8, 8))
        np.add.at(f, edof[:, 1::3].ravel(), np.repeat(per_elem_xyz[:, 1] / 8, 8))
        np.add.at(f, edof[:, 2::3].ravel(), np.repeat(per_elem_xyz[:, 2] / 8, 8))
        return f

    A = h * h
    G_ = np.zeros((ne, 3)); G_[:, 2] = -w_self
    Q_ = np.zeros((ne, 3)); Q_[walk, 2] = -loads["live_kN_m2"] * 1e3 * A
    S_ = np.zeros((ne, 3)); S_[roof, 2] = -loads["snow_kN_m2"] * 1e3 * A
    W_ = np.zeros((ne, 3)); W_[windward, 0] = loads["wind_kN_m2"] * 1e3 * A
    # membranes: prestress and wind arrive as point forces along their clamped edges
    Mp = np.zeros((ne, 3))
    if extra_forces is not None and len(extra_forces[0]):
        cent = np.stack([axes[a][eidx[:, a]] for a in range(3)], 1)
        from scipy.spatial import cKDTree
        _, near = cKDTree(cent).query(extra_forces[0])
        np.add.at(Mp, near, extra_forces[1])
        log(f"FE: {len(near)} membrane edge forces, resultant {np.round(Mp.sum(0) / 1e3, 1)} kN")
    cases = {
        "ULS crowd": 1.35 * G_ + 1.5 * Q_ + 0.75 * S_ + Mp,
        "ULS snow": 1.35 * G_ + 1.5 * S_ + 1.05 * Q_ + Mp,
        "ULS wind": 1.0 * G_ + 1.5 * W_ + Mp,
        "SLS": G_ + Q_ + 0.5 * S_ + Mp / 1.5,
    }

    # strain-displacement at the element centre for stresses
    corners = offs * 2 - 1
    Bc = b_matrix(corners, np.zeros(3), h)
    f = extras["f"]
    Rcurv = gaussian_curvature_radius(f, h)[tuple(eidx.T)]
    t_loc = extras["thickness"][tuple(eidx.T)]
    fm = material
    f_in = fm["f_t_inplane_MPa"] * fm["k_mod"] / fm["gamma_M"]
    f_z = fm["f_t_interlayer_MPa"] * fm["k_mod"] / fm["gamma_M"]
    f_c = fm["f_c_MPa"] * fm["k_mod"] / fm["gamma_M"]
    # local buckling: curved shell (classical value, 0.2 knock-down for imperfections) or, where
    # the minimal surface is nearly flat (its flat points, K -> 0), a plate panel of ~1/5 cell
    # (simply supported, 0.5 knock-down); whichever is higher governs
    b_panel = min(cfg["cell_m"][:2]) / 5
    sig_shell = 0.605 * fm["E_MPa"] * t_loc / np.minimum(Rcurv, 50.0) * 0.2
    sig_plate = 0.5 * 4 * math.pi ** 2 * fm["E_MPa"] / (12 * (1 - fm["nu"] ** 2)) * (t_loc / b_panel) ** 2
    sig_cr = np.maximum(sig_shell, sig_plate)

    results = {}
    util_env = np.zeros(ne)
    reason = np.zeros(ne, dtype=np.int8)  # 1 in-plane, 2 interlayer, 3 compression, 4 buckling
    total_load = {}
    for name, P in cases.items():
        fvec = vector(P)
        total_load[name] = [round(float(P[:, k].sum()) / 1e3, 1) for k in range(3)]
        u = np.zeros(ndof)
        res = []
        u[free] = ml.solve(fvec[free], tol=1e-7, maxiter=300, accel="cg", residuals=res)
        log(f"FE: {name}: {len(res)} iterations, residual {res[-1] / res[0]:.1e}")
        ue = u[edof]
        eps = ue @ Bc.T
        sig = eps @ D.T / 1e6                  # MPa
        if name == "SLS":
            disp = np.linalg.norm(u.reshape(-1, 3), axis=1)
            kmax = int(np.argmax(disp))
            results[name] = {"max_displacement_mm": round(float(disp.max() * 1e3), 1),
                             "at_m": [round(float(axes[0][0] - h / 2 + node_ijk[kmax, 0] * h), 2),
                                      round(float(axes[1][0] - h / 2 + node_ijk[kmax, 1] * h), 2),
                                      round(float(z_node[kmax]), 2)],
                             "node_disp": disp}
            continue
        S3 = np.zeros((ne, 3, 3))
        S3[:, 0, 0], S3[:, 1, 1], S3[:, 2, 2] = sig[:, 0], sig[:, 1], sig[:, 2]
        S3[:, 0, 1] = S3[:, 1, 0] = sig[:, 3]
        S3[:, 1, 2] = S3[:, 2, 1] = sig[:, 4]
        S3[:, 0, 2] = S3[:, 2, 0] = sig[:, 5]
        pr = np.linalg.eigvalsh(S3)
        s1, s3 = pr[:, 2], pr[:, 0]
        u_in = np.maximum(s1, 0) / f_in
        u_z = np.maximum(sig[:, 2], 0) / f_z           # tension across horizontal print layers
        u_c = np.maximum(-s3, 0) / f_c
        u_b = np.maximum(-s3, 0) / np.maximum(sig_cr, 1e-3)
        U = np.stack([u_in, u_z, u_c, u_b], 1)
        ue_max = U.max(1)
        better = ue_max > util_env
        util_env = np.where(better, ue_max, util_env)
        reason = np.where(better, U.argmax(1) + 1, reason)
        results[name] = {
            "max_principal_tension_MPa": round(float(np.percentile(s1, 99.9)), 2),
            "max_compression_MPa": round(float(-np.percentile(s3, 0.1)), 2),
            "max_interlayer_tension_MPa": round(float(np.percentile(sig[:, 2], 99.9)), 2),
            "utilisation_p99_9": round(float(np.percentile(ue_max, 99.9)), 2),
            "sigma": sig,
        }

    # ---- joints between print parts (cut planes of the build-volume grid)
    lo = np.array([axes[0][0], axes[1][0], 0.0])
    bv = np.array(cfg["printer"]["build_volume_m"])
    f_joint = f_in * fm["joint_efficiency"]
    joints = []
    for ax in range(3):
        cuts = lo[ax] + bv[ax] * np.arange(1, 40)
        cuts = cuts[cuts < axes[ax][-1]]
        for cpos in cuts:
            k = int(round((cpos - axes[ax][0]) / h))  # voxel layer just above the cut
            if k <= 0 or k >= F.shape[ax]:
                continue
            sl_a = [slice(None)] * 3; sl_b = [slice(None)] * 3
            sl_a[ax] = k - 1; sl_b[ax] = k
            ea = eid[tuple(sl_a)]; eb = eid[tuple(sl_b)]
            both = (ea >= 0) & (eb >= 0)
            if not both.any():
                continue
            ids = np.concatenate([ea[both], eb[both]])
            worst = 0.0; worst_case = None; force = 0.0
            for name in ("ULS crowd", "ULS snow", "ULS wind"):
                sn = results[name]["sigma"][ids, ax]   # normal stress across the cut
                tens = np.percentile(sn, 99) if len(sn) > 20 else sn.max()
                if tens > worst:
                    worst, worst_case = float(tens), name
                    force = float(np.clip(sn, 0, None).sum() * 1e6 * h * h / 2 / 1e3)
            # position of the most stressed contact (by plan quadrant)
            idx = np.argwhere(both)
            centre = [float(axes[a][int(np.median(idx[:, j]))]) for j, a in enumerate([a for a in range(3) if a != ax])]
            joints.append({"axis": "xyz"[ax], "at_m": round(float(cpos), 2),
                           "contact_area_m2": round(float(both.sum() * h * h), 3),
                           "max_tension_MPa": round(worst, 2), "case": worst_case,
                           "tension_force_kN": round(force, 1),
                           "utilisation": round(worst / f_joint, 2), "centre_m": [round(c, 2) for c in centre]})
    joints.sort(key=lambda j: -j["utilisation"])

    # ---- critical zones: cluster elements with utilisation > 0.8
    names = {1: "tension along layers", 2: "tension across print layers", 3: "compression", 4: "local buckling"}
    # the most loaded 0.3% of the material, or anything above 0.8, whichever is more
    hot = util_env > min(0.8, float(np.percentile(util_env, 99.7)))
    zones = []
    if hot.any():
        from scipy import ndimage
        grid = np.zeros(F.shape, bool)
        grid[tuple(eidx[hot].T)] = True
        lab, nlab = ndimage.label(ndimage.binary_dilation(grid, iterations=3))
        for L_ in range(1, nlab + 1):
            m = hot & (lab[tuple(eidx.T)] == L_)
            if m.sum() < 3:
                continue
            p = eidx[m]
            j = np.argmax(np.where(m, util_env, 0))
            if m.sum() < 8:
                continue
            zones.append({"elements": int(m.sum()),
                          "centre_m": [round(float(axes[a][int(p[:, a].mean())]), 2) for a in range(3)],
                          "max_utilisation": round(float(util_env[j]), 2),
                          "governing": names[int(reason[j])]})
        zones.sort(key=lambda z: -z["max_utilisation"] * math.log(1 + z["elements"]))

    # ---- overturning (wind on +x, about the leeward footprint edge)
    xs_e = axes[0][eidx[:, 0]]
    weight = w_self * ne
    x_lee = xs_e[eidx[:, 2] == eidx[:, 2].min()].max()
    x_g = (xs_e * w_self).sum() / weight
    wind_force = W_[:, 0].sum()
    z_w = (zc * W_[:, 0]).sum() / max(wind_force, 1e-9)
    m_over = 1.5 * wind_force * z_w
    m_rest = 0.9 * weight * (x_lee - x_g)

    sls = results["SLS"]
    span = max(cfg["cell_m"][0], cfg["cell_m"][1])   # limit span/250 on the cell size
    out = {
        "elements": ne, "unknowns": int(len(free)), "voxel_m": h,
        "material": material, "loads": loads,
        "total_load_kN": total_load,
        "cases": {k: {kk: vv for kk, vv in v.items() if kk not in ("sigma", "node_disp")} for k, v in results.items()},
        "max_displacement_mm": sls["max_displacement_mm"],
        "max_displacement_at_m": sls["at_m"],
        "displacement_limit_mm": round(span / 250 * 1e3, 1),
        "design_strength_MPa": {"in_plane_tension": round(f_in, 2), "interlayer_tension": round(f_z, 2),
                                "compression": round(f_c, 2), "joint_tension": round(f_joint, 2)},
        "utilisation": {"p99": round(float(np.percentile(util_env, 99)), 2),
                        "p99_9": round(float(np.percentile(util_env, 99.9)), 2),
                        "max": round(float(util_env.max()), 2),
                        "share_over_1": round(float((util_env > 1).mean()), 4)},
        "buckling_radius_m": {"min": round(float(Rcurv.min()), 2), "median": round(float(np.median(Rcurv)), 2)},
        "buckling_stress_MPa": {"p1": round(float(np.percentile(sig_cr, 1)), 1),
                                "median": round(float(np.median(sig_cr)), 1)},
        "overturning": {"wind_force_kN": round(wind_force / 1e3, 1), "lever_m": round(z_w, 2),
                        "M_overturning_kNm": round(m_over / 1e3, 1), "M_restoring_kNm": round(m_rest / 1e3, 1),
                        "ok": bool(m_rest > m_over)},
        "weight_kN": round(weight / 1e3, 1),
        "zones": zones[:12],
        "joints": joints,
    }
    field_out = {"eidx": eidx, "util": util_env, "reason": reason}
    return out, field_out


# =========================================================================== #

def vertex_values(mesh, axes, eidx, values, shape):
    """Carry per-element values to mesh vertices (max over the 3x3x3 neighbourhood)."""
    from scipy import ndimage
    grid = np.zeros(shape, np.float32)
    grid[tuple(eidx.T)] = values
    grid = ndimage.maximum_filter(grid, size=3)
    ix, iy, iz = grid_index(axes, mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.vertices[:, 2])
    return grid[ix, iy, iz]



def main_art(args):
    """Climb, mass and structure of the terraced sculpture."""
    cfg = merge(op.DEFAULT_CONFIG, op.ART)
    if args.config:
        cfg = merge(cfg, json.loads(Path(args.config).read_text()))
    cfg = merge(cfg, {"printer": {"density_kg_m3": MATERIAL["density_kg_m3"]}})
    out = Path(args.out if args.out != str(Path(__file__).parent / "out_analysis")
               else Path(__file__).parent / "out_sculpture")
    out.mkdir(parents=True, exist_ok=True)
    R = cfg["terraces"]["riser_m"]
    log = lambda m: print(m, flush=True)

    log("field (untreaded sheet) ...")
    plain = merge(cfg, {"terraces": {"enabled": False}})
    F0, axes = op.field(plain)
    mesh0 = op.to_mesh(F0, axes, cfg["voxel_m"])
    report = {"config": cfg, "mode": "sculpture"}
    avoid, mforces = None, None     # the climb passes the room through its door and hatch
    if args.membranes and (out / "membranes.json").exists():
        mf = np.load(out / "membrane_forces.npz")
        mforces = (mf["p"], mf["f"])
        report["membranes"] = json.loads((out / "membranes.json").read_text())
    report["climb"] = climb(cfg, F0, axes, mesh0, R, log=log, avoid=avoid)

    log("field (terraced) ...")
    F, axes, extras = op.field(cfg, extras=True)
    mesh = op.to_mesh(F, axes, cfg["voxel_m"])
    mesh.export(out / "sculpture.stl")
    dens = MATERIAL["density_kg_m3"]
    v_terr = float(mesh.volume - mesh0.volume)
    report["mesh"] = {"watertight": bool(mesh.is_watertight), "volume_m3": round(float(mesh.volume), 3),
                      "sheet_volume_m3": round(float(mesh0.volume), 3), "ledge_volume_m3": round(v_terr, 3),
                      "mass_kg": round(float(mesh.volume) * dens),
                      "mass_hollow_ledges_kg": round((float(mesh0.volume) + 0.3 * v_terr) * dens),
                      "bounding_box_m": (mesh.bounds[1] - mesh.bounds[0]).round(2).tolist()}
    parts, grid = op.segments(F, axes, cfg)
    report["print"] = {"parts": len(parts), "grid": grid,
                       "heaviest_kg": round(max(m.volume for _, m in parts) * dens),
                       "build_volume_m": cfg["printer"]["build_volume_m"]}
    log(f"mesh: {report['mesh']}, {len(parts)} print parts")
    if not args.no_fe:
        st, fo = structure(cfg, F, axes, extras, [], log=log, crowd_everywhere=True, extra_forces=mforces)
        report["structure"] = st
        np.save(out / "vertex_utilisation.npy", vertex_values(mesh, axes, fo["eidx"], fo["util"], F.shape).astype(np.float32))
    np.save(out / "vertex_class.npy", ledge_classes(mesh, mesh0))
    (out / "report.json").write_text(json.dumps(report, indent=1, default=float))
    if "structure" in report:
        print(json.dumps({k: v for k, v in report["structure"].items()
                          if k in ("utilisation", "max_displacement_mm", "overturning", "weight_kN")}, default=float))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", help="organic_pavilion config (default: two-storey tower)")
    ap.add_argument("--target", default=None, help="pod to climb to, e.g. 0,0,1 (default: highest pod)")
    ap.add_argument("--no-fe", action="store_true")
    ap.add_argument("--art", action="store_true", help="the climbable terraced gyroid sculpture")
    ap.add_argument("--membranes", action="store_true", help="with the membrane-closed room (run membranes.py first)")
    ap.add_argument("--out", default=str(Path(__file__).parent / "out_analysis"))
    args = ap.parse_args()

    if args.art:
        return main_art(args)
    cfg = merge(op.DEFAULT_CONFIG, TOWER)
    if args.config:
        cfg = merge(cfg, json.loads(Path(args.config).read_text()))
    cfg = merge(cfg, {"printer": {"density_kg_m3": MATERIAL["density_kg_m3"]}})
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    print("field ...", flush=True)
    F, axes, extras = op.field(cfg, extras=True)
    mesh = op.to_mesh(F, axes, cfg["voxel_m"])
    pods = op.pods(cfg)
    target = tuple(int(v) for v in args.target.split(",")) if args.target else max(pods, key=lambda p: p[2])
    report = {"config": cfg, "rules": STAIR_RULES}

    rows = []
    if target[2] > 0:
        print("route ...", flush=True)
        route = find_route(cfg, F, axes, mesh, target, log=lambda m: print(m, flush=True))
        if route is None:
            report["stairs"] = {"feasible": False, "reason": "no walkable route to the upper floor"}
        else:
            pts, s, kind = resample_route(F, axes, route["points"], route["segments"])
            fit = fit_stair(pts, s, kind, route["floor_m"])
            F, openings = carve_openings(F, axes, fit["steps"])
            widths = [lateral_width(F, axes, st) for st in fit["steps"]]
            rows, flights, issues = check_stair(fit, widths)
            report["stairs"] = {"feasible": True, "floor_m": route["floor_m"], "length_m": round(float(s[-1]), 2),
                                **{k: v for k, v in fit.items() if k != "steps"},
                                "flights": flights, "issues": issues, "openings": openings, "steps": rows,
                                "route": np.round(pts, 3).tolist(), "route_kind": kind.tolist()}
            F = stair_solids(F, axes, rows)
            F = np.maximum(F, -axes[2][None, None, :])          # keep the flat cut at the ground
            mesh = op.to_mesh(F, axes, cfg["voxel_m"])
            # the stairs alone, for drawings and the report
            alone = stair_solids(np.full(F.shape, 1.0, np.float32), axes, rows)
            op.to_mesh(np.maximum(alone, -axes[2][None, None, :]), axes, cfg["voxel_m"]).export(out / "stairs.stl")
            print(f"stairs: {fit['risers']} risers of {fit['riser_m'] * 100:.1f} cm, tread "
                  f"{fit['tread_m'] * 100:.0f} cm, flights {flights}, issues {issues}, openings {openings}")
    mesh.export(out / "pavilion_with_stairs.stl")
    report["mesh"] = {"watertight": bool(mesh.is_watertight), "volume_m3": round(float(mesh.volume), 3),
                      "mass_kg": round(float(mesh.volume) * MATERIAL["density_kg_m3"]),
                      "bounding_box_m": (mesh.bounds[1] - mesh.bounds[0]).round(2).tolist()}

    if not args.no_fe:
        st, fo = structure(cfg, F, axes, extras, rows, log=lambda m: print(m, flush=True))
        report["structure"] = st
        util_v = vertex_values(mesh, axes, fo["eidx"], fo["util"], F.shape)
        np.save(out / "vertex_utilisation.npy", util_v.astype(np.float32))
    (out / "report.json").write_text(json.dumps(report, indent=1, default=float))
    print(json.dumps({k: v for k, v in report.get("structure", {}).items()
                      if k not in ("material", "loads", "joints", "zones")}, indent=1, default=float))


if __name__ == "__main__":
    main()
