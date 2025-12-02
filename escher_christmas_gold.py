"""
Project: Escher Christmas (Gold Master)
Logic: Ported from Peeps Gold Master
  - Topology Repair (Vertex Clustering)
  - Perfect 4-Color Map
  - Sawtooth Tree Shape
  - Layered Rendering (Decorations on Top)
"""

import cmath
import math
import collections
import sys
import random

# --- CONFIGURATION ---
FILENAME = "output/escher_christmas_gold.svg"
SIZE = 2000             # 2K Resolution
MAX_TILES = 12000       # Full disk
VERTEX_TOLERANCE = 0.005
BG_COLOR = "#051505"    # Deep Pine Night

# Holiday Palette
PALETTE = [
    "#27AE60", # Emerald Green
    "#C0392B", # Velvet Red
    "#F1C40F", # Gold
    "#ECF0F1"  # Snow White
]

# --- MATH ENGINE ---

def to_screen(z):
    if abs(z) >= 0.99999: z = 0.99999 * (z / abs(z))
    x = (z.real + 1) * (SIZE / 2)
    y = (1 - z.imag) * (SIZE / 2)
    return x, y

def pts_to_svg(points):
    return " ".join([f"{x:.1f},{y:.1f}" for x, y in [to_screen(p) for p in points]])

def get_rotation_map(center, angle):
    c_conj = center.conjugate()
    rot = cmath.exp(1j * angle)
    def func(z):
        num = z - center
        den = 1 - c_conj * z
        w = (num / den) * rot
        return (w + center) / (1 + c_conj * w)
    return func

# --- GEOMETRY ---

def get_base_geometry():
    # (4,3,3) Geometry
    alpha = math.pi / 4
    beta = math.pi / 3
    num = math.cos(beta) + math.cos(alpha)*math.cos(beta)
    den = math.sin(alpha)*math.sin(beta)
    r = math.tanh(math.acosh(num / den) / 2)
    # A=Top, B=Right, C=Bottom, D=Left
    return 0j, r+0j, r*cmath.exp(1j*alpha), r*1j

# --- TREE SCULPTING ---

def build_pine_edge(start, end, bumps=3, depth=0.25):
    """Generates a jagged sawtooth edge."""
    vec = end - start
    perp = vec * -1j
    pts = [start]

    for i in range(1, bumps + 1):
        # Linear progression
        t_tip  = (i - 0.5) / bumps
        t_notch = i / bumps

        # Tip sticks out
        # Taper the depth so the top branches are smaller
        local_depth = depth * (1.0 - (i/bumps)*0.2)
        tip = start + vec * t_tip + perp * local_depth
        pts.append(tip)

        # Notch cuts back in
        if i < bumps:
            notch = start + vec * t_notch + perp * 0.05
            pts.append(notch)

    pts.append(end)
    return pts

def build_tree_shape(A, B, C, D):
    # 1. Right Side (Branches) - Sawtooth
    edge_AB = build_pine_edge(A, B, bumps=3, depth=0.3)

    # 2. Left Side (Socket) - Rotated 90 deg
    edge_AD = [z * 1j for z in edge_AB]

    # 3. Bottom Right (Trunk) - Curved
    def bez(p0, p1, p2, steps=4):
        return [(1-t/steps)**2 * p0 + 2*(1-t/steps)*(t/steps) * p1 + (t/steps)**2 * p2 for t in range(steps+1)]

    vec_BC = C - B
    cp_trunk = B + vec_BC*0.5 + (vec_BC*-1j)*0.15
    edge_BC = bez(B, cp_trunk, C)

    # 4. Bottom Left (Trunk) - Mirror
    edge_DC = [1j * z.conjugate() for z in edge_BC]

    poly = []
    poly.extend(edge_AB)
    poly.extend(edge_BC[1:])
    poly.extend(edge_DC[::-1][1:])
    poly.extend(edge_AD[::-1][1:])

    # --- DECORATIONS ---
    # We define them relative to the shape so they transform with it
    decos = []

    # Star (Top Center)
    decos.append({"pos": A + (C-A)*0.15, "type": "star", "col": 2}) # Gold

    # Bulbs on Branch Tips (Right)
    # Tips are indices 1, 3, 5 in the sawtooth list
    decos.append({"pos": edge_AB[1] - (edge_AB[1]-A)*0.1, "type": "bulb", "col": 0})
    decos.append({"pos": edge_AB[3] - (edge_AB[3]-A)*0.1, "type": "bulb", "col": 1})

    # Bulbs on Branch Tips (Left)
    decos.append({"pos": edge_AD[1] - (edge_AD[1]-A)*0.1, "type": "bulb", "col": 3})
    decos.append({"pos": edge_AD[3] - (edge_AD[3]-A)*0.1, "type": "bulb", "col": 0})

    return poly, decos

# --- STEP 1: GENERATE TILES ---

def generate_tiles():
    print(f"1. Generating {MAX_TILES} Trees...")
    A, B, C, D = get_base_geometry()
    base_poly, base_decos = build_tree_shape(A, B, C, D)
    base_verts = [A, B, C, D]

    nodes = []
    spatial_map = {}
    def get_hash(z): return (int(z.real*50), int(z.imag*50))

    def add_node(poly, verts, decos):
        center = sum(poly)/len(poly)
        k = get_hash(center)
        if k in spatial_map:
            for idx in spatial_map[k]:
                if abs(nodes[idx]['c'] - center) < 0.001: return idx
        idx = len(nodes)
        nodes.append({
            'id': idx, 'c': center, 'poly': poly, 'verts': verts, 'decos': decos,
            'neighbors': set(), 'color': -1
        })
        if k not in spatial_map: spatial_map[k] = []
        spatial_map[k].append(idx)
        return idx

    # Init Center
    queue = []
    rot_90 = lambda z: z * 1j
    curr_p, curr_v, curr_d = base_poly, base_verts, base_decos

    for i in range(4):
        idx = add_node(curr_p, curr_v, curr_d)
        queue.append(idx)
        curr_p = [rot_90(p) for p in curr_p]
        curr_v = [rot_90(p) for p in curr_v]
        curr_d = [{"pos": rot_90(d["pos"]), "type": d["type"], "col": d["col"]} for d in curr_d]

    gates = [
        {"v_idx": 1, "angle": 2*math.pi/3},
        {"v_idx": 2, "angle": 2*math.pi/3},
        {"v_idx": 3, "angle": 2*math.pi/3}
    ]

    while len(nodes) < MAX_TILES and queue:
        pidx = queue.pop(0)
        parent = nodes[pidx]
        if abs(parent['c']) > 0.999: continue

        for gate in gates:
            pivot = parent['verts'][gate['v_idx']]
            rot = get_rotation_map(pivot, gate['angle'])

            nidx = add_node(
                [rot(p) for p in parent['poly']],
                [rot(v) for v in parent['verts']],
                [{"pos": rot(d["pos"]), "type": d["type"], "col": d["col"]} for d in parent['decos']]
            )
            if nidx >= len(nodes) - 1:
                queue.append(nidx)
    return nodes

# --- STEP 2: BUILD GRAPH (TOPOLOGY REPAIR) ---

def build_graph(nodes):
    print("2. Rebuilding Connectivity...")
    v_map = collections.defaultdict(list)
    for n in nodes:
        for v in n['verts']:
            v_map[(round(v.real,3), round(v.imag,3))].append(n['id'])

    connections = 0
    for n in nodes:
        potential = collections.Counter()
        for v in n['verts']:
            for oid in v_map[(round(v.real,3), round(v.imag,3))]:
                if oid != n['id']: potential[oid] += 1
        for oid, count in potential.items():
            if count >= 2: # Share edge (2 verts)
                n['neighbors'].add(oid)
                nodes[oid]['neighbors'].add(n['id'])
                connections += 1
    print(f"   Graph built. {connections//2} edges.")

# --- STEP 3: COLOR SOLVER ---

def solve_colors(nodes):
    print("3. Solving Colors...")
    for n in nodes:
        used = {nodes[nbr]['color'] for nbr in n['neighbors'] if nodes[nbr]['color'] != -1}
        c = 0
        while c in used: c += 1
        n['color'] = c

    conflicts = 0
    for n in nodes:
        for nbr in n['neighbors']:
            if n['color'] == nodes[nbr]['color']: conflicts += 1

    if conflicts == 0: print("   SUCCESS: Perfect Coloring.")
    else: print(f"   WARNING: {conflicts} Conflicts.")

# --- STEP 4: RENDER ---

def render(nodes):
    print(f"4. Rendering to {FILENAME}...")
    with open(FILENAME, 'w') as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SIZE} {SIZE}">')
        f.write(f'<circle cx="{SIZE/2}" cy="{SIZE/2}" r="{SIZE/2}" fill="{BG_COLOR}"/>')

        nodes.sort(key=lambda n: abs(n['c']))

        # Split into Visible lists for Layering
        visible = []
        for n in nodes:
            p0 = to_screen(n['poly'][0])
            pMid = to_screen(n['poly'][len(n['poly'])//2])
            size = math.sqrt((p0[0]-pMid[0])**2 + (p0[1]-pMid[1])**2)
            n['size'] = size
            if size > 0.8: visible.append(n)

        # PASS 1: Bodies
        for n in visible:
            col = PALETTE[n['color'] % len(PALETTE)]
            pts = pts_to_svg(n['poly'])
            f.write(f'<polygon points="{pts}" fill="{col}" stroke="{col}" stroke-width="3.5" stroke-linejoin="round"/>')
            f.write(f'<polygon points="{pts}" fill="none" stroke="#002200" stroke-width="1.0" stroke-linejoin="round"/>')

        # PASS 2: Decorations (On Top)
        for n in visible:
            for d in n['decos']:
                cx, cy = to_screen(d['pos'])

                # Scale logic to prevent artifacts at edge
                r_base = n['size'] * 0.08

                if d['type'] == 'star':
                    r = r_base * 1.5
                    f.write(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="#F1C40F" stroke="#FFF" stroke-width="{r*0.2:.1f}"/>')
                else:
                    # Bulb: Use a color from palette that ISN'T the tree color
                    # Offset logic: (TreeColor + BulbOffset) % 4
                    bulb_idx = (n['color'] + d['col'] + 1) % len(PALETTE)
                    col = PALETTE[bulb_idx]
                    f.write(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r_base:.1f}" fill="{col}" />')

        f.write('</svg>')
    print("Done.")

if __name__ == "__main__":
    sys.setrecursionlimit(20000)
    nodes = generate_tiles()
    build_graph(nodes)
    solve_colors(nodes)
    render(nodes)
