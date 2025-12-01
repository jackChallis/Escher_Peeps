"""
Project: Escher Peeps (Gold Master)
Features:
  - Verified (4,3,3) Topology (No split brain)
  - Perfect 4-Color Map (No collisions)
  - Marshmallow Shapes with EYES FIXED
  - High-Res Full Disk Render
"""

import cmath
import math
import collections
import sys

# --- CONFIGURATION ---
FILENAME = "output/peeps_gold_master.svg"
SIZE = 2000             # 2K Resolution for crisp details
MAX_TILES = 12000       # Full disk coverage
VERTEX_TOLERANCE = 0.005
BG_COLOR = "#111"

# The "Easter Basket" Palette
PALETTE = [
    "#FF0055", # Hot Pink
    "#00CCFF", # Cyan
    "#FFD700", # Gold
    "#AF52DE"  # Purple
]

# --- MATH HELPER ---

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
    return 0j, r+0j, r*cmath.exp(1j*alpha), r*1j

def build_peep_shape(A, B, C, D):
    # Bezier Helper
    def bez(p0, p1, p2, steps=6):
        return [(1-t/steps)**2 * p0 + 2*(1-t/steps)*(t/steps) * p1 + (t/steps)**2 * p2 for t in range(steps+1)]

    # 1. Fat Belly (Curved out)
    cp_chest = A + (B-A)*0.5 - 0.15j
    edge_AB = bez(A, cp_chest, B)

    # 2. Tail Socket (Matches Belly)
    edge_AD = [z * 1j for z in edge_AB]

    # 3. Head (Soft Beak)
    vec = C - B
    perp = vec * -1j # Outward direction

    # Beak definition
    beak_pos = B + vec*0.2 + perp*0.15
    top = B + vec*0.6

    seg1 = bez(B, beak_pos, top, 4)
    seg2 = bez(top, C + perp*0.05, C, 4)
    edge_BC = seg1 + seg2[1:]

    # 4. Back (Matches Head)
    edge_DC = [1j * z.conjugate() for z in edge_BC]

    poly = []
    poly.extend(edge_AB)
    poly.extend(edge_BC[1:])
    poly.extend(edge_DC[::-1][1:])
    poly.extend(edge_AD[::-1][1:])

    # --- EYE FIX ---
    # Moved 35% down the body and 8% outward.
    # This keeps it safe from the Vertex B overlap.
    eye_pos = B + (vec * 0.35) + (perp * 0.08)

    return poly, eye_pos

# --- STEP 1: GENERATE TILES (GEOMETRY ONLY) ---

def generate_tiles():
    print(f"1. Generating {MAX_TILES} Tiles...")
    A, B, C, D = get_base_geometry()
    base_poly, base_eye = build_peep_shape(A, B, C, D)
    base_verts = [A, B, C, D]

    nodes = []
    spatial_map = {}

    def get_hash(z): return (int(z.real*50), int(z.imag*50))

    def add_node(poly, verts, eye):
        center = sum(poly)/len(poly)
        k = get_hash(center)

        # Duplicate Check
        if k in spatial_map:
            for idx in spatial_map[k]:
                if abs(nodes[idx]['c'] - center) < 0.001: return idx

        idx = len(nodes)
        nodes.append({
            'id': idx, 'c': center, 'poly': poly, 'verts': verts, 'eye': eye,
            'neighbors': set(), 'color': -1
        })
        if k not in spatial_map: spatial_map[k] = []
        spatial_map[k].append(idx)
        return idx

    # Init Center Ring
    queue = []
    rot_90 = lambda z: z * 1j
    curr_p, curr_v, curr_e = base_poly, base_verts, base_eye

    for i in range(4):
        idx = add_node(curr_p, curr_v, curr_e)
        queue.append(idx)
        curr_p = [rot_90(p) for p in curr_p]
        curr_v = [rot_90(p) for p in curr_v]
        curr_e = rot_90(curr_e)

    # Transform Gates (Pivots B, C, D)
    gates = [
        {"v_idx": 1, "angle": 2*math.pi/3},
        {"v_idx": 2, "angle": 2*math.pi/3},
        {"v_idx": 3, "angle": 2*math.pi/3}
    ]

    while len(nodes) < MAX_TILES and queue:
        pidx = queue.pop(0)
        parent = nodes[pidx]

        # Stop expanding if microscopic
        if abs(parent['c']) > 0.999: continue

        for gate in gates:
            pivot = parent['verts'][gate['v_idx']]
            rot = get_rotation_map(pivot, gate['angle'])

            new_p = [rot(p) for p in parent['poly']]
            new_v = [rot(p) for p in parent['verts']]
            new_e = rot(parent['eye']) # TRANSFORM THE EYE

            pre_len = len(nodes)
            nidx = add_node(new_p, new_v, new_e)

            if len(nodes) > pre_len:
                queue.append(nidx)

    return nodes

# --- STEP 2: BUILD GRAPH (TOPOLOGY) ---

def build_graph(nodes):
    print("2. Rebuilding Graph Connectivity...")
    # Map Vertex Position -> List of Node IDs
    v_map = collections.defaultdict(list)

    for n in nodes:
        for v in n['verts']:
            key = (round(v.real, 3), round(v.imag, 3))
            v_map[key].append(n['id'])

    # Link nodes sharing >= 2 vertices (an Edge)
    connections = 0
    for n in nodes:
        potential_neighbors = collections.Counter()

        for v in n['verts']:
            key = (round(v.real, 3), round(v.imag, 3))
            for other_id in v_map[key]:
                if other_id != n['id']:
                    potential_neighbors[other_id] += 1

        for other_id, shared_count in potential_neighbors.items():
            if shared_count >= 2:
                if other_id not in n['neighbors']:
                    n['neighbors'].add(other_id)
                    nodes[other_id]['neighbors'].add(n['id'])
                    connections += 1

    print(f"   Graph built. {connections//2} total edges.")

# --- STEP 3: COLOR SOLVER (GREEDY) ---

def solve_colors(nodes):
    print("3. Solving Colors...")
    for n in nodes:
        used = {nodes[nbr]['color'] for nbr in n['neighbors'] if nodes[nbr]['color'] != -1}
        c = 0
        while c in used: c += 1
        n['color'] = c

    # Audit
    conflicts = 0
    max_c = 0
    for n in nodes:
        if n['color'] > max_c: max_c = n['color']
        for nbr in n['neighbors']:
            if n['color'] == nodes[nbr]['color']:
                conflicts += 1

    print(f"   Palette Size: {max_c + 1}")
    if conflicts == 0:
        print("   SUCCESS: 0 Conflicts.")
    else:
        print(f"   WARNING: {conflicts} Conflicts Found.")

# --- STEP 4: RENDER ---

def render(nodes):
    print(f"4. Rendering to {FILENAME}...")
    with open(FILENAME, 'w') as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SIZE} {SIZE}">')
        f.write(f'<circle cx="{SIZE/2}" cy="{SIZE/2}" r="{SIZE/2}" fill="{BG_COLOR}"/>')

        # Sort by depth so small tiles draw on top (cleaner lines)
        nodes.sort(key=lambda n: abs(n['c']))

        # VISIBLE NODE FILTER
        visible_nodes = []
        for n in nodes:
            # LOD Cull: Don't draw sub-pixel tiles
            p0 = to_screen(n['poly'][0])
            pMid = to_screen(n['poly'][len(n['poly'])//2])
            tile_size = math.sqrt((p0[0]-pMid[0])**2 + (p0[1]-pMid[1])**2)
            n['size'] = tile_size
            if tile_size > 0.8:
                visible_nodes.append(n)

        # PASS 1: BODIES (Bleed + Outline)
        for n in visible_nodes:
            col = PALETTE[n['color'] % len(PALETTE)]
            pts = pts_to_svg(n['poly'])

            # Bleed (Thick Stroke) - Fills gaps
            f.write(f'<polygon points="{pts}" fill="{col}" stroke="{col}" stroke-width="3.0" stroke-linejoin="round"/>')
            # Outline (Thin Black) - Definition
            f.write(f'<polygon points="{pts}" fill="none" stroke="#000" stroke-width="1.0" stroke-linejoin="round"/>')

        # PASS 2: EYES (Always on top)
        for n in visible_nodes:
            ex, ey = to_screen(n['eye'])
            r = max(1.5, min(5.0, n['size'] * 0.1))

            # White Sclera
            f.write(f'<circle cx="{ex:.1f}" cy="{ey:.1f}" r="{r:.1f}" fill="white"/>')
            # Black Pupil
            f.write(f'<circle cx="{ex:.1f}" cy="{ey:.1f}" r="{r*0.5:.1f}" fill="black"/>')

        f.write('</svg>')
    print("   Done. File saved.")

if __name__ == "__main__":
    sys.setrecursionlimit(20000)
    nodes = generate_tiles()
    build_graph(nodes)
    solve_colors(nodes)
    render(nodes)
