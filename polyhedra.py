from math import pi, sqrt

import numpy as np


def cube():
    vertices = np.array([
        [-1, -1, -1],
        [-1, -1, 1],
        [-1, 1, 1],
        [-1, 1, -1],
        [1, -1, -1],
        [1, -1, 1],
        [1, 1, 1],
        [1, 1, -1],
    ])
    edges = select_shortest_edges(vertices)
    faces = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 1, 4, 5],
        [2, 3, 6, 7],
        [0, 3, 4, 7],
        [1, 2, 5, 6],
    ]
    return vertices, edges, faces


phi = (sqrt(5) + 1)/2  # golden ratio

icosahedron_circumradius_per_side = sqrt(phi*sqrt(5))/2
icosahedron_inradius_per_side = phi**2 / (2 * sqrt(3))
icosahedron_circumradius_per_inradius = icosahedron_circumradius_per_side / icosahedron_inradius_per_side

def icosahedron(circumradius=None, inradius=None):
    if not circumradius and not inradius:
        circumradius = 1

    if inradius and not circumradius:
        circumradius = icosahedron_circumradius_per_inradius * inradius

    vertices = np.array([
        [0, phi, 1],
        [0, phi, -1],
        [0, -phi, 1],
        [0, -phi, -1],
        [1, 0, phi],
        [1, 0, -phi],
        [-1, 0, phi],
        [-1, 0, -phi],
        [phi, 1, 0],
        [phi, -1, 0],
        [-phi, 1, 0],
        [-phi, -1, 0],
    ])
    vertices = vertices * circumradius / np.linalg.norm(vertices[0,:])
    edges = select_shortest_edges(vertices)

    faces = select_triangles(edges)

    return vertices, edges, faces


def icosahedron_face_transform(fid, verts):
    # assumes a regular polygon, centered at the origin
    #
    # given a face id and the 2d vertex positions,
    # compute 2d translation and rotations necessary to
    # move the face, and any shapes it contains, into the
    # appropriate 2d position to produce a polyhedron net.
    #
    # this is a crappy ad-hoc solution, and it's tied to the specific
    # way that the icosahedron is defined, but it's good enough for now.
    #
    R = np.linalg.norm(verts[0,:])
    L = R/icosahedron_circumradius_per_side  # side length

    p3 = pi/3
    r3_6 = sqrt(3)/6

    triangle_inradius = L/2 * np.tan(pi/6)
    triangle_circumradius = L/(2*np.cos(pi/6))
    triangle_width = L
    triangle_height = L * sqrt(3)/2

    x0, y0 = L, triangle_circumradius  # adjust to tight fit in quadrant I

    transmap = {
        # face_id: [x, y, angle]
        # for clarity, angle is decomposed into
        # face_alignment_angle + shape_alignment_angle
        # with explicit zero values.
        # note that faces with no shape content do not necessarily
        # have the correct shape_alignment_angle.
        # south cap faces
        19: [0, 0, 1*p3 + 2*p3],
        15: [1, 0, 1*p3 + 0*p3],
        13: [2, 0, 1*p3 + 4*p3],
        5:  [3, 0, 1*p3 + 0*p3],
        7:  [4, 0, 1*p3 + 4*p3],
        # south equator faces
        18: [0, 2*r3_6, 0*p3 + 0*p3],
        9:  [1, 2*r3_6, 0*p3 + 0*p3],
        14: [2, 2*r3_6, 0*p3 + 4*p3],
        6:  [3, 2*r3_6, 0*p3 + 2*p3],
        1:  [4, 2*r3_6, 0*p3 + 0*p3],
        # north equator faces
        4:  [0-.5, 3*r3_6, 1*p3 + 4*p3],
        12: [1-.5, 3*r3_6, 1*p3 + 0*p3],
        8:  [2-.5, 3*r3_6, 1*p3 + 4*p3],
        17: [3-.5, 3*r3_6, 1*p3 + 2*p3],
        0:  [4-.5, 3*r3_6, 1*p3 + 0*p3],
        # north cap faces
        2:  [0-.5, 5*r3_6, 0*p3 + 4*p3],
        10: [1-.5, 5*r3_6, 0*p3 + 2*p3],
        11: [2-.5, 5*r3_6, 0*p3 + 4*p3],
        16: [3-.5, 5*r3_6, 0*p3 + 0*p3],
        3:  [4-.5, 5*r3_6, 0*p3 + 2*p3],
    }

    verts = verts - np.mean(verts, axis=0)
    angle = np.arctan2(verts[0, 0], verts[0, 1])
    x, y, a = transmap[fid]
    return x0 + L * x, y0 + L * y, -angle+a


truncated_icosahedron_circumradius_per_side = 1/2 * sqrt(1 + 9 * phi**2)  # 2.478

def truncated_icosahedron(circumradius=None, inradius=None):
    if not circumradius and not inradius:
        circumradius = 1

    if inradius and not circumradius:
        circumradius = truncated_icosahedron_circumradius_per_inradius * inradius

    v0 = [
        [0, 1, 3*phi],
        [0, 1, -3*phi],
        [0, -1, 3*phi],
        [0, -1, -3*phi],

        [1, 2+phi, 2*phi],
        [1, 2+phi, -2*phi],
        [1, -2-phi, 2*phi],
        [1, -2-phi, -2*phi],
        [-1, 2+phi, 2*phi],
        [-1, 2+phi, -2*phi],
        [-1, -2-phi, 2*phi],
        [-1, -2-phi, -2*phi],

        [phi, 2, 2*phi+1],
        [phi, 2, -2*phi-1],
        [phi, -2, 2*phi+1],
        [phi, -2, -2*phi-1],
        [-phi, 2, 2*phi+1],
        [-phi, 2, -2*phi-1],
        [-phi, -2, 2*phi+1],
        [-phi, -2, -2*phi-1],
    ]
    # generate all even permutations of each of ^
    v1 = [[x[2], x[0], x[1]] for x in v0]
    v2 = [[x[1], x[2], x[0]] for x in v0]
    vertices = np.array(v0 + v1 + v2)
    vertices = vertices * circumradius / np.linalg.norm(vertices[0,:])

    edges = select_shortest_edges(vertices)

    faces = select_planar_sets(vertices, edges)
    return vertices, edges, faces


class Unfolder(object):
    def __init__(self):
        pass

    def unfold(self, xyz, face_id):
        M = np.eye(3)
        dxy = [0, 0]
        xyz_transformed = (xyz @ M)[:,0:2] + dxy
        return xyz_transformed


class GridUnfolder(Unfolder):
    def __init__(self):
        pass

    def unfold(self, xyz, face_id):
        M = np.eye(3)
        dxy = [0, 0]
        xyz_transformed = (xyz @ M)[:,0:2] + dxy
        return xyz_transformed


class HardcodedUnfolderU22(Unfolder):
    name = 'icosahedron'
    fullname = 'convex-regular-icosahedron'
    num_faces = 20

    def __init__(self):
        pass

    def unfold(self, xyz, face_id):
        # TODO move face_transform functions here...
        M = np.eye(3)
        dxy = [0, 0]
        xyz_transformed = (xyz @ M)[:,0:2] + dxy
        return xyz_transformed


class HardcodedUnfolderU25(Unfolder):
    name = 'truncated-icosahedron'
    fullname = 'truncated-icosahedron'
    num_faces = 32

    def __init__(self):
        pass

    def unfold(self, xyz, face_id):
        # TODO move face_transform functions here...
        M = np.eye(3)
        dxy = [0, 0]
        xyz_transformed = (xyz @ M)[:,0:2] + dxy
        return xyz_transformed


class TaggedEdgeUnfolder(Unfolder):
    # fully generic unfolder
    # TODO:
    # tag each edge as eitehr a "cut" or a "fold", and write some
    # code to implement the "unfolding" based on these tags.
    # http://philogb.github.io/page/myriahedral/
    def __init__(self):
        pass

    def unfold(self, xyz, face_id, tagged_edge_list):
        pass


def truncated_icosahedron_face_transform(fid, verts):
    # assumes a regular polygon, centered at the origin
    #
    # given a face id and the 2d vertex positions,
    # compute 2d translation and rotations necessary to
    # move the face, and any shapes it contains, into the
    # appropriate 2d position to produce a polyhedron net.
    #
    # this is a crappy ad-hoc solution, and it's tied to the specific
    # way that the icosahedron is defined, but it's good enough for now.
    R = np.linalg.norm(verts[0,:])
    L = R/truncated_icosahedron_circumradius_per_side  # side length

    p5 = pi/5
    p3 = pi/3
    inr5 = 1/10 * sqrt(25 + 10*sqrt(5))# inradius of pentagon
    r3_2 = sqrt(3)/2
    pentagon_inradius = L / (2*np.tan(p5))
    pentagon_circumradius = L / (2*np.sin(p5))
    pentagon_height = pentagon_inradius + pentagon_circumradius
    hexagon_width = L*2
    hexagon_height = L*2*r3_2

    x0, y0 = hexagon_width/2, hexagon_height/2 + pentagon_height  # adjust to tight fit in quadrant I

    transmap = {
        # face_id: [x, y, angle]
        # note that faces with no shape content do not necessarily
        # have the correct shape alignment angle.

        # south pole pentagon
        5:  [0, -r3_2-inr5, p5 + 4*p5],

        # antarctic hexagons
        10: [ 0, 0, 4*p3],
        17: [ 3, 0, 4*p3],
        3:  [ 6, 0, 0*p3],
        4:  [ 9, 0, 4*p3],
        11: [12, 0, 3*p3],

        # capricorn pentagons
        30: [ 1.5, 2*r3_2-inr5, p5 + 6*p5],
        27: [ 4.5, 2*r3_2-inr5, p5 + 0*p5],
        7:  [ 7.5, 2*r3_2-inr5, p5 + 4*p5],
        26: [10.5, 2*r3_2-inr5, p5 + 0*p5],
        28: [13.5, 2*r3_2-inr5, p5 + 2*p5],

        # south equator hexagons
        9:  [ 0, 2*r3_2, 0*p3],
        24: [ 3, 2*r3_2, 1*p3],
        19: [ 6, 2*r3_2, 0*p3],
        14: [ 9, 2*r3_2, 5*p3],
        20: [12, 2*r3_2, 0*p3],

        # north equator hexagons
        16: [ 1.5, 3*r3_2, 5*p3],
        23: [ 4.5, 3*r3_2, 0*p3],
        15: [ 7.5, 3*r3_2, 3*p3],
        22: [10.5, 3*r3_2, 2*p3],
        8:  [13.5, 3*r3_2, 2*p3],

        # cancer pentagons
        1:  [0 , 3*r3_2 + inr5, 0*p5],
        25: [3 , 3*r3_2 + inr5, 0*p5],
        31: [6 , 3*r3_2 + inr5, 2*p5],
        29: [9 , 3*r3_2 + inr5, 8*p5],
        21: [12, 3*r3_2 + inr5, 6*p5],

        # arctic hexagons
        2:  [ 1.5, 5*r3_2, 5*p3],
        18: [ 4.5, 5*r3_2, 1*p3],
        12: [ 7.5, 5*r3_2, 1*p3],
        13: [10.5, 5*r3_2, 0*p3],
        0:  [13.5, 5*r3_2, 2*p3],

        # north pole pentagon
        6:  [1.5, 6*r3_2 + inr5, 6*p5],
    }
    verts = verts - np.mean(verts, axis=0)
    angle = np.arctan2(verts[0, 0], verts[0, 1])
    if len(verts) == 6:
        angle += np.pi/6

    x, y, a = transmap[fid]
    return x0 + L * x, y0 + L * y, -angle+a


def select_planar_sets(v, edges):
    # given vertices v, edges e
    # find all planar sets of vertices (faces)

    r_edges = [[x[1], x[0]] for x in edges]  # reversed edges
    d_edges = edges + r_edges  # directed edges

    faces = []

    while d_edges:  # face loop
        # print('face #%d' % (len(faces)))
        face = list(d_edges[0])  # initialize a face with an arbitrarily chosen item from the list of directed edges
        # print('  edge 0: 0: [%d, %d]' % (face[0], face[1]))
        del(d_edges[0])  # remote it from the list
        next = False

        while True:  # edge loop
            # find all other edges that start at the end of the current face path (except the backtrack)
            ei = [k for k, edge in enumerate(d_edges) if edge[0] == face[-1] and edge[1] != face[-2]]

            # choose the clockwisest edge to follow.
            # 1. compute the cross product of the previous edge and the next edge
            # 2. compute the dot product of that and the "vertex normal"
            # 3. choose the edge with the biggest dot product
            # this might only work for polyhedra with exactly 3 edges per node.
            choices = []
            for n, eid in enumerate(ei):
                cc = np.cross(v[face[-2]]-v[face[-1]], v[d_edges[eid][0]]-v[d_edges[eid][1]])
                mm = np.dot(cc, v[face[-1]])
                turn_angle = 0
                # TODO instead of this dot product,  find the turn angle, in a coordinate system where:
                # 1. face[-2] -> face[-1] is the x axis
                # 2. vertex normal is the z axis
                choices.append((n, mm))

            choices.sort(key=lambda x: -x[1])
            n = choices[0][0]

            face.append(d_edges[ei[n]][1])  # add it to the face

            # print('  edge %d: %d: [%d, %d]' % (len(face)-2, n, face[-2], face[-1]))
            if len(face) > 7:
                import ipdb; ipdb.set_trace()

            del(d_edges[ei[n]])  # remove it from the list
            if face[0] == face[-1]:
                # face cycle is complete, move on
                faces.append(face[0:-1])
                # print('  done: %s' % face)
                break

    return faces


def coplanar(v, e1, e2):
    # p = p1 + x*d1
    # p = p2 + y*d2
    # (p1-p2) . (d1 x d2) = 0  # coplanar condition

    v11, v12 = v[e1[0]], v[e1[1]]
    v21, v22 = v[e2[0]], v[e2[1]]
    z = np.dot(v11 - v21, np.cross(v12-v11, v22-v21))
    return abs(z) < 1e-6


def select_smallest_faces(edges, L):
    # TODO: will need this for any regular(ish) polyhedra with non-triangular faces
    # generalization of select_triangles
    # - find all cycles of length L
    # - find perimeter of each cycle
    # - find the minimum perimeter
    # - select all cycles with minimum perimeter (plus epsilon)
    pass

def select_triangles(edges):
    tris = []
    for i, e1 in enumerate(edges[0:(len(edges)-1)]):
        for e2 in edges[(i+1):len(edges)]:
            if e1[0] == e2[0] and (e1[1], e2[1]) in edges:
                tris.append([e1[0], e1[1], e2[1]])
    return tris


def select_shortest_edges(vertices):
    def dist2(x, y):
        return sum([float(a-b)**2 for a, b in zip(x, y)])

    # compute all edge lengths, and the min length
    N = len(vertices)
    all_edges = []
    min_d = 10000
    for n1 in range(N-1):
        for n2 in range(n1+1, N):
            d = dist2(vertices[n1], vertices[n2])
            all_edges.append((n1, n2, d))
            if d < min_d:
                min_d = d

    # select all edges that are floating-point-equal to the min
    thresh = min_d + 0.000001
    valid_edges = [e[0:2] for e in all_edges if e[2] <= thresh]
    return valid_edges
