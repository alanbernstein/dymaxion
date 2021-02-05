from math import pi, cos, radians

import numpy as np

def sphere_polyarea(lon, lat):
    # returns the area of a polygon defined in longitude,latitude,
    # via simple equal-area projection onto a unit sphere
    lat_dist = pi / 180.0

    y = [la * lat_dist for la in lat]
    x = [lo * lat_dist * cos(radians(la)) for la, lo in zip(lat, lon)]

    a = sum([x[i] * (y[i+1] - y[i-1]) for i in range(-1, len(x)-1)])
    return abs(a)/2.0


def rotation_matrix_from_src_dest_vecs(src, dest):
    # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    v = np.cross(src, dest)
    s = np.linalg.norm(v)
    c = np.dot(src, dest)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + np.dot(vx,vx) * (1-c)/(s**2)


def rotation_matrix_from_euler(**kwargs):
    # args are x,y,z, but accepted in any order, and composed in that order
    R = np.eye(3)
    for k, v in kwargs.items():
        c, s = np.cos(v), np.sin(v)
        if k == 'x':
            R = R @ [[1, 0, 0], [0, c, -s], [0, s, c]]
        if k == 'y':
            R = R @ [[c, 0, s], [0, 1, 0], [-s, 0, c]]
        if k == 'z':
            R = R @ [[c, -s, 0], [s, c, 0], [0, 0, 1]]
    return R


class DymaxionProjection(object):
    # designed to work for a regular icosahedron
    # should work for any platonic solid, any archimedean solid
    # and also another class of polyhedra that i don't know a name for:
    # those that are equivalent to spherical voronoi diagrams - that is,
    # any point on a face is closer to the central bisector of that face,
    # than to the central bisector of any other face.
    # that assumes "central bisector" is well-defined.
    # TODO: research this ^
    def __init__(self, vertices, edges, faces):
        self.vertices = vertices
        self.edges = edges
        self.faces = faces

        self.face_centers = np.array([np.mean(vertices[f], axis=0) for f in self.faces])
        self.face_center_mags = np.linalg.norm(self.face_centers, axis=1)
        self.face_unit_normals = self.face_centers / self.face_center_mags[:,None]

        self.projection = 'simple'  # 'simple', 'predistort'

    def set_projection(self, pstring):
        self.projection = pstring

    def project(self, xyz):
        if self.projection == 'simple':
            return self.project_simple(xyz)
        if self.projection == 'predistort':
            return self.project_predistort(xyz)

    def project_predistort(self, xyz):
        # plot shapes projected onto polyhedron of nonzero thickness,
        # but in such a way that after sanding the polyhedron to a sphere,
        # the shapes match what they should be for the sphere
        #
        # specifically: for each point `pt` in path,
        # 1. project onto insphere (call this `pts`)
        # 2. find intersection of
        #    a) line through pts, with direction vector = face_normal
        #    b) plane of the face


        R_ci = 1.258408572364819 # ratio of icosahedron circumradius/inradius
        # gives insphere radius relative to vertex magnitude

        pxyz = []
        best_faces = []
        for pt in xyz:
            face_axis_angles = np.arccos(np.sum(pt/np.linalg.norm(pt) * self.face_unit_normals, axis=1))
            best_face_id = np.argmin(face_axis_angles)
            best_faces.append(best_face_id)

            fc = self.face_centers[best_face_id] # arbitrary point on plane
            fn = self.face_unit_normals[best_face_id] # normal to plane

            # project pt onto insphere first...
            pts = np.linalg.norm(self.vertices[0])/R_ci * pt/np.linalg.norm(pt)

            # ...then project that point, through a line perpendicular to the face,
            # onto the face
            s = np.dot(fc-pts, fn) / np.dot(fn, fn)
            projected = pts + fn * s

            # TODO - this is missing something -
            # paths that cross polyhedron edges are no longer contiguous

            """
            (p-fc).fn = 0             # (p-p0).n = 0
            p = pt + fn*d             # p = l0 + l*d
            ((pt + fn*d)-fc).fn = 0   # ((l0+l*d)-p0).n = 0
            (pt + fn*d - fc).fn = 0
            pt.fn + d*fn.fn - fc.fn = 0
            d*fn.fn = fc.fn - pt.fn
            d = (fc-pt).fn/(fn.fn)
            p = pt + fn*d
            p = pt + fn * (fc-pt).fn/(fn.fn)
            """

            pxyz.append(projected)

        return np.array(pxyz), best_faces


    def project_simple(self, xyz):
        # for each point in the shape:
        # - find which face-line has smallest angle, select that plane
        # - use equation of plane to project point onto polyhedron.
        #   note this is a "dymaxion projection", which requires finding
        #   the intersection of a ray and a plane, /not/ the projection of
        #   that ray onto the plane, as you might mistakenly assume
        #   thanks to the overloaded term "projection".

        # specifically: for each point `pt` in path,
        # find intersection of
        # a) line through `pt` with direction vector = `pt`
        # b) the plane of the face

        pxyz = []
        best_faces = []
        for pt in xyz:
            # figure out which face this point belongs to
            face_axis_angles = np.arccos(np.sum(pt/np.linalg.norm(pt) * self.face_unit_normals, axis=1))
            best_face_id = np.argmin(face_axis_angles)
            best_faces.append(best_face_id)

            # find intersection of line (O, p) with plane with point fc and normal fn
            fc = self.face_centers[best_face_id] # arbitrary point on plane
            fn = self.face_unit_normals[best_face_id] # normal to plane
            s = np.dot(fc, fn)/np.dot(pt, fn)
            projected = pt * s

            """
            https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form
            (p-fc).fn = 0          # (p-p0).n = 0
            p = 0 + pt*d           # p = l0 + l*d
            ((pt*d)-fc).fn = 0     # ((l0+l*d)-p0).n = 0
            (pt*d).fn - fc.fn = 0
            d * pt.fn - fc.fn = 0
            d = fc.fn/pt.fn
            p = pt * d
            p = pt * fc.fn/pt.fn
            """

            """
            # project p onto the face
            vv = p - fc
            v_parallel = np.dot(vv, fn) / np.linalg.norm(fn) ** 2 * fn  # component parallel to the normal
            v_perp = vv - v_parallel  # component perpendicular to the normal
            projected = fc + v_perp
            """

            pxyz.append(projected)


        return np.array(pxyz), best_faces


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




p = (np.sqrt(5) + 1)/2  # golden ratio

icosahedron_circumradius_per_side = np.sqrt(p*np.sqrt(5))/2
icosahedron_inradius_per_side = p**2 / (2 * np.sqrt(3))
icosahedron_circumradius_per_inradius = icosahedron_circumradius_per_side / icosahedron_inradius_per_side

def icosahedron(circumradius=None, inradius=None):
    if not circumradius and not inradius:
        circumradius = 1

    if inradius and not circumradius:
        circumradius = icosahedron_circumradius_per_inradius * inradius

    vertices = np.array([
        [0, p, 1],
        [0, p, -1],
        [0, -p, 1],
        [0, -p, -1],
        [1, 0, p],
        [1, 0, -p],
        [-1, 0, p],
        [-1, 0, -p],
        [p, 1, 0],
        [p, -1, 0],
        [-p, 1, 0],
        [-p, -1, 0],
    ])
    vertices = vertices * circumradius / np.linalg.norm(vertices[0,:])
    edges = select_shortest_edges(vertices)

    faces = select_triangles(edges)

    return vertices, edges, faces


def truncated_icosahedron(circumradius=None):
    v0 = [
        [0, 1, 3*p],
        [0, 1, -3*p],
        [0, -1, 3*p],
        [0, -1, -3*p],

        [1, 2+p, 2*p],
        [1, 2+p, -2*p],
        [1, -2-p, 2*p],
        [1, -2-p, -2*p],
        [-1, 2+p, 2*p],
        [-1, 2+p, -2*p],
        [-1, -2-p, 2*p],
        [-1, -2-p, -2*p],

        [p, 2, 2*p+1],
        [p, 2, -2*p-1],
        [p, -2, 2*p+1],
        [p, -2, -2*p-1],
        [-p, 2, 2*p+1],
        [-p, 2, -2*p-1],
        [-p, -2, 2*p+1],
        [-p, -2, -2*p-1],
    ]
    # generate all even permutations of each of ^
    v1 = [[x[2], x[0], x[1]] for x in v0]
    v2 = [[x[1], x[2], x[0]] for x in v0]
    vertices = np.array(v0 + v1 + v2)

    edges = select_shortest_edges(vertices)

    faces = []
    return vertices, edges, faces


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
