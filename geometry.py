from math import pi, cos, radians, sqrt

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
            if len(self.faces) in [4, 6, 8, 12, 20]:
                return self.project_simple_platonic(xyz)
            else:
                return self.project_simple_archimedean(xyz)
        if self.projection == 'predistort-90':
            return self.project_predistort_90(xyz)
        if self.projection == 'predistort-45':
            return self.project_predistort_45(xyz)

    def project_predistort_90(self, xyz):
        # plot shapes projected onto polyhedron of nonzero thickness,
        # but in such a way that after sanding the polyhedron to a sphere,
        # the shapes match what they should be for the sphere.
        # the "90" refers to using an end mill (90 degrees from the plane)
        # to cut out the shapes, which means the extrusion of the shapes
        # through the thick face is perpendicular to the face.
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

    def project_predistort_45(self, xyz):
        # similar concept to predistort_90, but here the "45" refers
        # to using a v-groove bit with a 45 angle off the plane,
        # so the exxtrusion of the shapes through the thick face is
        # at a 45-degree angle. should also be able to generalize this
        # to arbitrary other angles (but really just 30, 60)
        raise NotImplementedError

    def project_simple_archimedean(self, xyz):
        # for each point in the shape:
        # - select corresponding face
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
        k = 2
        for pt in xyz:
            # figure out which face this point belongs to.
            #
            # instead of just choosing the smallest face_axis_angle, look at the
            # two smallest, find the corresponding point-on-face for both, and
            # choose the one that is closer to the centroid of the polyhedron.
            #
            # this is a quick and dirty solution to the problem that the method
            # of just choosing the smallest angle only works for platonic solids.
            # there may be a way that is more "correct", but this seems like it will
            # be about as fast as we can hope for?
            #
            # another idea: find a multiplicative factor to apply to the computed
            # angles, based on the size of the face. not sure if this makes sense.
            #
            # see project_simple_platonic for geometry explanation

            face_axis_angles = np.arccos(np.sum(pt/np.linalg.norm(pt) * self.face_unit_normals, axis=1))

            # argpartition sorts the input array, but only enough that the
            # first k values are in ascending order, which is all we need
            idxs = np.argpartition(face_axis_angles, k)

            min_rad = 10 * np.linalg.norm(self.vertices[0])
            min_projected = []

            for i in range(k):
                fc = self.face_centers[idxs[i]] # arbitrary point on plane
                fn = self.face_unit_normals[idxs[i]] # normal to plane
                s = np.dot(fc, fn)/np.dot(pt, fn)
                projected = pt * s
                rad = np.linalg.norm(projected)
                if rad < min_rad:
                    min_rad = rad
                    min_projected = projected
                    best_face_id = idxs[i]

            best_faces.append(best_face_id)
            pxyz.append(min_projected)


        return np.array(pxyz), best_faces


    def project_simple_platonic(self, xyz):
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
    # given a face id and the 2d vertex positions,
    # compute 2d translation and rotations necessary to
    # move the face, and any shapes it contains, into the
    # appropriate 2d position to produce a polyhedron net.
    #
    # this is a crappy ad-hoc solution, and it's tied to the specific
    # way that the icosahedron is defined, but it's good enough for now.
    #
    # TODO: probably want to automate this for larger polyhedra.
    R = np.linalg.norm(verts[0,:])
    L = R/icosahedron_circumradius_per_side  # side length

    p3 = pi/3
    r3_6 = sqrt(3)/6

    x0, y0 = 11/L, 7/L  # SVG doesn't like negative coordinates

    transmap = {
        # face_id: [x, y, angle]
        # for clarity, angle is decomposed into
        # face_alignment_angle + shape_alignment_angle
        # with explicit zero values.
        # south cap faces
        19: [x0 + 0, y0 + 0, 1*p3 + 2*p3],
        15: [x0 + 1, y0 + 0, 1*p3 + 0*p3],
        13: [x0 + 2, y0 + 0, 1*p3 + 4*p3],
        5:  [x0 + 3, y0 + 0, 1*p3 + 0*p3],
        7:  [x0 + 4, y0 + 0, 1*p3 + 4*p3],
        # south equator faces
        18: [x0 + 0, y0 + 2*r3_6, 0*p3 + 0*p3],
        9:  [x0 + 1, y0 + 2*r3_6, 0*p3 + 0*p3],
        14: [x0 + 2, y0 + 2*r3_6, 0*p3 + 4*p3],
        6:  [x0 + 3, y0 + 2*r3_6, 0*p3 + 2*p3],
        1:  [x0 + 4, y0 + 2*r3_6, 0*p3 + 0*p3],
        # north equator faces
        4:  [x0 + 0-.5, y0 + 3*r3_6, 1*p3 + 4*p3],
        12: [x0 + 1-.5, y0 + 3*r3_6, 1*p3 + 0*p3],
        8:  [x0 + 2-.5, y0 + 3*r3_6, 1*p3 + 4*p3],
        17: [x0 + 3-.5, y0 + 3*r3_6, 1*p3 + 2*p3],
        0:  [x0 + 4-.5, y0 + 3*r3_6, 1*p3 + 0*p3],
        # north cap faces
        2:  [x0 + 0-.5, y0 + 5*r3_6, 0*p3 + 4*p3],
        10: [x0 + 1-.5, y0 + 5*r3_6, 0*p3 + 2*p3],
        11: [x0 + 2-.5, y0 + 5*r3_6, 0*p3 + 4*p3],
        16: [x0 + 3-.5, y0 + 5*r3_6, 0*p3 + 0*p3],
        3:  [x0 + 4-.5, y0 + 5*r3_6, 0*p3 + 2*p3],
    }

    verts = verts - np.mean(verts, axis=0)
    angle = np.arctan2(verts[0, 0], verts[0, 1])
    x, y, a = transmap[fid]
    return L * x, L * y, -angle+a


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


def truncated_icosahedron_face_transform(fid, verts):
    R = np.linalg.norm(verts[0,:])
    L = R/truncated_icosahedron_circumradius_per_side  # side length

    x = 2 * L * (fid % 8)
    y = 2 * L * (fid // 8)

    # return x, y, angle
    return x, y, 0


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
