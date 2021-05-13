import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon

from geometry import rotation_matrix_from_src_dest_vecs

from polyhedra import (
    icosahedron,
    truncated_icosahedron,
    icosahedron_face_transform,
    truncated_icosahedron_face_transform,
)

# a "map projection class" should be capable of accepting a set of parameters, a list of lat-lon shapes,
# and returning a corresponding list of shapes in the output/projected space.
# for conventional planar projections, this is a simpler concept, because the output space is simply a 2D plane.
# dymaxion projections for CNC introduce several complications:
# - "interruptions" can cause input shapes to result in multiple output shapes,
#   so the input is a list of shapes but the output is a list of lists of shapes.
#   (https://en.wikipedia.org/wiki/Interruption_(map_projection))
# - there are two output spaces:
#   - the surface of the polyhedron
#   - the 2d plane after unfolding the polyhedron
# - the unfolding step is itself nontrivial (see polyhedra.py)
# - any shape in the 2d plane that crosses a polyhedron face boundary
#   (even if not interrupted in the 2d layout after unfolding)
#   needs to be represented as a closed loop that coincides with that boundary.
#   this is the CNC requirement, and this is what necessitates this "partial-identity azimuthal"
#   intermediate "warp" projection.

"""
DymaxionProjection provides methods to compute a variety of dymaxion-like projections,
from the sphere to a polyhedron.

init methods:
  __init__ accepts a polyhedron, either by name, or by literally specified vertices, edges, and faces
    _init_from_name and _init_from_data handle the actual initialization of these
  _init_2d calculates and caches some values that are used to transform faces, and their corresponding projected shapes, to the final 2d layout

simple API methods:
  set_projection
  project

projection implementations:

  project_predistort_90
  project_predistort_45
  project_simple_closed - returns closed-loop shapes in the XY plane, in proper unrolled-polyhedron position+orientation. should be called once for each (face, shape) pair. "simple" indicates "no predistortion"
  project_simple_archimedean_face
  project_simple_archimedean
  project_simple_platonic

azimuthal_warp_projection - helper function for project_simple_closed



"""


class DymaxionProjection(object):
    # designed to work for a regular icosahedron
    # should work for any platonic solid, any archimedean solid
    # and also another class of polyhedra that i don't know a name for:
    # those that are equivalent to spherical voronoi diagrams - that is,
    # any point on a face is closer to the central bisector of that face,
    # than to the central bisector of any other face.
    # that assumes "central bisector" is well-defined.
    # TODO: research this "voronoi polyhedron" concept
    def __init__(self, *args, **kwargs):
        if 'polyhedron' in kwargs:
            self._init_from_name(**kwargs)
        elif 'vertices' in kwargs and 'edges' in kwargs and 'faces' in kwargs:
            self._init_from_data(**kwargs)
        else:
            print('Dymaxion initialization error')
            import ipdb; ipdb.set_trace()

        self.projection = 'simple'  # 'simple', 'predistort'

    def _init_from_name(self, polyhedron, mat=None):
        # TODO: also accept an "unfolder" which might be a class satisfying a certain interface:
        if polyhedron in ['icosahedron', '20', 'icosa']:
            pv, pe, pf = icosahedron(circumradius=1)
            self.face_transform = icosahedron_face_transform
        elif polyhedron in ['truncated-icosahedron', '32', 'soccerball']:
            pv, pe, pf = truncated_icosahedron(circumradius=1)
            self.face_transform = truncated_icosahedron_face_transform
        else:
            raise ValueError

        if mat is not None:
            pv = pv @ mat

        self._init_from_data(pv, pe, pf)

        self._init_2d()

    def _init_from_data(self, vertices, edges, faces):
        # set verts, edges, faces
        # calculate simple geometry values based on v, e, f
        self.vertices = vertices
        self.edges = edges
        self.faces = faces

        self.face_centers = np.array([np.mean(vertices[f], axis=0) for f in self.faces])
        self.face_center_mags = np.linalg.norm(self.face_centers, axis=1)
        self.face_unit_normals = self.face_centers / self.face_center_mags[:,None]

    def _init_2d(self):
        # calculate some linear transform values, and cache them.
        # these are used for:
        # 1. rotating the 3d faces (and shapes projected onto them) into the XY plane
        # 2. rotating and translating the results of (1) into the proper unfolded-polyhedron position
        # note that these are currently only used internally by project_simple_closed,
        # but they can be applied externally to the 3d results of any of the other
        # projection methods. project_simple_closed depends on the 2d intersection functions
        # of the Shapely library, so it has to operate in the 2d plane, whereas the other
        # projection methods produce 3d results.
        self.face_vertices_2d = {}
        self.face_transforms_3d = {}
        self.face_transforms_2d = {}

        for face_id in range(len(self.faces)):
            fn = self.face_unit_normals[face_id]                    # face normal
            fv = self.vertices[self.faces[face_id]]                 # face vertices
            M3 = rotation_matrix_from_src_dest_vecs(fn, [0, 0, 1])  # rotation matrix to bring face into xy plane
            fv2 = fv @ M3.T                                         # 3d face vertices rotated to xy plane
            fx, fy, fr = self.face_transform(face_id, fv2)          # 2d transformation parameters
            c, s = np.cos(fr), np.sin(fr)
            M2 = np.array([[c, -s], [s, c]])                        # build 2d rotation matrix
            fv2_oriented = fv2[:, 0:2] @ M2 + [fx, fy]              # apply 2d transform to face

            self.face_vertices_2d[face_id] = fv2_oriented
            self.face_transforms_3d[face_id] = M3
            self.face_transforms_2d[face_id] = (M2, [fx, fy])

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

    def project_simple_closed(self, xyz, face_id):
        # project a shape onto a single, specified face,
        # but return a closed shape rather than an open one,
        # in the case when the shape extends beyond the face.
        #
        # this requires "clipping", i.e. computing the intersection of
        # the face polygon with the shape polygon, which is nontrivial.
        #
        # several approaches are generally possible:
        #
        # 1. project entire shape onto the face, then compute planar intersection.
        #    doesn't work because some shapes span too much of the globe, so the
        #    projection explodes
        # 2. compute intersection of spherical polygons.
        #    haven't yet found a decent library to do this, not worth
        #    implementing myself unless necessary
        # 3. compute planar intersection within some proper azimuthal projection
        #
        # 3b. use some ad-hoc azimuthal projection that works well enough

        # given a face_id, get the corresponding face and:
        # - 3d-rotate it and its corresponding shapes onto XY plane,
        # - use the face_transform function to adjust the layout within the XY plane
        # - compute intersection of face with shapes, to clip them properly.

        # given the shape path `xyz`:
        # - project it from its sphere-surface xyz 3d points to an intermediate projection
        #     this retains the shape as is, in the vicinity of the face, but condenses the rest of it,
        #     so that the polyhedron-face projection doesn't blow up
        # - project the intermediate projection onto the (single known face of the) polyhedron
        # - transform the fully-projected shape to the xy plane, in correct poly-net orientation
        # - compute intersection of shape and face


        projected = {}

        fn = self.face_unit_normals[face_id]             # face normal
        M3 = self.face_transforms_3d[face_id]            # rotation matrix to bring face into xy plane
        M2, dxy = self.face_transforms_2d[face_id]       # transform for bringing face/shape to final position+orientation
        fv2 = self.face_vertices_2d[face_id]
        poly_face = Polygon(fv2)

        def poly_intersection(face, shape):
            geometry_error = False

            clipped = None
            try:
                clipped = shape.intersection(face)
            except Exception as exc:
                print('    invalid geometry')
                geometry_error = True

            paths = []
            if type(clipped) == Polygon:
                # convert single polygon to multipolygon
                clipped = MultiPolygon([clipped])

            if type(clipped) == MultiPolygon:
                for poly in clipped:
                    paths.append(np.vstack(poly.exterior.coords.xy).T)

            return paths, geometry_error


        # the same sequence of operations is done on both the original shape and the "warped" shape
        # because geometry flaws can break either (or both) of the resulting intersections
        proj3d = self.project_simple_archimedean_face(xyz, face_id)
        proj2d = proj3d @ M3.T
        proj2d_oriented = proj2d[:, 0:2] @ M2 + dxy
        poly_shape = Polygon(proj2d_oriented)


        warped = azimuthal_warp_projection(xyz, fn)

        proj3d_warped = self.project_simple_archimedean_face(warped, face_id)
        proj2d_warped = proj3d_warped @ M3.T
        proj2d_warped_oriented = proj2d_warped[:, 0:2] @ M2 + dxy
        poly_shape_warped = Polygon(proj2d_warped_oriented)

        projected['unwarped'], gerror_unwarped = poly_intersection(poly_face, poly_shape)
        projected['warped'], gerror_warped = poly_intersection(poly_face, poly_shape_warped)

        if not gerror_unwarped:
            projected['final'] = projected['unwarped']
        elif not gerror_warped:
            projected['final'] = projected['warped']
        else:
            projected['final'] = []
            # known to happen on face 5, with antarctica's "main" section
            import matplotlib.pyplot as plt
            plt.plot(proj2d_oriented[:,0], proj2d_oriented[:,1],'r')
            plt.plot(fv2[:,0], fv2[:,1],'g')
            import ipdb; ipdb.set_trace()

        return projected

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
        print('WARNING: project_predistort_90 is not implemented generically') # TODO

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

            # TODO - why are paths that cross polyhedron edges no longer contiguous
            # that might be expected... hard to be sure

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

    def project_simple_archimedean_face(self, xyz, face_id):
        # project a shape entirely onto a single, specified face.
        # used when computing the proper 2D intersection of the shape
        # with the face polygon - only the points that get projected onto
        # the face should result from that intersection

        pxyz = []
        fc = self.face_centers[face_id] # arbitrary point on plane
        fn = self.face_unit_normals[face_id] # normal to plane
        num = np.dot(fc, fn)
        for pt in xyz:
            s = num/np.dot(pt, fn)
            # TODO if too far away from the face, then...?
            pxyz.append(pt * s)

        return np.array(pxyz)

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


def azimuthal_warp_projection(xyz, c):
    # project Nx3 shape `xyz` to the plane with normal `c`, (TODO: what plane? what projection?)
    # using a goofy azimuthal projection that is identity under some limit,
    # and asymptotically approaches pi/2 above the limit. this maps the
    # entire sphere to the hemisphere centered on `c`, which simplifies computation
    # of intersections of oversized spherical polyhedra.

    x0, y0, y1 = np.pi * 0.14, np.pi * 0.14, np.pi*0.5  # limit for truncated icosahedron is tan(1/2.478) = pi*0.1350
    h, k, m = x0+y0-y1, y1, -(y1-y0)**2
    f = lambda x: (m + k*(x-h))/(x-h)
    # this is a piecewise, asymptotic function designed such that
    # f(x) = x for x < x0    (this case is handled outside of the lambda)
    # f(x0) = y0
    # f'(x0) = 1
    # f(inf) -> y1

    proj = []
    for pt in xyz:
        if np.linalg.norm(pt) < 1e-15:
            # dumb glitch
            continue

        # TODO: comment this
        a = np.arccos(np.dot(pt, c)/(np.linalg.norm(pt)*np.linalg.norm(c)))

        # TODO: comment this, separate into function
        # https://en.wikipedia.org/wiki/Slerp
        omega = np.arccos(np.dot(pt, c))
        #if omega == 0:
        #    db()
        t = 1
        if a >= x0:
            t = f(a)/a
        v = np.sin((1-t)*omega)/np.sin(omega) * c + np.sin(t*omega)/np.sin(omega) * pt

        # TODO: comment this
        axis = np.cross(c, pt)
        proj.append(v)

    return np.array(proj)
