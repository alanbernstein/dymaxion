import numpy as np


class DymaxionProjection(object):
    # designed to work for a regular icosahedron
    # should work for any platonic solid, any archimedean solid
    # and also another class of polyhedra that i don't know a name for:
    # those that are equivalent to spherical voronoi diagrams - that is,
    # any point on a face is closer to the central bisector of that face,
    # than to the central bisector of any other face.
    # that assumes "central bisector" is well-defined.
    # TODO: research this "voronoi polyhedron" concept
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
