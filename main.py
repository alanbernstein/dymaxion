from collections import defaultdict
import json
from pprint import pprint as pp
import sys
from ipdb import set_trace as db
from debugx import pm

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors as mcolors
import numpy as np
from shapely.geometry import Point, Polygon, MultiPolygon
from spherical_geometry.polygon import SphericalPolygon

from geography import load_geojson, latlon2xyz

from geometry import (
    rotation_matrix_from_euler,
    rotation_matrix_from_src_dest_vecs,
    rotate_by_axis_angle,
)

from polyhedra import (
    icosahedron,
    truncated_icosahedron,
    icosahedron_face_transform,
    truncated_icosahedron_face_transform,
)

from dymaxion import DymaxionProjection

from svg import write_svg


# TODO: select other islands (like hawaii) from world-110m
# TODO: major lakes and seas
# TODO: use a proper 3d plotting library
#       alternatively: https://stackoverflow.com/questions/41699494/how-to-obscure-a-line-behind-a-surface-plot-in-matplotlib


def get_rotation(cfg):
    poly = cfg['projection']['polyhedron']
    name = cfg['projection']['rotation']

    # TODO: define rotations based on two earth reference points; north pole and ???
    if name in ['', 'identity', 'default']:
        return np.eye(3)
    if poly == 'icosahedron' and name == 'poles-at-vertices':
        return rotation_matrix_from_euler(y=np.pi*0.175, z=np.pi*0.0)
    if poly == 'icosahedron' and name == 'australia-face':
        return rotation_matrix_from_euler(x=np.pi*-0.03)
    if poly == 'truncated-icosahedron' and name == 'poles-at-pentagons':
        return rotation_matrix_from_euler(x=np.pi*0.17620819117478337)
    # if poly == 'icosahedron' and name == 'poles-at-faces':
    # if poly == 'icosahedron' and name == 'minimal-land-disruption':

    raise NotImplementedError(poly + "/" + name)

@pm
def main():
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    else:
        config_file = "configs/truncated-icosahedron-simple.json"

    with open(config_file, "r") as f:
        cfg = json.load(f)

    R = cfg['projection']['circumradius_in']

    # load border data
    shapes2d = load_geojson(cfg['map_data_spec'][0])
    shapes3d = latlon2xyz(1, shapes2d)

    # fixes self-intersecting geometry in continents.geo.json
    # TODO check geojson filename before doing this
    deletes = {
        1: [1909],  # north america
        7: [22, 23],  # antarctica
    }
    for shape_id, shape in enumerate(shapes3d):
        if shape_id in deletes:
            for d in deletes[shape_id][::-1]:
                shapes3d[shape_id] = np.vstack((shape[:d, :], shape[(d+1):, :]))


    # define polyhedron and projection
    if cfg['projection']['polyhedron'] in ['icosahedron', '20', 'icosa']:
        polyhedron = 'icosahedron'
        pv, pe, pf = icosahedron(circumradius=1)
        ft = icosahedron_face_transform
    if cfg['projection']['polyhedron'] in ['truncated-icosahedron', '32', 'soccerball']:
        polyhedron = 'truncated-icosahedron'
        pv, pe, pf = truncated_icosahedron(circumradius=1)
        ft = truncated_icosahedron_face_transform

        # compute rotation angle for poles-at-pentagons
        #verts = pv[pf[6]]
        #center_of_edge = np.mean(pv[[0, 2],:], axis=0)
        #center_of_pent = np.mean(verts,axis=0)
        #theta = np.arccos(np.dot(center_of_edge, center_of_pent)/(np.linalg.norm(center_of_edge) * np.linalg.norm(center_of_pent)))


    Rot = get_rotation(cfg)
    pv = pv @ Rot
    dym = DymaxionProjection(pv, pe, pf)
    dym.set_projection(cfg['projection']['method'])

    cnc_layout, W, H = generate_cnc_layout(shapes3d, dym, ft, R)

    # db()

    # json.dump(cnc_layout_simple[1]['paths'][26].tolist(), open('border-sample-australia.json', 'w'))

    # for direct comparison:
    # dym.set_projection('predistort-90')
    # cnc_layout_predistort = generate_cnc_layout(shapes3d, dym, ft)
    # cnc_layout_predistort[1]['plot_kwargs']['color'] = 'g'

    if cfg.get('plot2d'):
        # 2d plots
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # ax = fig.add_subplot(111, projection='3d')
        # plot_map_latlon(ax, shapes2d)

        plot_layers(ax, cnc_layout)
        # db()
        # plot_layers(ax, cnc_layout_predistort)
        ax.invert_yaxis()  # match svg coordinate system
        ax.set_aspect('equal')

    if cfg.get('svg'):
        # 96 pts per inch ???
        layout_svg, _, _ = adjust_layers(cnc_layout, scale=20)
        write_svg(layout_svg, cfg['svg_filename'], W*20, H*20)

    if cfg.get('dxf'):
        write_dxf(cnc_layout, cfg['dxf_filename'])
        # TODO: write dxf https://pypi.org/project/ezdxf/0.6.2/

    if cfg.get('plot3d'):
        # for 3d plots
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # globe stuff
        # plot_globe_sphere(ax, shapes3d)
        plot_globe_polyhedron(ax, shapes3d, dym)
        #dym.set_projection('predistort-90')
        #plot_globe_polyhedron(ax, shapes3d, dym)
        plot_polyhedron(ax, pv, pe)
        plot_polyhedron_labels(ax, pv, pf)

        # auxiliary stuff
        plot_polar_axis(ax, 1)
        ax.view_init(elev=-20, azim=130)
        # plot_latlon_grid(ax, 1, d=30)
        plot_hidden_cube(ax, 1)
        plt.xlabel('x')
        plt.xlabel('y')

        fig.tight_layout()
        fig.subplots_adjust(left=-0.5, right=1.5)

    plt.show()


def adjust_layers(layers, origin=False, invert=False, scale=1):
    # apply some simple transformations to all paths/points in a set of layers
    # - re-zero to the origin
    # - invert y-axis
    # - scale

    # get range of values
    lo = [1000000, 1000000]    # xmin, ymin
    hi = [-1000000, -1000000]  # xmax, ymax
    for l in layers:
        if l['type'] == 'polyline':
            for p in l['paths']:
                lo = np.min(np.vstack((p, lo)), axis=0)
                hi = np.max(np.vstack((p, hi)), axis=0)

    W, H = scale * (hi - lo)

    for n, l in enumerate(layers):
        # print(l['desc'])
        if l['type'] == 'text':
            p = l['pts']

            # scale, shift to origin, invert
            p = scale * p
            if origin:
                p = p-scale*lo
            if invert:
                p[:,1] = H - p[:,1]

            layers[n]['pts'] = p
        if l['type'] == 'polyline':
            for m, p in enumerate(l['paths']):

                # scale, shift to origin, invert
                p = scale * p
                if origin:
                    p = p-scale*lo
                if invert:
                    p[:,1] = H - p[:,1]

                layers[n]['paths'][m] = p

    return layers, W, H


def kludge_projection(xyz, c):
    # project Nx3 shape `xyz` to the plane with normal `c`,
    # using a goofy azimuthal projection that is identity under some limit,
    # and asymptotically approaches pi/2 above the limit. this maps the
    # entire sphere to the hemisphere centered on `c`, which should make
    # computing intersections reasonable.
    x0, y0, y1 = np.pi * 0.14, np.pi * 0.14, np.pi*0.5  # limit for truncated icosahedron is tan(1/2.478) = pi*0.1350
    h, k, m = x0+y0-y1, y1, -(y1-y0)**2
    f = lambda x: (m + k*(x-h))/(x-h)
    # this is a piecewise, asymptotic function designed such that
    # f(x) = x for x < x0
    # f(x0) = y0
    # f'(x0) = 1
    # f(inf) -> y1

    proj = []
    for pt in xyz:
        if np.linalg.norm(pt) < 1e-15:
            # dumb glitch
            continue

        a = np.arccos(np.dot(pt, c)/(np.linalg.norm(pt)*np.linalg.norm(c)))

        # https://en.wikipedia.org/wiki/Slerp
        omega = np.arccos(np.dot(pt, c))
        #if omega == 0:
        #    db()
        t = 1
        if a >= x0:
            t = f(a)/a
        v = np.sin((1-t)*omega)/np.sin(omega) * c + np.sin(t*omega)/np.sin(omega) * pt

        axis = np.cross(c, pt)
        proj.append(v)

    return np.array(proj)



def generate_cnc_layout(shapes3d, dym, face_transform, R):
    # plot the polyhedron net in 2d, with the corresponding projected shapes,
    # clipped by the faces they belong to.

    # dict of lists, indicating which shapes are included on each face.
    # some shapes will be included on multiple faces.
    # {face_id: []shape_id}
    shapes_on_face = {n: [] for n in range(len(dym.faces))}
    for shape_id, shape in enumerate(shapes3d):
        pxyz, faces = dym.project(shape)
        for face_id in list(set(faces)):
            shapes_on_face[face_id].append(shape_id)

    # for each face:
    # - 3d-rotate it and its corresponding shapes onto XY plane,
    # - use the face_transform function to adjust the layout within the XY plane
    # - compute intersection of face with shapes, to clip them properly.
    print('---- compute layout')
    edge_paths = []
    label_locs = []
    label_texts = []
    border_paths_dict = defaultdict(dict)

    for face_id, shape_ids in shapes_on_face.items():
        print('\nface id = %d' % face_id)
        fn = dym.face_unit_normals[face_id]                        # face normal
        fv_open = dym.vertices[dym.faces[face_id]]                 # face vertices
        fv = np.vstack((fv_open, fv_open[0,:]))                    # face vertices (closed shape)
        Rot = rotation_matrix_from_src_dest_vecs(fn, [0, 0, 1])    # rotation matrix to bring face into xy plane
        fv2 = fv @ Rot.T                                           # 3d face vertices rotated to xy plane
        fx, fy, fr = face_transform(face_id, fv2)                  # get 2d transformation parameters
        fRot = np.array([[np.cos(fr), -np.sin(fr)], [np.sin(fr), np.cos(fr)]])  # rotation matrix to adjust orientation of face within xy plane
        fv2_oriented = fv2[:, 0:2] @ fRot                          # 2d face vertices rotated to proper poly-net orientation

        edge_paths.append(fv2_oriented + [fx, fy])                 # 2d face vertices transformed to proper poly-net position

        label_locs.append([fx, fy])
        label_texts.append('%s' % face_id)

        # for each shape:
        # - project it from its sphere-surface xyz 3d points to an intermediate projection
        #     this retains the shape as is, in the vicinity of the face, but condenses the rest of it,
        #     so that the polyhedron-face projection doesn't blow up
        # - project the intermediate projection onto the (single known face of the) polyhedron
        # - transform the fully-projected shape to the xy plane, in correct poly-net orientation
        # - compute intersection of shape and face

        for shape_id in shape_ids:
            print('  shape id = %d' % shape_id)

            # NOTE: the best way to clip the shape to the face is to compute their
            # intersection as two polygons.
            # 1. project entire shape onto the face, then compute planar intersection.
            #    doesn't work because some shapes span too much of the globe, so the
            #    projection explodes
            # 2. compute intersection of spherical polygons.
            #    haven't yet found a decent library to do this, not worth
            #    implementing myself unless necessary
            # 3. compute planar intersection within some proper azimuthal projection
            #
            # 3b. use some ad-hoc azimuthal projection that works well enough

            warped = kludge_projection(shapes3d[shape_id], fn)  # 3b it is!

            # TODO generalize this face projection call
            proj3d = dym.project_simple_archimedean_face(shapes3d[shape_id], face_id)
            proj2d = proj3d @ Rot.T
            proj2d_oriented = proj2d[:, 0:2] @ fRot

            proj3d_warped = dym.project_simple_archimedean_face(warped, face_id)
            proj2d_warped = proj3d_warped @ Rot.T
            proj2d_warped_oriented = proj2d_warped[:, 0:2] @ fRot


            poly_shape = Polygon(proj2d_oriented)
            poly_shape_warped = Polygon(proj2d_warped_oriented)
            poly_face = Polygon(fv2_oriented)

            geometry_error = False

            poly_proj2d_clipped = None
            try:
                poly_proj2d_clipped = poly_shape.intersection(poly_face)
            except Exception as exc:
                print('    invalid geometry')
                geometry_error = True

            poly_proj2d_warped_clipped = None
            try:
                poly_proj2d_warped_clipped = poly_shape_warped.intersection(poly_face)
            except Exception as exc:
                print('    invalid geometry (warped)')
                geometry_error = True

            paths = []
            if type(poly_proj2d_clipped) == Polygon:
                # convert single polygon to multipolygon
                poly_proj2d_clipped = MultiPolygon([poly_proj2d_clipped])

            if type(poly_proj2d_clipped) == MultiPolygon:
                for poly in poly_proj2d_clipped:
                    paths.append(np.vstack(poly.exterior.coords.xy).T + [fx, fy])

            border_paths_dict[(face_id, shape_id)]['unwarped'] = paths


            paths = []
            if type(poly_proj2d_warped_clipped) == Polygon:
                poly_proj2d_warped_clipped = MultiPolygon([poly_proj2d_warped_clipped])

            if type(poly_proj2d_warped_clipped) == MultiPolygon:
                for poly in poly_proj2d_warped_clipped:
                    paths.append(np.vstack(poly.exterior.coords.xy).T + [fx, fy])

            border_paths_dict[(face_id, shape_id)]['warped'] = paths

            # debug info
            t1, t2 = 'MultiPolygon', 'MultiPolygon'
            if poly_proj2d_clipped is None:
                t1 = 'None'
            if poly_proj2d_warped_clipped is None:
                t2 = 'None'

            print('    clipped shape types = %s   %s' % (t1, t2))


    # assemble output
    layers = [
        {
            'desc': 'face-edges',
            'paths': edge_paths,
            'action': 'cut',
            'svg_kwargs': {
                'stroke': 'black',
                'fill-opacity': 0,
                'stroke_width': 0.1,
            },
            'plot_kwargs': {'color': 'k'},
            'type': 'polyline',
        },
        {
            'desc': 'face-labels',
            'pts': np.array(label_locs),
            'labels': label_texts,
            'type': 'text',
            'plot_kwargs': {'color': 'r'},
        }
    ]

    for (fid, sid), data in border_paths_dict.items():
        color = colorizer(sid)
        layers.append({
            'desc': 'continent-borders-warped-%02d-%02d' % (fid, sid),
            'paths': data['warped'],
            'action': 'cut',
            'svg_kwargs': {
                'stroke': floatcolor2hex(color),
                'fill-opacity': 0,
                'stroke_width': 0.1,
            },
            'plot_kwargs': {'color': color, 'linewidth': 1,},
            'type': 'polyline',
        })
        layers.append({
            'desc': 'continent-borders-unwarped-%02d-%02d' % (fid, sid),
            'paths': data['unwarped'],
            'action': 'cut',
            'svg_kwargs': {
                'stroke': floatcolor2hex(color),
                'fill-opacity': 0,
                'stroke_width': 0.1,
            },
            'plot_kwargs': {'color': color, 'linewidth': 1,},
            'type': 'polyline',
        })


    adjusted, W, H = adjust_layers(layers, origin=True, invert=True, scale=R)

    adjusted.append({
            'desc': 'inch-scale',
            'paths': [np.array([[1.0, 0], [0, 0], [0, 1]])],
            'action': 'draw',
            'svg_kwargs': {
                'stroke': 'black',
                'fill-opacity': 0,
                'stroke_width': 0.1,
            },
            'plot_kwargs': {'color': 'k'},
            'type': 'polyline',
    })

    return adjusted, W, H


def floatcolor2hex(floatcolor):
    return '#%02x%02x%02x' % tuple(int(255*x) for x in floatcolor)


def plot_layers(ax, layers):
    # print('plotting layout layers')
    for l in layers:
        # print('  %s' % l['desc'])
        if l['type'] == 'polyline':
            for p in l['paths']:
                # print('    len = %d' % len(p))
                ax.plot(p[:,0], p[:,1], **l['plot_kwargs'])
        if l['type'] == 'text':
            for p, txt in zip(l['pts'], l['labels']):
                ax.text(p[0], p[1], txt, **l['plot_kwargs'])


colormap = {}
def colorizer_hue(idx=None):
    if idx not in colormap:
        hue = idx/8
        colormap[idx] = mcolors.hsv_to_rgb([hue, 1.0, 1.0])
        # print('new color: %s -> %s' % (idx, colormap[idx]))
    return colormap[idx]


def colorizer_random(idx=None):
    # maintains a global map of arbitrary indexes to random colors
    # for repeatable, distinguishable colors
    if idx not in colormap:
        hue = np.random.random()
        colormap[idx] = mcolors.hsv_to_rgb([hue, 1.0, 1.0])
        # print('new color: %s -> %s' % (idx, colormap[idx]))
    return colormap[idx]


colorizer = colorizer_hue


def plot_globe_polyhedron(ax, shapes3d, dym):
    # plot border shapes projected onto polyhedron
    S = 1.0
    for n, s in enumerate(shapes3d):
        color = colorizer(n)
        pxyz, best_faces = dym.project(s)
        ax.plot(S*pxyz[:,0], S*pxyz[:,1], S*pxyz[:,2], '-', markersize=1, color=color, linewidth=1)


R_ci = 1.258408572364819 # ratio of icosahedron circumradius/inradius
def plot_polyhedron(ax, pv, pe):
    # plot wireframe of polyhedron
    for e in pe:
        v0, v1 = pv[e[0]], pv[e[1]]
        ax.plot(*zip(v0, v1), 'k-', linewidth=1, alpha=1)
        # ax.plot(*zip(v0/R_ci, v1/R_ci), 'k-', linewidth=1, alpha=1)


def plot_polyhedron_labels(ax, pv, pf):
    # text labels of face indexes
    for n, f in enumerate(pf):
        fc = np.mean(pv[f], axis=0)
        ax.text(fc[0], fc[1], fc[2], n, color='r')


def plot_polar_axis(ax, R):
    # plot earth axis for visual reference
    ax.plot([0, 0], [0, 0], [-R, R], '.b-')


def plot_hidden_cube(ax, R):
    # plot bounding cube to force better aspect ratio
    ax.plot(
        [-R, -R, -R, -R, R, R, R, R],
        [-R, -R, R, R, -R, -R, R, R],
        [-R, R, -R, R, -R, R, -R, R],
        'w.',
        markersize=1,
    )


def plot_latlon_grid(ax, R, d=30):
    # plot lat/lon lines for visual reference
    t = np.linspace(0, 2*np.pi, 64)
    for ll in np.arange(-90, 91, d):
        ax.plot(*sphere2cart(R, t/d2r, ll), 'b--', linewidth=1, alpha=0.5)  # constant latitude
        ax.plot(*sphere2cart(R, ll, t/d2r), 'b--', linewidth=1, alpha=0.5)  # constant longitude


def plot_globe_sphere(ax, shapes3d):
    # plot border shapes projected on sphere
    for n, s in enumerate(shapes3d):
        color = colorizer(n)
        ax.plot(s[:,0]/R_ci, s[:,1]/R_ci, s[:,2]/R_ci, '-', color=color, linewidth=1)


def plot_map_latlon(ax, shapes):
    # plot shapes in 2d with no projection
    for shape in shapes:
        lon, lat = zip(*shape)
        lon = np.array(lon)
        lat = np.array(lat)
        ax.plot(lon, lat, 'k-', linewidth=1)


if __name__ == '__main__':
    main()
