import json
from pprint import pprint as pp
from ipdb import set_trace as db

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors as mcolors
import numpy as np


from geography import filter_geojson
from geometry import (
    icosahedron,
    DymaxionProjection,
    rotation_matrix_from_euler,
    rotation_matrix_from_src_dest_vecs,
    icosahedron_circumradius_per_side,
)
from svg import write_svg


# TODO: select other islands (like hawaii) from world-110m
# TODO: major lakes and seas
# TODO: use a proper 3d plotting library
#       alternatively: https://stackoverflow.com/questions/41699494/how-to-obscure-a-line-behind-a-surface-plot-in-matplotlib

# data_spec = {'fname': 'world-110m.geo.json', pmap: {'name': 'name', 'continent': 'continent'}}
data_spec = {'fname': 'continents.geo.json', 'pmap': {'name': 'CONTINENT', 'continent': 'CONTINENT'}}

SVG, PLOT2D, PLOT3D = True, True, True

svg_filename = 'icosahedron.svg'

d2r = np.pi / 180
r_ball_in = 10  # planned radius of CNC ball
R = r_ball_in

# count_thresh, area_thresh, DS = -1, 10000, 10
count_thresh, area_thresh, DS = 0, 100000, 1

name_whitelist = []
# name_whitelist = ['Canada', 'United States', 'Australia', 'China', 'Egypt', 'Brazil', 'Germany', 'Russia', 'Antarctica']


# @pm
def main():
    # load border data
    with open(data_spec['fname'], 'r') as f:
        geoj = json.load(f)
    shapes2d = filter_geojson(geoj, data_spec['pmap'], name_whitelist, count_thresh, area_thresh, DS)
    shapes3d = latlon2xyz(shapes2d)


    # define polyhedron and projection
    pv, pe, pf = icosahedron(circumradius=R)
    # Rot = np.eye(3)  # default
    # TODO: define these rotations based on the location of the north pole, and some other lat/lon reference point
    # Rot = rotation_matrix_from_euler(y=np.pi*0.175, z=np.pi*0.0)  # align two icosahedron vertices with earth poles
    Rot = rotation_matrix_from_euler(x=np.pi*-0.03)  # align australia to be contained in a face
    # Rot = rotation_matrix_from_euler(???)  # TODO: aligned poles to centers of faces
    pv = pv @ Rot
    dym = DymaxionProjection(pv, pe, pf)
    dym.set_projection('simple')
    # dym.set_projection('predistort')


    cnc_layout = generate_cnc_layout(shapes3d, dym)
    if PLOT2D:
        # 2d plots
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # ax = fig.add_subplot(111, projection='3d')
        # plot_map_latlon(ax, shapes2d)

        plot_layers(ax, cnc_layout)
        ax.set_aspect('equal')

    if SVG:
        write_svg(cnc_layout, svg_filename)
        # TODO: write dxf https://pypi.org/project/ezdxf/0.6.2/


    if PLOT3D:
        # for 3d plots
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # globe stuff
        # plot_globe_sphere(ax, shapes3d)
        plot_globe_polyhedron(ax, shapes3d, dym)
        plot_polyhedron(ax, pv, pe)
        plot_polyhedron_labels(ax, pv, pf)

        # auxiliary stuff
        plot_polar_axis(ax)
        ax.view_init(elev=-20, azim=130)
        # plot_latlon_grid(ax)
        plot_hidden_cube(ax)
        plt.xlabel('x')
        plt.xlabel('y')

    plt.show()



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
    L = R/icosahedron_circumradius_per_side  # side length

    p3 = np.pi/3
    r3_6 = np.sqrt(3)/6

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


def generate_cnc_layout(shapes3d, dym):
    # plot the polyhedron net in 2d, with the corresponding projected shapes

    # split up shapes into segments, based on which face they're on
    face_segments_map = {n: [] for n in range(20)}
    # this is a dict of lists, with each list element representing a segment
    # of a border-shape, that is contained within a single face of the polyhedron.
    # the list element is a tuple (vertex_list, face_id)

    for n, s in enumerate(shapes3d):
        # TODO: these are no longer closed loops, fix this. this could normally be
        # solved by appending x[0] to x, but the faceted projection of shapes
        # with variable point spacing means that some shapes can get cut through
        # long line segments. need to compute true start and end points by
        # intersecting with the polyhedron edge.
        pxyz, faces = dym.project(s)
        starts, lens, values = rlencode(faces)
        for s, l, v in zip(starts, lens, values):
            face_segments_map[v].append((pxyz[s:(s+l),:], n))

    # for each face, rotate it and its corresponding segments onto XY plane,
    # then use the face_transform function to adjust the layout
    edge_paths = []
    border_paths = []
    label_locs = []
    label_texts = []
    for face_idx, segment_list in face_segments_map.items():
        fn = dym.face_unit_normals[face_idx]
        Rot = rotation_matrix_from_src_dest_vecs(fn, [0, 0, 1])

        fv = dym.vertices[dym.faces[face_idx]]
        fv2 = fv @ Rot.T

        fx, fy, fr = icosahedron_face_transform(face_idx, fv2)
        fRot = np.array([[np.cos(fr), -np.sin(fr)], [np.sin(fr), np.cos(fr)]])
        fv2_oriented = fv2[:, 0:2] @ fRot

        edge_paths.append(np.vstack((
            fx + fv2_oriented[[0, 1, 2, 0], 0],
            fy + fv2_oriented[[0, 1, 2, 0], 1],
        )).T)

        label_locs.append([fx, fy])
        label_texts.append('%s' % face_idx)

        for segment3d, shape_idx in segment_list:
            color = colorizer(shape_idx)
            segment2d = segment3d @ Rot.T
            segment2d_oriented = segment2d[:, 0:2] @ fRot
            border_paths.append(np.vstack((
                fx + segment2d_oriented[:,0],
                fy + segment2d_oriented[:,1],
            )).T)


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
            'desc': 'continent-borders',
            'paths': border_paths,
            'action': 'cut',
            'svg_kwargs': {
                'stroke': 'red',
                'fill-opacity': 0,
                'stroke_width': 0.1,
            },
            'plot_kwargs': {'color': 'r'},
            'type': 'polyline',
        },
        {
            'desc': 'face-labels',
            'pts': label_locs,
            'labels': label_texts,
            'type': 'text',
            'plot_kwargs': {'color': 'b'},
        }
    ]

    return layers


def plot_layers(ax, layers):
    for l in layers:
        # db()
        if l['type'] == 'polyline':
            for p in l['paths']:
                ax.plot(p[:,0], p[:,1], **l['plot_kwargs'])
        if l['type'] == 'text':
            for p, txt in zip(l['pts'], l['labels']):
                ax.text(p[0], p[1], txt, **l['plot_kwargs'])


def rlencode(x):
    # https://gist.github.com/nvictus/66627b580c13068589957d6ab0919e66
    where = np.flatnonzero
    x = np.asarray(x)
    n = len(x)

    starts = np.r_[0, where(~np.isclose(x[1:], x[:-1], equal_nan=True)) + 1]
    lengths = np.diff(np.r_[starts, n])
    values = x[starts]

    return starts, lengths, values


colormap = {}
def colorizer(idx=None):
    # maintains a global map of arbitrary indexes to random colors
    # for repeatable, distinguishable colors
    if idx not in colormap:
        hue = np.random.random()
        colormap[idx] = mcolors.hsv_to_rgb([hue, 1, .6])
    return colormap[idx]


def plot_globe_polyhedron(ax, shapes3d, dym):
    # plot border shapes projected onto polyhedron
    S = 1.0
    for n, s in enumerate(shapes3d):
        color = colorizer(n)
        pxyz, best_faces = dym.project(s)
        ax.plot(S*pxyz[:,0], S*pxyz[:,1], S*pxyz[:,2], '-', color=color, linewidth=1)


def plot_polyhedron(ax, pv, pe):
    # plot wireframe of polyhedron
    R_ci = 1.258408572364819 # ratio of icosahedron circumradius/inradius
    for e in pe:
        v0, v1 = pv[e[0]], pv[e[1]]
        ax.plot(*zip(v0, v1), 'k-', linewidth=1, alpha=1)
        # ax.plot(*zip(v0*R_ci, v1*R_ci), 'k-', linewidth=1, alpha=1)


def plot_polyhedron_labels(ax, pv, pf):
    # text labels of face indexes
    for n, f in enumerate(pf):
        fc = np.mean(pv[f], axis=0)
        ax.text(fc[0], fc[1], fc[2], n, color='r')


def plot_polar_axis(ax):
    # plot earth axis for visual reference
    ax.plot([0, 0], [0, 0], [-R, R], '.b-')


def latlon2xyz(shapes2d):
    # apply spherical->cartesian transform for all shapes in list
    shapes3d = []
    for shape in shapes2d:
        lon, lat = zip(*shape)
        xyz = sphere2cart(R, np.array(lon), np.array(lat))
        shapes3d.append(np.vstack(xyz).T)

    return shapes3d


def sphere2cart(R, lon, lat):
    # transform spherical polar coordinates to cartesian
    x = R * np.cos(lon * d2r) * np.cos(lat * d2r)
    y = R * np.sin(lon * d2r) * np.cos(lat * d2r)
    z = R * np.sin(lat * d2r)
    return x, y, z


def plot_hidden_cube(ax):
    # plot bounding cube to force better aspect ratio
    ax.plot(
        [-R, -R, -R, -R, R, R, R, R],
        [-R, -R, R, R, -R, -R, R, R],
        [-R, R, -R, R, -R, R, -R, R],
        'w.',
        markersize=1,
    )


def plot_latlon_grid(ax, d=30):
    # plot lat/lon lines for visual reference
    t = np.linspace(0, 2*np.pi, 64)
    for ll in np.arange(-90, 91, d):
        ax.plot(*sphere2cart(R, t/d2r, ll), 'b--', linewidth=1, alpha=0.5)  # constant latitude
        ax.plot(*sphere2cart(R, ll, t/d2r), 'b--', linewidth=1, alpha=0.5)  # constant longitude


def plot_globe_sphere(ax, shapes3d):
    # plot border shapes projected on sphere
    for s in shapes3d:
        ax.plot(s[:,0], s[:,1], s[:,2], 'k-', linewidth=1)


def plot_map_latlon(ax, shapes):
    # plot shapes in 2d with no projection
    for shape in shapes:
        lon, lat = zip(*shape)
        lon = np.array(lon)
        lat = np.array(lat)
        ax.plot(lon, lat, 'k-', linewidth=1)


if __name__ == '__main__':
    main()
