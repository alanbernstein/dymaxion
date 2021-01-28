import json

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors as mcolors
import numpy as np

from geography import filter_geojson
from geometry import (
    icosahedron,
    DymaxionProjection,
    rotation_matrix_from_euler,
)

# TODO: select other islands (like hawaii) from world-110m
# TODO: major lakes and seas
# TODO: use a proper 3d plotting library
#       alternatively: https://stackoverflow.com/questions/41699494/how-to-obscure-a-line-behind-a-surface-plot-in-matplotlib

# data_spec = {'fname': 'world-110m.geo.json', pmap: {'name': 'name', 'continent': 'continent'}}
data_spec = {'fname': 'continents.geo.json', 'pmap': {'name': 'CONTINENT', 'continent': 'CONTINENT'}}


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
    Rot = np.eye(3)
    # Rot = rotation_matrix_from_euler(y=np.pi*0.175, z=np.pi*0.0)  # align two icosahedron vertices with poles
    Rot = rotation_matrix_from_euler(x=np.pi*-0.03)  # align australia to be contained in a face
    pv = pv @ Rot
    dym = DymaxionProjection(pv, pe, pf)
    dym.set_projection('simple')
    # dym.set_projection('predistort')


    fig = plt.figure()
    # ax = fig.add_subplot(111)  # for 2d plots
    ax = fig.add_subplot(111, projection='3d')  # for 3d plots

    # globe stuff
    # plot_map_latlon(ax, shapes2d)
    # plot_globe_sphere(ax, shapes3d)
    plot_globe_polyhedron(ax, shapes3d, dym)
    plot_polyhedron(ax, pv, pe)

    plot_cnc_layout(ax, shapes3d, dym)

    # auxiliary stuff
    plot_polar_axis(ax)
    ax.view_init(elev=-20, azim=130)
    # plot_latlon_grid(ax)
    plt.xlabel('x')
    plt.xlabel('y')

    plt.show()


def plot_cnc_layout(ax, shapes3d, dym):
    # plot the polyhedron net in 2d, with the corresponding projected shapes
    for s in shapes3d:
        color = random_color()
        # ax.plot(s[:,0], s[:,1], s[:,2], '-', color=color, linewidth=1)
        pxyz, faces = dym.project(s)
        ax.plot(1.05*pxyz[:,0], 1.05*pxyz[:,1], 1.05*pxyz[:,2], '-', color=color, linewidth=1)


def plot_globe_polyhedron(ax, shapes3d, dym):
    # plot border shapes projected onto polyhedron
    for s in shapes3d:
        color = random_color()
        # ax.plot(s[:,0], s[:,1], s[:,2], '-', color=color, linewidth=1)
        pxyz, faces = dym.project(s)
        ax.plot(1.05*pxyz[:,0], 1.05*pxyz[:,1], 1.05*pxyz[:,2], '-', color=color, linewidth=1)


def random_color():
    hue = np.random.random()
    return mcolors.hsv_to_rgb([hue, 1, .6])
    # return np.random.random((1, 3))


def plot_polyhedron(ax, pv, pe):
    # plot wireframe of polyhedron
    R_ci = 1.258408572364819 # ratio of icosahedron circumradius/inradius
    for e in pe:
        v0, v1 = pv[e[0]], pv[e[1]]
        ax.plot(*zip(v0, v1), 'k-', linewidth=1, alpha=1)
        # ax.plot(*zip(v0*R_ci, v1*R_ci), 'k-', linewidth=1, alpha=1)


def plot_polar_axis(ax):
    # plot earth axis for visual reference
    ax.plot([0, 0], [0, 0], [-R, R], '.b-')


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


def plot_map_latlon(ax, shapes):
    # plot shapes in 2d with no projection
    for shape in shapes:
        lon, lat = zip(*shape)
        lon = np.array(lon)
        lat = np.array(lat)
        ax.plot(lon, lat, 'k-', linewidth=1)


if __name__ == '__main__':
    main()
