import json

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors as mcolors
import numpy as np

from geometry import (
    geo_polyarea,
    icosahedron,
    DymaxionProjection,
)

# TODO: select other islands (like hawaii) from world-110m
# TODO: major lakes and seas

# data_spec = {'fname': 'world-110m.geo.json', pmap: {'name': 'name', 'continent': 'continent'}}
data_spec = {'fname': 'continents.geo.json', 'pmap': {'name': 'CONTINENT', 'continent': 'CONTINENT'}}


d2r = np.pi / 180
r_earth_mi = 3958.8  # radius of actual earth
r_ball_in = 10  # planned radius of CNC ball
R = r_ball_in

verbosity = 3

# count_thresh, area_thresh, DS = -1, 10000, 10
count_thresh, area_thresh, DS = 0, 100000, 1

name_whitelist = []
# name_whitelist = ['Canada', 'United States', 'Australia', 'China', 'Egypt', 'Brazil', 'Germany', 'Russia', 'Antarctica']


# @pm
def main():

    with open(data_spec['fname'], 'r') as f:
        geoj = json.load(f)

    shapes = filter_country_geojson(geoj, name_whitelist, count_thresh, area_thresh, DS)

    fig = plt.figure()

    # ax = fig.add_subplot(111)  # for 2d plots
    ax = fig.add_subplot(111, projection='3d')  # for 3d plots


    # plot_2d_latlon(ax, shapes)
    # plot_3d_sphere(ax, shapes)
    plot_3d_icosahedron(ax, shapes)

    # plot_geo_grid(ax)

    plt.show()


def random_color():
    hue = np.random.random()
    return mcolors.hsv_to_rgb([hue, 1, .6])
    # return np.random.random((1, 3))


def plot_3d_icosahedron(ax, shapes):
    pv, pe, pf = icosahedron(circumradius=R)
    dym = DymaxionProjection(pv, pe, pf)

    # plot polyhedron
    R_ci = 1.258408572364819 # ratio of icosahedron circumradius/inradius
    # TODO: rotate icosahedron to achieve contiguous landmass dymaxion projection
    #       and/or rotate so that icosahedron vertices are at earth poles

    for shape in shapes:
        lon, lat = zip(*shape)
        lon = np.array(lon)
        lat = np.array(lat)
        x = R * np.cos(lon * d2r) * np.cos(lat * d2r)
        y = R * np.sin(lon * d2r) * np.cos(lat * d2r)
        z = R * np.sin(lat * d2r)
        color = random_color()
        # ax.plot(x, y, z, '-', color=color, linewidth=1)

        pxyz, faces = dym.project_cartesian(np.vstack((x, y, z)).T)
        ax.plot(1.05*pxyz[:,0], 1.05*pxyz[:,1], 1.05*pxyz[:,2], '-', color=color, linewidth=1)


    for e in pe:
        v0, v1 = pv[e[0]], pv[e[1]]
        ax.plot(*zip(v0, v1), 'k-', linewidth=1, alpha=1)
        # ax.plot(*zip(v0*R_ci, v1*R_ci), 'k-', linewidth=1, alpha=1)


def plot_geo_grid(ax, d=30):
    # plot lat/lon lines for visual reference
    t = np.linspace(0, 2*np.pi, 64)
    for ll in np.arange(-90, 91, d):
        # lines of constant latitude
        x = R*np.cos(t)*np.cos(ll*d2r)
        y = R*np.sin(t)*np.cos(ll*d2r)
        z = R*np.sin(ll*d2r) * np.ones(t.shape)
        ax.plot(x, y, z, 'b--', linewidth=1, alpha=0.5)

        # lines of constant longitude
        x = R*np.cos(ll*d2r)*np.cos(t)
        y = R*np.sin(ll*d2r)*np.cos(t)
        z = R*np.sin(t) * np.ones(t.shape)
        ax.plot(x, y, z, 'b--', linewidth=1, alpha=0.5)

    # ax.plot([0, 0], [0, 0], [-R, R], '.b-')  # polar axis


def plot_3d_sphere(ax, shapes):
    for shape in shapes:
        lon, lat = zip(*shape)
        lon = np.array(lon)
        lat = np.array(lat)
        x = R * np.cos(lon * d2r) * np.cos(lat * d2r)
        y = R * np.sin(lon * d2r) * np.cos(lat * d2r)
        z = R * np.sin(lat * d2r)

        ax.plot(x, y, z, 'k-', linewidth=1)


def plot_2d_latlon(ax, shapes):
    for shape in shapes:
        lon, lat = zip(*shape)
        lon = np.array(lon)
        lat = np.array(lat)
        ax.plot(lon, lat, 'k-', linewidth=1)


def filter_country_geojson(geojson, name_whitelist=None, count_thresh=0, area_thresh=100000, DS=10):
    # returns a list of shapes (list of lists of points).
    #
    # geojson file contains a list of administrative areas (countries, but also antarctica, etc)
    # each country contains a list of shapes, one for each contiguous land mass that belongs to the admin
    #
    # high resolution file: 544740 points / 4001 shapes / 246 countries
    # low resolution file : 10565 points /  285 shapes / 176 countrie
    #
    # so some amount of filtering needs to be done

    counts_total = {
        'countries': len(geojson['features']),
        'shapes': 0,
        'points': 0,
    }
    counts_filtered = {
        'shapes': 0,
        'points': 0,
    }

    filtered_shapes = []
    for country in geojson['features']:
        # map inconsistent geojson property keys to a consistent format
        props = {}
        for k, v in data_spec['pmap'].items():
            props[k] = country['properties'][v]
        name = props['name']
        continent = props['continent']

        # apply whitelist filter
        if name_whitelist and name not in name_whitelist:
            continue

        # force same structure for all shape types
        typ = country['geometry']['type']
        # typ=Polygon -> list of lists of points
        # typ=MultiPolygon -> list of lists of lists of points
        shapes = country['geometry']['coordinates']
        if typ == 'Polygon':
            shapes = [shapes]

        num_shapes = len(shapes)
        points_per_shape = [len(p[0]) for p in shapes]
        counts_total['shapes'] += num_shapes
        counts_total['points'] += sum(points_per_shape)

        if verbosity > 1:
            print('%s %s %s[%d]' % (name.ljust(30), continent.ljust(15), typ, num_shapes))
        elif verbosity > 3:
            print('%s %s %s[%d] %s' % (name.ljust(30), continent.ljust(15), typ, num_shapes, points_per_shape))

        # compute area, sort by area
        areas_real_km2 = []  # area in km^2 on earth
        for shape in shapes:
            areas_real_km2.append(geo_polyarea(*zip(*shape[0])) / 1e6)

        shapes_by_area = [(s, a) for a, s in sorted(zip(areas_real_km2, shapes), reverse=True)]

        # add acceptable shapes to result list
        for n, (shape, area) in enumerate(shapes_by_area):
            # filter logic
            if n > count_thresh and area < area_thresh:
                # for each country, keep the N largest shapes, and all shapes above an area threshold
                continue

            if verbosity > 2:
                print('  %s  %7.2f km^2  %s' % (n, area, len(shape[0])))

            shape_ds = shape[0][::DS]  # downsample
            shape_ds.append(shape_ds[0]) # ensure closed loop
            filtered_shapes.append(shape_ds)
            counts_filtered['shapes'] += 1
            counts_filtered['points'] += len(shape_ds)

    if verbosity > -1:
        print('total   %6d points / %4d shapes / %3d countries' % (counts_total['points'], counts_total['shapes'], counts_total['countries']))
        print('fitered %6d points / %4d shapes / %3d countries' % (counts_filtered['points'], counts_filtered['shapes'], counts_total['countries']))

    return filtered_shapes



if __name__ == '__main__':
    main()
