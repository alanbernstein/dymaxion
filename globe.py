import json

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from geometry import geo_polyarea, icosahedron


# TODO: russia mainland isn't contiguous
# TODO: combine into continents (metadata has this info)

borders_file = 'world-110m.geo.json'

d2r = np.pi / 180
r_earth_mi = 3958.8  # radius of actual earth
r_ball_in = 10  # planned radius of CNC ball
R = r_ball_in

verbosity = 2

# DS = downsample rate
# plot_mode, DS = '3d-sphere', 10
plot_mode, DS = '3d-icosahedron', 1
# plot_mode, DS = '2d-latlon', 30
PLOT_GEO_GRID = False
country_whitelist = []
country_whitelist = ['Canada', 'United States', 'Australia', 'China', 'Egypt', 'Brazil', 'Germany', 'Russia', 'Antarctica']


# @pm
def main():

    with open(borders_file, 'r') as f:
        geoj = json.load(f)

    fig = plt.figure()
    if plot_mode.startswith('3d'):
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.add_subplot(111)


    counts_total = {
        'countries': len(geoj['features']),
        'shapes': 0,
        'points': 0,
    }
    counts_filtered = {
        'shapes': 0,
        'points': 0,
    }


    for country in geoj['features']:
        name = country['properties']['name']
        continent = country['properties']['continent']

        if country_whitelist and name not in country_whitelist:
            continue

        pop_est = country['properties']['pop_est']
        typ = country['geometry']['type']
        # typ=Polygon -> list of lists of points
        # typ=MultiPolygon -> list of lists of lists of points
        shapes = country['geometry']['coordinates']
        if typ == 'Polygon':
            # force same shapes structure for all types
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

        for n, (shape, area) in enumerate(shapes_by_area):
            # filter logic
            if n > 0 and area < 100000:
                # for each country, keep the N largest shapes, and all shapes above an area threshold
                continue

            if verbosity > 2:
                print('  %s  %7.2f km^2  %s' % (n, area, len(shape[0])))

            # downsample and close-loop the path
            lon, lat = zip(*shape[0])
            lon_ds = np.array(lon[::DS])
            lat_ds = np.array(lat[::DS])
            lon_ds = np.hstack((lon_ds, lon_ds[0]))
            lat_ds = np.hstack((lat_ds, lat_ds[0]))
            counts_filtered['shapes'] += 1
            counts_filtered['points'] += len(lat_ds)

            # plot
            if plot_mode == '2d-latlon':
                ax.plot(lon_ds, lat_ds, 'k-', linewidth=1)

            if plot_mode == '3d-sphere':
                x = R * np.cos(lon_ds * d2r) * np.cos(lat_ds * d2r)
                y = R * np.sin(lon_ds * d2r) * np.cos(lat_ds * d2r)
                z = R * np.sin(lat_ds * d2r)

                ax.plot(x, y, z, 'k-', linewidth=1)

            if plot_mode == '3d-icosahedron':
                x = R * np.cos(lon_ds * d2r) * np.cos(lat_ds * d2r)
                y = R * np.sin(lon_ds * d2r) * np.cos(lat_ds * d2r)
                z = R * np.sin(lat_ds * d2r)
                ax.plot(x, y, z, 'k-', linewidth=1)

                if PLOT_GEO_GRID:
                    # plot lat/lon lines for visual reference
                    t = np.linspace(0, 2*np.pi, 64)
                    d = 30
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
                    # ax.plot([0, 0], [0, 0], [-R, R], '.b-')  # poles



                # plot polyhedron
                pv, pe = icosahedron(circumradius=R)
                # TODO: rotate icosahedron to achieve contiguous landmass dymaxion projection
                #       and/or rotate so that icosahedron vertices are at earth poles
                for e in pe:
                    v0, v1 = pv[e[0]], pv[e[1]]
                    ax.plot(*zip(v0, v1), 'r--', linewidth=1, alpha=0.5)

                # TODO project country shapes onto polyhedron
                # for each face of polyhedron:
                # - define frustum bounds check function
                # - compute equation of plane
                # (these steps should be applied to a generic polyhedron,
                # not computed analytically for the predefined vertices,
                # so we can rotate the polyhedron arbitrarily first)

                # for each point in the shape:
                # - use frustum bounds checks to identify corresponding face
                # - use equation of plane to project point onto plane

                # draw the projected points instead

    if verbosity > -1:
        print('total   %6d points / %4d shapes / %3d countries' % (counts_total['points'], counts_total['shapes'], counts_total['countries']))
        print('fitered %6d points / %4d shapes / %3d countries' % (counts_filtered['points'], counts_filtered['shapes'], counts_total['countries']))



    plt.show()



if __name__ == '__main__':
    main()
