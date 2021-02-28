import json
from math import pi, cos, radians
import numpy as np
verbosity = 3


def load_geojson(cfg):
    # returns a list of shapes (list of lists of points).
    #
    # geojson file contains a list of administrative areas (countries, but also antarctica, etc)
    # each country contains a list of shapes, one for each contiguous land mass that belongs to the admin
    #
    # high resolution file: 544740 points / 4001 shapes / 246 countries
    # low resolution file : 10565 points /  285 shapes / 176 countrie
    #
    # so some amount of filtering needs to be done, at least for plotting

    fname = cfg['fname']
    prop_map = cfg['pmap']
    name_whitelist = cfg['name_whitelist']
    count_thresh = cfg['count_threshold']
    area_thresh = cfg['area_threshold_km2']
    DS = cfg['downsample_rate']

    with open(fname, 'r') as f:
        geojson = json.load(f)

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
        for k, v in prop_map.items():
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
        areas = []  # area in km^2 on earth
        for shape in shapes:
            areas.append(sphere_polyarea(*zip(*shape[0])) / 1e6)

        shapes_by_area = [(s, a) for a, s in sorted(zip(areas, shapes), reverse=True)]

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


def sphere_polyarea(lon, lat):
    # returns the area of a polygon defined in longitude,latitude,
    # via simple equal-area projection onto a unit sphere
    lat_dist = pi / 180.0

    y = [la * lat_dist for la in lat]
    x = [lo * lat_dist * cos(radians(la)) for la, lo in zip(lat, lon)]

    a = sum([x[i] * (y[i+1] - y[i-1]) for i in range(-1, len(x)-1)])
    return abs(a)/2.0


def latlon2xyz(R, shapes2d):
    # apply spherical->cartesian transform for all shapes in list
    shapes3d = []
    for shape in shapes2d:
        lon, lat = zip(*shape)
        xyz = sphere2cart(R, np.array(lon), np.array(lat))
        shapes3d.append(np.vstack(xyz).T)

    return shapes3d


d2r = np.pi / 180


def sphere2cart(R, lon, lat):
    # transform spherical polar coordinates to cartesian
    x = R * np.cos(lon * d2r) * np.cos(lat * d2r)
    y = R * np.sin(lon * d2r) * np.cos(lat * d2r)
    z = R * np.sin(lat * d2r)
    return x, y, z
