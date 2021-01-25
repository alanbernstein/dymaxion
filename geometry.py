from math import pi, cos, radians
import numpy as np

earth_radius = 6371009 # in meters

def geo_polyarea(lon, lat):
    # returns the area (in m^2) of a polygon defined in longitude,latitude,
    # via simple equal-area projection onto earth-sized sphere
    lat_dist = pi * earth_radius / 180.0

    y = [la * lat_dist for la in lat]
    x = [lo * lat_dist * cos(radians(la)) for la, lo in zip(lat, lon)]

    a = sum([x[i] * (y[i+1] - y[i-1]) for i in range(-1, len(x)-1)])
    return abs(a)/2.0



p = (np.sqrt(5) + 1)/2  # golden ratio

def icosahedron(circumradius=1):
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
    return vertices, edges


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
