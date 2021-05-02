import numpy as np


def rotate_by_axis_angle(v, k, t):
    # https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    return v * np.cos(t) + np.cross(k, v) * np.sin(t) + k * np.dot(k, v) * (1 - np.cos(t))


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
