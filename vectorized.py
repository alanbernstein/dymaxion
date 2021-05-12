import numpy as np
import svgwrite
"""
Vectorized is a python library for CAD and geometric design.
The primary functionality is:
- a set of classes representing vector objects, providing methods for both
  - plotting via python-based plotting libraries (such as matplotlib)
  - drawing into vectorized file formats (such as SVG, via svgwrite).

  The point is, I like to create designs for CNC projects by writing python code.
  I design using matplotlib, and I need an SVG file to send to a CNC machine.
  I want to generate both the plot, and the SVG file, from the same source vector objects.
  This library accomplishes that, with a concise API.

- a set of classes that provide "control point" functionality, simplifying
  geometric design tasks.


Other helper tools are also included:
- matplotlib helper tools
- plotly helper tools (?)
- geometric design helper tools
  - ControlGroup
- simple Color class

"""

# TODO: might want to think about incorporating SDFs as well.
# or maybe that's a separate class that just generates paths?
# or maybe this has another output format, which is go/sdfx code?


class Vectorized(object):
    """
    A Vectorized drawing is structured like this:
    drawing = []Vectorized
        children = []Vectorized

    Each leaf node in the tree can represent:
    - a list of paths
    - a list of text objects

    A path represents a single shape (assumed closed, which may not be important)

    This base class is nestable, which allows it to represent arbitrary
    grouping hierarchies. These can be SVG groups, logical layers, or
    any other arbitrary component.

    This design should enable writing to SVG as groups, instead of
    a flat list of paths. This is not yet implemented.

    This design also allows for calling bounding_box() and transform()
    on the entire hierarchy of objects, which is useful for the adjust()
    method, as well as general transformation actions.
    """
    def __init__(self, children, name=None):
        self.children = children or []
        self.name = name or ''

    def plot(self, ax):
        for c in self.children:
            c.plot(ax)

    def write_svg_file(self, filename):
        minxy, maxxy = self.bounding_box()
        w, h = maxxy - minxy  # TODO: this assumes the bounding box has one corner at the origin
        dwg = svgwrite.Drawing(filename, profile='tiny', height='%fin' % h, width='%fin' % w)
        self.write_svg(dwg)
        dwg.save()

    def write_svg(self, drawing):
        for c in self.children:
            c.write_svg(drawing)

    def write_dxf_file(self, filename):
        # TODO: write dxf https://pypi.org/project/ezdxf/0.6.2/
        raise NotImplementedError

    def transform(self, dx, dy, mat):
        for c in self.children:
            c.transform(dx, dy, mat)

    def bounding_box(self):
        minxy = [np.inf, np.inf]
        maxxy = [-np.inf, -np.inf]
        for n, c in enumerate(self.children):
            minxy_c, maxxy_c = c.bounding_box()
            minxy = np.nanmin(np.vstack((minxy_c, minxy)), axis=0)
            maxxy = np.nanmax(np.vstack((maxxy_c, maxxy)), axis=0)

        return minxy, maxxy

    def adjust(self, scale=1, invert=True, origin=True):
        # utility for adjusting a drawing to scale, and fit within the viewport
        # of an SVG definition.
        # scale: multiply every point in every path by this
        # invert: if true, invert the y dimension, and add h to the y-offset, because SVG places origin at top-left
        # origin: if true, translate entire drawing to a bounding box with one corner at the origin
        minxy, maxxy = self.bounding_box()

        left, bottom = minxy
        right, top = maxxy
        w, h = maxxy - minxy

        dx, dy = 0, 0
        mat = scale * np.eye(2)

        if origin:
            dx, dy = -left, -bottom

        if invert:
            dy -= h
            mat = mat @ [[1, 0], [0, -1]]

        # import ipdb; ipdb.set_trace()
        self.transform(dx, dy, mat)

    def dimensions(self):
        # TODO: avoid calling bounding_box twice when caller does adjust() and then dimensions()
        minxy, maxxy = self.bounding_box()
        w, h = maxxy - minxy
        return w, h


class TextGroup(Vectorized):
    def __init__(self, path, texts, name=None, color=None, draw_svg=False):
        self.paths = [path]  # just a list to keep consistent with other Layer types
        self.texts = texts or ['%d' for n in range(len(points))]
        self.name = name or ''
        if color is None:
            color = 'k'
        self.color = Color(color)
        self.draw_svg = draw_svg

        self.plot_kwargs = {
            'color': self.color.rgb(),
        }

    def plot(self, ax, **plot_kwargs):
        kw = {k: v for k, v in self.plot_kwargs.items()}
        if plot_kwargs is not None:
            kw.update(plot_kwargs)
        for p, t in zip(self.paths[0], self.texts):
            ax.text(p[0], p[1], t, **kw)

    def write_svg(self, drawing, **svg_kwargs):
        if not self.draw_svg:
            return
        kw = {k: v for k, v in self.svg_kwargs.items()}
        if svg_kwargs is not None:
            kw.update(svg_kwargs)
        raise NotImplementedError

    def transform(self, dx, dy, mat):
        for n in range(len(self.paths)):
            self.paths[n] += [dx, dy]
            self.paths[n] = self.paths[n] @ mat

    def bounding_box(self):
        if not self.draw_svg:
            return [np.inf, np.inf], [-np.inf, -np.inf]
        raise NotImplementedError


class PolylineGroup(Vectorized):
    def __init__(self, paths, name=None, color=None, action=None):
        # TODO: support unpack_lineformat
        # TODO: normalize paths as an Nx2 numpy array
        # TODO: accept several path formats, including:
        #       Nx1 complex array,
        #       pair of lists
        #       list of pairs
        self.paths = paths  # TODO: does this need to be plural?
        self.name = name or ''
        if color is None:
            color = 'k'
        self.color = Color(color)
        self.action = action or 'cut'

        self.plot_kwargs = {
            'color': self.color.rgb(),
            'linewidth': 1,
        }

        # default laser-appropriate settings
        self.svg_kwargs = {
            'stroke': self.color.to_hex(),
            'fill-opacity': 0,
            'stroke_width': 0.1,
            # epilog lasers require "hairline-width" strokes, which is not
            # actually a well-defined concept in SVG, but 0.1 is thin enough.
        }

    def bounding_box(self):
        minxy = [np.inf, np.inf]
        maxxy = [-np.inf, -np.inf]
        for p in self.paths:
            minxy = np.nanmin(np.vstack((p, minxy)), axis=0)
            maxxy = np.nanmax(np.vstack((p, maxxy)), axis=0)

        return minxy, maxxy

    def transform(self, dx, dy, mat):
        # this is intended for scaling and moving a drawing to keep it
        # in a tight viewport for reasonable SVG rendering.
        # for this use, applying translation before matrix is more sensible
        for n in range(len(self.paths)):
            self.paths[n] += [dx, dy]
            self.paths[n] = self.paths[n] @ mat

    def plot(self, ax, **plot_kwargs):
        # accepts at matplotlib.pyplot axis
        kw = {k: v for k, v in self.plot_kwargs.items()}
        if plot_kwargs is not None:
            kw.update(plot_kwargs)
        if 'label' not in kw and self.name != '':
            kw['label'] = self.name
        for p in self.paths:
            ax.plot(p[:,0], p[:,1], **kw)

    def write_svg(self, drawing, **svg_kwargs):
        # accepts an svgwrite.Drawing
        # TODO: write into a group
        kw = {k: v for k, v in self.svg_kwargs.items()}
        if svg_kwargs is not None:
            kw.update(svg_kwargs)

        for path in self.paths:
            rounded = np.round(path, decimals=4)
            el = svgwrite.shapes.Polyline(rounded, **kw)
            drawing.add(el)

        drawing.save()


class CircleGroup(Vectorized):
    def __init__(self, paths, radii, name=None, color=None):
        pass

    def bounding_box(self):
        pass

    def transform(self, dx, dy, mat):
        pass

    def plot(self, ax, **plot_kwargs):
        pass

    def write_svg(self, drawing, **svg_kwargs):
        pass


class ControlShape(object):
    # represents a set of control points
    # they can be referenced, and interpolated between
    # for easy representation of positions on a drawing
    pass


class ControlPolyline(ControlShape):
    # represents an arbitrary polyline
    def __init__(self, path):
        # path should be an Nx2 np.array
        # TODO make np.array if not already
        self.path = path

    def __getitem__(self, key):
        return self.path[key, :]

    def lerp(self, i0, i1, proportion):
        # returns a point Linearly intERPolated between the points `i0` and `i1`,
        # based on `proportion`, which should be in [0, 1], corresponding to [i0, i1]
        return self.path[i0] + proportion*(self.path[i1]-self.path[i0])

    def dist(self, i0, i1, dist):
        # returns a point at a distance `dist` along the line from `i0` to `i1`
        d = self.path[i1]-self.path[i0]
        d /= np.linalg.norm(d)
        return self.path[i0] + dist * d


class ControlCircle(ControlShape):
    # represents a circle
    # interpolation should work based on a supplied angle parameter
    def __init__(self, center, r):
        self.center = center
        self.r = r

        self.path = self.arc(0, 360)

    def __call__(self, degrees=None, radians=None):
        # instead of indexing a specific point, use the angle as the index
        radians = radians or degrees * np.pi/180 or 0
        return np.array([
            self.center[0] + self.r*np.cos(radians),
            self.center[1] + self.r*np.sin(radians),
        ])

    def arc(self, deg0, deg1, R=None):
        N = 128  # TODO choose this based on arc size
        rad0, rad1 = deg0*np.pi/180, deg1*np.pi/180
        radius = R or self.r
        t = np.linspace(rad0, rad1, N)
        return np.vstack((self.center[0] + radius*np.cos(t), self.center[1] + radius*np.sin(t))).T


class ControlEllipse(ControlShape):
    # represents an arbitrary ellipse
    # interpolation should work based on a supplied angle parameter
    def __init__(self, center, r1, r2=None, theta=None):
        self.center = center
        self.r1 = r1
        self.r2 = r2 or r1
        self.theta = theta or 0

    def __call__(self, radians=None, degrees=None):
        # instead of indexing a specific point, use the angle as the index
        radians = radians or degrees * np.pi/180
        raise NotImplementedError


class ControlConic(ControlShape):
    pass


class ControlSpline(ControlShape):
    pass



lineformat_parts = {
    'color': 'bgrcmykw',
    'linestyle': [':', '--', '-.', '-'],
    'marker': '.,ov^<>1234sp*hH+xDd|_',
}


def unpack_plot_kwargs(kwargs):
    # accept a dict of pyplot.plt kwargs
    # return a dict of pyplot.plt kwargs
    # break out "lineformat" value into the 3 values that are valid kwargs
    # the value of "lineformat" is a valid positional arg to plot(), but there doesn't seem to be
    # a corresponding keyword for using it as a kwarg. so instead, have to pass kwargs dict to this
    # function and then pass that to plot()

    # input is a dict with plot_kwargs, plus a 'lineformat' key, output is a plot_kwargs dict

    # the order of both the keys and the value elements is important here, to ensure substring checking works properly
    if 'lineformat' in kwargs:
        fmt = kwargs['lineformat']
        for k, v in lineformat_parts.items():
            for e in v:
                if e in fmt:
                    kwargs[k] = e
                    fmt = fmt.replace(e, '')
        del(kwargs['lineformat'])

        if 'marker' in kwargs and 'linestyle' not in kwargs:
            # if marker, but no linestyle, then don't draw any lines
            kwargs['linestyle'] = ''

    return kwargs


def unpack_lineformat(fmt):
    # input is a lineformat string, output is a plot_kwargs dict
    kwargs = {}
    for k, v in lineformat_parts.items():
        for e in v:
            if e in fmt:
                kwargs[k] = e
                fmt = fmt.replace(e, '')

    return kwargs


def is_linfeformat_string(s):
    return False

def test_unpack_plot_kwargs():
    cases = [
        ('r', {'color': 'r'}),
        ('rb', {'color': 'r'}),  # actually, this should error, but don't really care about that
        ('ro-', {'color': 'r', 'linestyle': '-', 'markerstyle': 'o'}),
        ('k-.p', {'color': 'k', 'linestyle': '-.', 'markerstyle': 'p'}),
        ('b--.', {'color': 'b', 'linestyle': '--', 'markerstyle': '.'}),
    ]
    for inp, exp in cases:
        outp = unpack_plot_kwargs({'lineformat': inp})
        if outp != exp:
            print(inp, exp, outp)


class Color(object):
    """
    Color is a simple class for concisely defining colors in a variety of formats,
    and then generating the equivalent representation in other formats.

    The canonical format is RGB, float in [0, 1]. Integers in [0, 255] would be closer
    to the hardware, but floats have more precision, which... might be useful?

    fully-specified:
    # TODO: implement
    Color(r=float, g=float, b=float)   _from_rgb_float
    Color(r=int, g=int, b=int)         _from_rgb_int
    Color(h=float, s=float, v=float)   _from_hsv
    Color(hex=#HHHHHH)                 _from_hex
    Color(hex=#HHH)                    _from_hex_short
    Color(name=string)                 _from_name
    Color(name=string)                 _from_char

    shortcuts:
    Color(float, float, float)    _from_rgb_float
    Color([float, float, float])  _from_rgb_float
    Color(int, int, int)          _from_rgb_int
    Color([int, int, int])        _from_rgb_int
    Color(#string)                _from_hex
    Color(string)                 _from_name
    """

    name_map = {
        'k': [0, 0, 0],
        'r': [1, 0, 0],
        'g': [0, 1, 0],
        'b': [0, 0, 1],
        'c': [0, 1, 1],
        'm': [1, 0, 0],
        'y': [1, 1, 0],
        'w': [1, 1, 1],
    }

    def __init__(self, *args, **kwargs):
        if len(args) == 3 and type(args[0]) is float:
            self._from_rgb_float(*args)
        elif len(args) == 1 and type(args[0]) is list and type(args[0][0]) is float:
            self._from_rgb_float(*args[0])
        elif len(args) == 1 and type(args[0]) is np.ndarray and args[0].dtype == float:
            self._from_rgb_float(*args[0])
        elif len(args) == 3 and type(args[0]) is int:
            self._from_rgb_int(*args)
        elif len(args) == 1 and type(args[0]) is list and type(args[0][0]) is int:
            self._from_rgb_int(*args[0])
        elif len(args) == 1 and type(args[0]) is np.ndarray and args[0].dtype == int:
            self._from_rgb_int(*args[0])

        elif len(args) == 1 and type(args[0]) is str and type(args[0][0]) == '#':
            self._from_hex(*args)
        elif len(args) == 1 and type(args[0]) is str:
            self._from_name(*args)
        else:
            print('Color initialization failed')
            import ipdb; ipdb.set_trace()

    def __repr__(self):
        return '<Color>'
        # return '<Color(r=%d, g=%d, b=%d)>' % (self.r, self.g, self.b)

    def _from_rgb_float(self, r, g, b):
        self.r, self.g, self.b = r, g, b

    def _from_rgb_int(self, r, g, b):
        self.r, self.g, self.b = r/255, g/255, b/255

    def _from_hex(self, hex_str):
        raise NotImplementedError

    def _from_name(self, name):
        if name not in self.name_map:
            import ipdb; ipdb.set_trace()
            raise ValueError
        self.r, self.g, self.b = self.name_map[name]

    def to_hex(self):
        return '#%02x%02x%02x' % (int(255*self.r), int(255*self.g), int(255*self.b))

    def rgb(self):
        # shorthand
        return self.to_rgb_float()

    def to_rgb_float(self):
        return [self.r, self.g, self.b]
