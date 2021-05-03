import numpy as np
import svgwrite

def write_svg(layers, fname, W, H):
    dwg = svgwrite.Drawing(fname, profile='tiny', height='%din' % H, width='%din' % W)

    for l in layers:
        # print(l['desc'])
        if l['type'] == 'polyline':
            for p in l['paths']:
                rounded = np.round(p, decimals=4)
                el = svgwrite.shapes.Polyline(rounded, **l['svg_kwargs'])
                dwg.add(el)
                # print('polyline(%d)' % len(p))
    dwg.save()
