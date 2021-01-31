import svgwrite

def write_svg(layers, fname, svg_kwargs={}, units=None):
    # svg_kwargs intended to be just {'width': '10in', 'height': '10in'} but can be used for whatever
    scale = 20
    dwg = svgwrite.Drawing(fname, profile='tiny', **svg_kwargs)
    if units:
        scale = scale_map[units]

    for l in layers:
        print(l['desc'])
        if l['type'] == 'polyline':
            for p in l['paths']:
                el = svgwrite.shapes.Polyline(points=scale*p, **l['svg_kwargs'])
                dwg.add(el)
                # print('polyline(%d)' % len(p))
    dwg.save()
