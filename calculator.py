from math import sqrt, pi
from ipdb import set_trace as db



######## practical parameters
# these two parameters fully determine the inner and outer icosahedra
T1 = wood_thickness = 0.75  # the wood stock we have to work with
T2 = corner_thickness = 0.125  # adjustable, must be >0 for proper appearance

######## practical limits
# physical bounds, which should be respected (ideally)
maximum_sphere_radius = 10.0  # lathe capacity
minimum_edge_thickness = 0.25  # for glue strength



######## icosahedron relationships
# immutable geometric relationships of the polyhedron
phi = (sqrt(5)+1)/2
A = icosahedron_circumradius_per_edge = 1/2 * sqrt(phi*sqrt(5))  # .951
B = icosahedron_inradius_per_edge = phi ** 2/(2*sqrt(3))  # .756
C = icosahedron_midradius_per_edge = phi/2  # 0.809
# this means:
# circumradius = A * edge_length
# inradius = B * edge_length


"""
here's some math:

1.
inradius, measured at center of face, is where wood thickness is measured, so:
outer_inradius = inner_inradius + wood_thickness
B * Lo = B * Li + T1

2.
after sanding down to achieve specified corner_thickness,
also ensuring minimal sanding, we know:
outer_inradius = inner_circumradius + corner_thickness
B * Lo = A * Li + T2

so, solve this system of equations:
B * Lo = B * Li + T1
B * Lo = A * Li + T2


0 = B * Li + T1 - (A * Li + T2)
0 = (B-A) * Li + T1 - T2
0 = (B-A) * Li + T1 - T2


(T2-T1)/(B-A) = Li                 or   (T1-T2)/(A-B) = Li
Lo = Li + T1/B
"""

####### system solution results
# these equations are the solution
inner_edge_length = (wood_thickness - corner_thickness) / (icosahedron_circumradius_per_edge - icosahedron_inradius_per_edge)
outer_edge_length = inner_edge_length + wood_thickness/icosahedron_inradius_per_edge


# other icosahedron measurements, derived from the solution
inner_inradius = inner_edge_length * icosahedron_inradius_per_edge
inner_midradius = inner_edge_length * icosahedron_midradius_per_edge
inner_circumradius = inner_edge_length * icosahedron_circumradius_per_edge
outer_inradius = outer_edge_length * icosahedron_inradius_per_edge
outer_midradius = outer_edge_length * icosahedron_midradius_per_edge
outer_circumradius = outer_edge_length * icosahedron_circumradius_per_edge



##### verify limits
thickness_at_edge = outer_midradius - inner_midradius
unsanded_diameter = outer_circumradius
sanded_diameter = outer_inradius


print('wood thickness: %5.3f    corner_thickness: %5.3f' % (wood_thickness, corner_thickness))
print('--------')
print('inner edge length:  %5.3f' % inner_edge_length)
print('inner inradius:     %5.3f' % inner_inradius)
print('inner midradius:    %5.3f' % inner_midradius)
print('inner circumradius: %5.3f' % inner_circumradius)

print('outer edge length:  %5.3f' % outer_edge_length)
print('outer inradius:     %5.3f' % outer_inradius)
print('outer midradius:    %5.3f' % outer_midradius)
print('outer circumradius: %5.3f' % outer_circumradius)

print('---')

print('thickness at midradius: %5.3f  (> %5.3f)' % (thickness_at_edge, minimum_edge_thickness))
print('unsanded_diameter:      %5.3f  (< %5.3f)' % (2*unsanded_diameter, maximum_sphere_radius))
print('sanded_diameter:        %5.3f' % (2*sanded_diameter))
