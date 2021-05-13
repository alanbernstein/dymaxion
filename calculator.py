from math import sqrt, pi
from ipdb import set_trace as db



######## practical parameters
# these two parameters fully determine the inner and outer icosahedra
T1 = wood_thickness = 0.75  # the wood stock we have to work with
T2 = corner_thickness = 0.125  # adjustable, must be >0 for proper appearance

######## practical limits
# physical bounds, which should be respected (ideally)
maximum_sphere_diameter = 20.0  # lathe capacity
minimum_edge_thickness = 0.25  # for glue strength



######## polyhedron relationships
# immutable geometric relationships of the polyhedron.
# these definitions are for:
# circumradius = circumradius_per_edge * edge_length
# etc.
# for non-platonic polyhedra, inradius here refers to the
# /smallest/ inradius, because need to sand down to that
# radius to achieve a sphere

phi = (sqrt(5)+1)/2
icosahedron_circumradius_per_edge = 1/2 * sqrt(phi*sqrt(5))  # .951
icosahedron_inradius_per_edge = phi ** 2/(2*sqrt(3))  # .756
icosahedron_midradius_per_edge = phi/2  # 0.809


# https://mathworld.wolfram.com/TruncatedIcosahedron.html
soccerball_circumradius_per_edge = 1/2 * sqrt(1 + 9 * phi**2)  # 2.478
soccerball_inradius5_per_edge = 1/2 * sqrt(1/10 * (125 + 41*sqrt(5)))
soccerball_inradius6_per_edge = 1/2 * sqrt(3/2 * (7+3*sqrt(5)))
soccerball_inradius_per_edge = min(soccerball_inradius5_per_edge, soccerball_inradius6_per_edge)
soccerball_midradius_per_edge = 3/2 * phi  # 2.427


polyhedron = 'soccerball'


if polyhedron == 'icosahedron':
    A = circumradius_per_edge = icosahedron_circumradius_per_edge
    B = inradius_per_edge = icosahedron_inradius_per_edge
    C = midradius_per_edge = icosahedron_midradius_per_edge
elif polyhedron == 'soccerball':
    A = circumradius_per_edge = soccerball_circumradius_per_edge
    B = inradius_per_edge = soccerball_inradius_per_edge
    C = midradius_per_edge = soccerball_midradius_per_edge


"""
calculator notes:

Corner Sanded thickness matters
Wood thickness             = outer_inradius minus inner_inradius
sanded diameter            = outer_inradius
unsanded envelope diameter = outer_circumradius
Edge length                = some function of circumradius (Wikipedia)

Solve this system for edge length as function of wood thickness and one of the outer diameters
"""

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
inner_edge_length = (wood_thickness - corner_thickness) / (circumradius_per_edge - inradius_per_edge)
outer_edge_length = inner_edge_length + wood_thickness/inradius_per_edge


# other polyhedron measurements, derived from the solution
inner_inradius = inner_edge_length * inradius_per_edge
inner_midradius = inner_edge_length * midradius_per_edge
inner_circumradius = inner_edge_length * circumradius_per_edge
outer_inradius = outer_edge_length * inradius_per_edge
outer_midradius = outer_edge_length * midradius_per_edge
outer_circumradius = outer_edge_length * circumradius_per_edge



##### verify limits
thickness_at_edge = outer_midradius - inner_midradius
unsanded_diameter = 2 * outer_circumradius
sanded_diameter = 2 * outer_inradius


print('wood thickness: %5.3f    corner_thickness: %5.3f' % (wood_thickness, corner_thickness))
print(polyhedron)
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
print('unsanded_diameter:      %5.3f  (< %5.3f)' % (unsanded_diameter, maximum_sphere_diameter))
print('sanded_diameter:        %5.3f' % (sanded_diameter))


"""
example:

wood thickness: 0.750    corner_thickness: 0.125
--------
inner edge length:  3.200
inner inradius:     2.419
inner midradius:    2.589
inner circumradius: 3.044
outer edge length:  4.193
outer inradius:     3.169
outer midradius:    3.392
outer circumradius: 3.987
---
thickness at midradius: 0.803  (> 0.250)
unsanded_diameter:      7.975  (< 20.000)
sanded_diameter:        6.337

"""
