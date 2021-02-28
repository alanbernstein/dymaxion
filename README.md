
Icosahedron [Dymaxion](https://en.wikipedia.org/wiki/Dymaxion_map) projection in 3d:
![Icosahedron Dymaxion 3D](images/icosahedron-dymaxion-3d.png)

Unfolded into 2d for printing or CNC:
![Icosahedron Dymaxion 2D](images/icosahedron-dymaxion-2d.png)

Truncated icosahedron projection:
![Truncated Icosahedron Dymaxion 3D](images/truncated-icosahedron-dymaxion-3d.png)

Unfolded into 2d for printing or CNC:
![Truncated Icosahedron Dymaxion 2D](images/truncated-icosahedron-dymaxion-2d.png)


# code
- main.py: entry point, includes all plot functions,
- dymaxion.py: implementation of the dymaxion projection.
- polyhedra.py: definitions of several polyhedra: vertices (generated from symmetry when possible), edges (automatic), faces (automatic), 2D nets (ad hoc)
- geometry.py: minor general-purpose 3d geometry functions
- geography.py: functions for dealing with geojson data.
- configs/*.json: config files to record combinations of settings. pass name as only argument to main.py.


# world map

## data
https://geojson-maps.ash.ms/

- shapefiles - clunky, old fashioned way of storing geo data
- geojson - human readable (so larger), easier to deal with


continents data: https://gist.github.com/hrbrmstr/91ea5cc9474286c72838?short_path=f3fde31
