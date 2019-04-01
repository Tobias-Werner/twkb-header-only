[![Build Status](https://travis-ci.org/Tobias-Werner/twkb-header-only.svg?branch=dev)](https://travis-ci.org/Tobias-Werner/twkb-header-only)
# Encoding TWKB geometries
This library encodes and decodes geometries by using the Tiny Well-Known Binary (TWKB) format. Neither prior compilation nor installation of dependencies required. Just include the file TWKB.h.

### Examples

##### Encoding 2D point
~~~c++
GeomFactory factory;

LocationXY pos(7.625752, 53.942254);
char precisionXY = 6;

bytes_t twkb = factory.makePoint(pos, precisionXY); 
~~~

##### Encoding 4D multiline (with integrated bounding box)
~~~c++
GeomFactory factory;

LocationXYZT pos1(7.625752, 53.942254, 10.175, 164.0);
LocationXYZT pos2(7.615752, 53.932254, 10.231, 165.0);
LocationXYZT pos3(7.532752, 53.915354, 10.335, 166.0);
LocationXYZT pos4(7.60679, 53.93277, 10.155, 171.0);
LocationXYZT pos5(7.61170, 53.93686, 10.158, 172.0);
LocationXYZT pos6(7.61577, 53.93547, 10.651, 173.0);

char precisionXY = 6;
char precisionZ = 3;
char precisionT = 0;
bool bbox = true;

vector<LocationXYZT> line1({pos1, pos2, pos3});
vector<LocationXYZT> line2({pos4, pos5, pos6});
vector<vector<LocationXYZT>> lines({line1, line2});

bytes_t twkb = factory.makeMultiLine(lines, precisionXY, precisionZ, precisionT, bbox);
~~~