# topology_cut

This package includes all the necessary libs to build a block identificaiton program using the topology method proposed by Lin et al. (1987). Some ideas of this package comes from the work of Lu (2002). Most of these codes are finished during my visit to UC Berkeley.

The foundmental arithmetic of these libs is supported by CGAL (www.cgal.org). So you need to download and config CGAL properly if you want to compile this package. The version I am using is CGAL 4.6, but I believe the versions that below 4.6 still work.

ps: You can refers to the demo cases provided by CGAL offical release if you have trouble setting up CGAL for this package.

Any comments are welcomed, ~~especially~~ including the comments about the code styles (I know some codes may look ugly, so please help me to improve it if you are interested). Please feel free to send me pull requests. I am happy to see this package being helpful to you.

## Reference
* Lin D, Fairhurst C, Starfield A M. Geometrical identification of three-dimensional rock block systems using topological techniques. International Journal of Rock Mechanics & Mining Sciences & Geomechanics Abstracts, 1987, 24(6):331-338.
* Lu J. Systematic identification of polyhedral rock blocks with arbitrary joints and faults. Computers & Geotechnics, 2002, 29(1):49-72.

---
# TODO
* Delete unnecessary trace lines on the surface of the block
* Optimize the computation of cross sections of each fracture
* Translate the comments in the code to English
