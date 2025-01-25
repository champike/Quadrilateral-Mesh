# Quadrilateral-Mesh
Generates a structured quadrilateral mesh of arbitrary order within any rectangular domain

This code generates a quadrilateral mesh of arbitrary order (linear, quadratic, cubic, etc.) within a rectangular domain. The mesh parameters are defined by the rectangle’s dimensions, horizontal and vertical partition sizes, and the desired element order.

Indexing in the code commences from the bottom left and increments counterclockwise.
## Methods
**Calculates a matrix X and Y coordinate of each degree of freedom.**
'''
globalDOFs()
'''
**Matrix of indices of degrees of freedom of each element.**
‘’’
ptrLocalDOFs()
‘’’
**Matrix of indices of vertices of each element.**
‘’’
ptrEleNodes()
‘’’
**Indices of elements on the boundary.**
‘’’
boundaryElements()
‘’’
**Indices of all nodes on the boundary.**
‘’’
boundaryNodes()
‘’’
**Tuple containing four vectors. Each vector holds boundary nodal indices on the bottom, top, left, and right edges, respectively.**
‘’’
boundarySpecificNodes()
‘’’
**Indices of edges of each element.**
‘’’
ElementEdges()



The mesh can be displayed using 
'''
displayMesh( . . . )
'''
function. 

## Dependencies
Qmesh utilizes some C++11 features and compiles successfully with both clang and gcc. To visualize the mesh, the gnuplot library (https://sourceforge.net/projects/gnuplot/) is required.

