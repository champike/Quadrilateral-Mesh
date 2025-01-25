# Quadrilateral-Mesh
Generates a structured quadrilateral mesh of arbitrary order within any rectangular domain


This code generates a quadrilateral mesh of arbitrary order (linear, quadratic, cubic, etc.) within a rectangular domain. The mesh parameters are defined by the rectangle’s dimensions, horizontal and vertical partition sizes, and the desired element order.

![plt1](https://github.com/user-attachments/assets/3ec347cc-bbad-48d8-8172-f61dbc8c7c29)

Indexing in the code commences from the bottom left and increments counterclockwise.
## Methods
* `globalDOFs()` Calculates a matrix X and Y coordinate of each degree of freedom.**

* `ptrLocalDOFs()` Matrix of indices of degrees of freedom of each element.**

* `ptrEleNodes()` Matrix of indices of vertices of each element.**

* `boundaryElements()` Indices of elements on the boundary.**

* `boundaryNodes()` Indices of all nodes on the boundary.**

* `boundarySpecificNodes()` Tuple containing four vectors. Each vector holds boundary nodal indices on the bottom, top, left, and right edges, respectively.**

* `ElementEdges()` Indices of edges of each element.**


The mesh can be displayed using 
`
displayMesh( . . . )
`
function. 

## Dependencies
Qmesh utilizes some C++11 features and compiles successfully with both clang and gcc. To visualize the mesh, the gnuplot library (https://sourceforge.net/projects/gnuplot/) is required.

