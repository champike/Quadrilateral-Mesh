//
//  QuadMesh.hpp
//  QuadrilateralFEM
//
//  Created by Champike A on 8/16/24.
//
#include <vector>
#include <tuple>
#include <iostream>
#include <random>
#ifndef QuadMesh_hpp
#define QuadMesh_hpp


class QuadMesh {
private:
    double leftX{0.0};
    double rightX{1.0};
    double bottomY{0.0};
    double topY{1.0};
    int nPartsX{5};
    int nPartsY{5};
    int order{1};
    double deviationX;
    double deviationY;
    int nElements;
    int deformType;
    
public:
    QuadMesh(int *Pparts, double *Pdims, int ord, double *divt, int flag);
    ~QuadMesh();
    // Creates a 2D vector of coordinates for all degrees of freedom in the quad mesh.
    std::vector <std::vector <double>> globalDOFs();
    
    // Creates a 2D vector containing each element's local degree of freedom indices.
    std::vector <std::vector <int>> ptrLocalDOFs();
    
    // Creates a 2D vector containing node indices of vertices for each quadrilateral element.
    std::vector <std::vector <int>> ptrEleNodes();
    
    // Creates a vector containing element indices for each quadrilateral element on the edges.
    std::vector <int> boundaryElements();
    
    // Creates a vector containing nodal indices for each quadrilateral element on the edges.
    std::vector<int> boundaryNodes();
    
    // Create four vectors containing nodal indices for each quadrilateral element on each  edge.
    std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>> boundarySpecificNodes();
    
    // Creates a 2D vector containing node indices for edges of  each quadrilateral element.
    std::vector<std::vector<int>> ElementEdges();
    
    
    std::vector<double> Translation(double x1, double eta, std::vector<std::vector <double>>& eleCrnrNd);

};
#endif /* QuadMesh_hpp */
