//
//  QuadMesh.cpp
//  QuadrilateralFEM
//
//  Created by Champike A on 8/16/24.
//

#include "QuadMesh.hpp"
#include "Display.hpp"
/**
 * @brief Constructs a QuadMesh object using the partition information and order given.
 *
 * This constructor initializes the QuadMesh object with the specified number of partitions in the x and y directions,
 * the domain boundaries, and the polynomial order of the basis functions.
 *
 * This constructor initializes the member variables of the QuadMesh object:
 * nPartsX: Number of partitions in the x-direction.
 * nPartsY: Number of partitions in the y-direction.
 * leftX, rightX: Left and right endpoints of the x-axis.
 * bottomY, topY: Bottom and top endpoints of the y-axis.
 * order: Polynomial order of the basis functions.
 *
 * @param Pparts A pointer to an integer array containing the number of partitions in the x and y directions.
 * @param Pdims A pointer to a double array containing the left and right endpoints of the x-axis and the bottom and top endpoints of the y-axis.
 * @param ord The polynomial order of the basis functions.
 *
 * @throws std::invalid_argument If any of the input parameters are invalid.
 */
QuadMesh::QuadMesh(int *Pparts, double *Pdims, int ord, double *divt, int flag) {
    
    // Check for invalid input parameters
    if (Pparts == nullptr || Pdims == nullptr) {
        throw std::invalid_argument("Null pointer passed to QuadMesh constructor.");
    }
    
    if (Pparts[0] <= 0 || Pparts[1] <= 0) {
        throw std::invalid_argument("Number of partitions must be positive.");
    }
    
    if (Pdims[0] >= Pdims[1] || Pdims[2] >= Pdims[3]) {
        throw std::invalid_argument("Invalid domain dimensions.");
    }
    
    if (ord < 0) {
        throw std::invalid_argument("Polynomial order must be non-negative.");
    }
    
    if (divt[0] < 0 || divt[0] >= 1) {
        throw std::invalid_argument("Element diviation term must be between 0 and 1.");
    }
    
    if (divt[1] < 0 || divt[1] >= 1) {
        throw std::invalid_argument("Element diviation term must be between 0 and 1.");
    }
    
    if ( flag != 0 && flag != 1) {
      //  throw std::invalid_argument("Flag must be either 0 or 1.");
    }
    // Initialize member variables
    nPartsX = Pparts[0];
    nPartsY = Pparts[1];
    leftX = Pdims[0];
    rightX = Pdims[1];
    bottomY = Pdims[2];
    topY = Pdims[3];
    order = ord;
    deviationX = divt[0];
    deviationY = divt[1];
    deformType = flag;
    nElements = nPartsX * nPartsY; // Calculate the total number of elements
}


/**
 * @brief Creates a 2D vector of coordinates for all degrees of freedom in the quad mesh.
 *
 * This method calculates the X and Y coordinates of each node in the mesh,
 * considering the number of elements in each direction, the order of the elements,
 * and the domain boundaries.
 *
 * @return nodes.  A 2D vector of size (nNodesX * nNodesY) x 2, where each inner vector contains
 *         the X and Y coordinates of a node.
 */
std::vector <std::vector <double>> QuadMesh::globalDOFs(){
    
    // Calculate the total number of nodes in each direction
    int nNodesX = (order*nPartsX + 1);
    int nNodesY = (order*nPartsY + 1);
    int count{0};
    
    // Initialize a 2D vector to store nodal coordinates
    std::vector<std::vector<double>> nodes(nNodesX*nNodesY, std::vector<double>(2, 0.0));
    
    // Calculate coordinates  in X,Y direction
    double dx = (rightX-leftX)/(nNodesX-1.0);
    double dy = (topY-bottomY)/(nNodesY-1.0);
    
    for (int i = 0; i < nNodesY; i++ ){
        for (int j = 0; j < nNodesX; j++){
            nodes[count][0] = (double)leftX + j*dx;
            nodes[count][1] = (double)bottomY + i*dy;
            count++;
        }
    }
    
    // Use deviation parameter to deform rectangular elements to quadrlateral
    // first consider the four corners
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
    // deform on X coordinate
    count = -1;
    int ndIdx = 0;
    double devV = dx*deviationX;
    for (int iy = 0; iy <= nPartsY; iy++){
        count = -count;
        for (int ix = 1; ix < nPartsX; ix++){
            ndIdx = ix*order + iy * nNodesX * order;
            nodes[ndIdx][0] = nodes[ndIdx][0] + deformType*count*dx*dis(gen) + (1-deformType)*count*devV;
        }
    }
    
    // deform on Y coordinate
    
    count = -1;
    devV = dy*deviationY;
    for (int iy = 0; iy < nPartsY-1; iy++){
        for (int ix = 0; ix <= nPartsX; ix++){
            count = -count;
            ndIdx = ix*order + (iy+1) * nNodesX * order;
            nodes[ndIdx][1] = nodes[ndIdx][1] + deformType*count*dy*dis(gen) + (1-deformType)*count*devV;
        }
    }

    
    // - - ONLY if order > 1
    /* Second consider nodes except on the
     -  four corners,
     - right edge
     - top edge
     to avoid repeat calculations.
     */
    //Based on the matrix created in the ptrLocalDOFs() create a vector with the required nodal indices.
    int reqLocalNd [order*order];
    count = 0;
    // zeroth element (node on the bottom left corner) included only for simplicity of the calculation
    for (int iy = 0; iy < order; iy++){
        for (int ix = 0; ix < order; ix++){
            reqLocalNd[count] = ix + iy*(order+1);
            count++;
        }
    }
    
    std::vector <std::vector <int>> cornerNodes = QuadMesh::ptrEleNodes();
    // Create nodes on the reference element
    double ref1D = 2.0/order; // distance between nodes on reference element
    double refNd[(order+1)*(order+1)][2]; // reference nodes
    count = 0;
    for (int iy = 0; iy < order+1; iy++){
        for (int ix = 0; ix < order+1; ix++){
            refNd[count][0] = -1.0 + double(ix)*ref1D;
            refNd[count][1] = -1.0 + double(iy)*ref1D;
            count++;
        }
    }
    std::vector<std::vector<int>> dofsL = ptrLocalDOFs();
    std::vector<double> transXY = {0.0,0.0};
    std::vector <std::vector <double>> CornerCoords(4, std::vector<double>(2,0.0));
    //loop over each element
    for (int ele = 0; ele < nElements; ele++){
        CornerCoords = {
            nodes[cornerNodes[ele][0]],
            nodes[cornerNodes[ele][1]],
            nodes[cornerNodes[ele][2]],
            nodes[cornerNodes[ele][3]]
        };
       
        for (int i = 1; i < order*order; i++){ // skip the first element at the botton left corner
            nodes[dofsL[ele][reqLocalNd[i]]] = Translation(refNd[reqLocalNd[i]][0],refNd[reqLocalNd[i]][1],CornerCoords);
        }
    }
    
    // Transform nodes on the right edge of the elements on the right edge
    int reqLocalNdEdge[order-1];
    for (int i = 0; i < order-1; i++){
         reqLocalNdEdge[i] = (order+1)*(i+2)-1;
    }
    
    int ele = 0;
    for (int i = 0; i < nPartsY; i++){
        ele = nPartsX*(i+1) - 1;
        CornerCoords = {
            nodes[cornerNodes[ele][0]],
            nodes[cornerNodes[ele][1]],
            nodes[cornerNodes[ele][2]],
            nodes[cornerNodes[ele][3]]
        };
        for (int i = 0; i < order-1; i++){ // skip the first element at the botton left corner
            nodes[dofsL[ele][reqLocalNdEdge[i]]] = Translation(refNd[reqLocalNdEdge[i]][0],refNd[reqLocalNdEdge[i]][1],CornerCoords);
        }
    }
    // Transform nodes on the top edge of the elements on the top edge
    for (int i = 0; i < order-1; i++){
         reqLocalNdEdge[i] = (order+1)*(order+1)-order+i;
    }
    
    for (int i = 0; i < nPartsX; i++){
        ele = nPartsX*nPartsY - nPartsX + i;
        CornerCoords = {
            nodes[cornerNodes[ele][0]],
            nodes[cornerNodes[ele][1]],
            nodes[cornerNodes[ele][2]],
            nodes[cornerNodes[ele][3]]
        };
        for (int i = 0; i < order-1; i++){ // skip the first element at the botton left corner
            nodes[dofsL[ele][reqLocalNdEdge[i]]] = Translation(refNd[reqLocalNdEdge[i]][0],refNd[reqLocalNdEdge[i]][1],CornerCoords);
        }
    }
    return nodes;
}


/**
 * @brief Creates a vector containing a transformed point's x- and y-coordinates.
 *
 * This method transforms a point from the given reference element to the specified quadrilateral.
 *
 * @param xi  x- coordinate of the point in the  reference element for the transformation.
 * @param eta  y- coordinate of the point in the  reference element for the transformation.
 *
 * @return xyCoords. A vector containing the x and y coordinates of the transformed point.
 */
std::vector<double> QuadMesh::Translation(double xi, double eta, std::vector<std::vector <double>>& vertices){
    std::vector<double> xyCoords = {0,0};
    // Shape functions
    std::vector<double> psi = {(1-xi)*(1-eta)/4,
        (1+xi)*(1-eta)/4,
        (1+xi)*(1+eta)/4,
        (1-xi)*(1+eta)/4};
    for (int i = 0; i < 4; i++){
        xyCoords[0] += psi[i]*vertices[i][0];
        xyCoords[1] += psi[i]*vertices[i][1];
    }
    return xyCoords;
}



/**
 * @brief Creates a 2D vector containing each element's local degree of freedom indices.

 * This method determines the indices of the local degrees of freedom for each element in the mesh.
 * The indices are arranged in a specific order, corresponding to the local node numbering within each element.

 * @return dofsL. A 2D vector of size (nElements) x (order+1)^2, where each inner vector contains the indices of the local degrees of freedom of an element.
 */
std::vector <std::vector <int>> QuadMesh::ptrLocalDOFs(){
    
    std::vector<std::vector<int>> dofsL(nPartsX*nPartsY, std::vector<int>((order+1)*(order+1), 0));
    
    int count{0};
    
    for (int iy = 0; iy < nPartsY; iy++){
        for (int ix = 0; ix < nPartsX; ix++){
            for(int j = 0; j <= order; j++){
                // calculate & insert nodes on the bottom raw of the element "count"
                dofsL[count][j] = iy*order*(order*nPartsX+1) + ix*order+j;
                for (int k = 1; k <= order; k++){
                    // calculate & insert nodes on the other raws of the element "count"
                    dofsL[count][j+k*(order+1)] = dofsL[count][j] + k*(nPartsX*order+1);
                }
            }
            count++;
        }
    }
    return dofsL;
}

/**
 * @brief Creates a 2D vector containing node indices for vertices of  each quadrilateral element.

 * This method determines the indices of the four corner nodes for each quadrilateral element in the mesh.
 * The indices are arranged counter-clockwise, starting from the bottom-left corner.

 * @return eleNodes. A 2D vector of size (nPartsX * nPartsY) x 4, where each inner vector contains the indices of the four corner nodes of an element.
 */

std::vector <std::vector <int>> QuadMesh::ptrEleNodes(){

    std::vector<std::vector<int>> eleNodes(nPartsX*nPartsY, std::vector<int>(4, 0));
    
    // Iterate over each element and calculate node indices
    int count{0};
    for (int iy = 0; iy < nPartsY; iy++){
        for (int ix = 0; ix < nPartsX; ix++){
            // Calculate the index of the bottom-left corner node
            eleNodes[count][0] = ix*order + iy*order*(order*nPartsX+1);
            // Calculate the indices of the other three corner nodes
            eleNodes[count][1] = eleNodes[count][0] + order;
            eleNodes[count][3] = eleNodes[count][0] + (order*nPartsX+1)*order;
            eleNodes[count][2] = eleNodes[count][3] + order;
            
            count++;
        }
    }
    return eleNodes;
}

/**
 * @brief Creates a vector containing element indices for each quadrilateral element on the edges.
 *
 * This method creates a vector with indices of elements which are
 * arranged counter-clockwise, on the bottom, left, right, and top edges respectively.
 *
 * @return bndyEles. A vector of size 2*nPartsX + 2*nPartsY - 4,  consisting of indices of elements on the bottom, left, right, and top edges respectively
 *
 */
std::vector<int> QuadMesh::boundaryElements() {
    
    // Create a vector to store the boundary element indices
    std::vector<int> bndyEles(2 * nPartsX + 2 * nPartsY - 4, 0);
    
    // Initialize a counter for the left edge elements
    int leftCounter = nPartsX + 2 * (nPartsY - 2);
    
    // Insert elements in the bottom and the top raw
    for (int i = 0; i < nPartsX; ++i) {
        bndyEles[i] = i;
        bndyEles[leftCounter + i] = nPartsX * (nPartsY - 1) + i;
    }
    
    // Insert the elements on the left and right boundary edges
    for (int i = 0; i < nPartsY - 2; ++i) {
        leftCounter = nPartsX + 2 * i;
        bndyEles[leftCounter] = (i + 1) * nPartsX; // Left boundary element
        bndyEles[leftCounter + 1] = bndyEles[leftCounter] + nPartsX - 1; // Right boundary element
    }
    
    return bndyEles;
}



/**
 * @brief Creates a vector containing nodal indices for each quadrilateral element on the edges.
 *
 * This method creates a vector with indices of nodes which are
 * arranged in the order  of bottom, left, right, and top edges respectively.
 *
 * @return A vector of size 2 * (order * nPartsX + 1) + 2 * (order * nPartsY - 1)4,
 * consisting of indices of elements on the bottom, left, right, and top edges respectively
 *
 */
std::vector<int> QuadMesh::boundaryNodes() {
    // Create a vector to store the boundary nodal indices
    std::vector<int> bndyNds(2 * (order * nPartsX + 1) + 2 * (order * nPartsY - 1), 0);
    
    // Initialize a counter for the nodes on the left edge
    int leftCounter = 2*order*(nPartsY - 1) + 2*(order-1) + order * nPartsX + 1;
    
    // Insert elements in the bottom and top  row respectively
    for (int i = 0; i < order * nPartsX + 1; ++i) {
        bndyNds[i] = i;
        bndyNds[leftCounter + i] = order*nPartsY * (order*nPartsX + 1) + i;
    }
    
    // Insert the nodes on the left and right boundary edges
    for (int i = 0; i < order*(nPartsY-1)+order-1; i++){
        leftCounter = order*nPartsX + 1 + 2*i;
        bndyNds[leftCounter] = (order*nPartsX + 1)*(i+1);
        bndyNds[leftCounter + 1] = bndyNds[leftCounter] + order*nPartsX;
        
    }
    return bndyNds;
}

/**
 * @brief Create four vectors containing nodal indices for each quadrilateral element on the specific  edges.
 *
 * This method creates four  vectors with indices of nodes which are
 * arranged in the order  of bottom, top, left, and right  edges respectively.
 *
 * @return A tuple with four vectors ),
 * consisting of indices of elements on the bottom, left, right, and top edges respectively
 *
 */
std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>> QuadMesh::boundarySpecificNodes(){
    std::vector<int> bndyNdsBottom(order * nPartsX + 1,0);
    std::vector<int> bndyNdsTop(order * nPartsX + 1,0);
    // Corner nodes counted into top and bottom raws
    std::vector<int> bndyNdsLeft(order * nPartsY - 1,0);
    std::vector<int> bndyNdsRight(order * nPartsY - 1,0);
    
    // Insert elements in the bottom and top  row respectively
    for (int i = 0; i < order * nPartsX + 1; ++i) {
        bndyNdsBottom[i] = i;
        bndyNdsTop[i] = order*nPartsY * (order*nPartsX + 1) + i;
    }
    
    // Insert the nodes on the left and right boundary edges
    for (int i = 0; i < order*(nPartsY-1)+order-1; i++){
        bndyNdsLeft[i] = (order*nPartsX + 1)*(i+1);
        bndyNdsRight[i] = bndyNdsLeft[i] + order*nPartsX;
    }
    return std::make_tuple(bndyNdsBottom,bndyNdsTop,bndyNdsLeft,bndyNdsRight);
}

/**
 * @brief Creates a 2D vector containing node indices for edges of  each quadrilateral element.
 
 * This method determines the indices of the four edges  for each quadrilateral element in the mesh.
 * The indices are arranged counter-clockwise, starting from the bottom.
 
 * @return A 2D vector of size (nPartsX * nPartsY) x 4, where each inner vector contains the indices of the four corner nodes of an element.
 */

std::vector<std::vector<int>> QuadMesh::ElementEdges(){
    std::vector<std::vector<int>> eleEdges(nPartsX*nPartsY, std::vector<int>(4, 0));
    int count = 0;
    for (int iy = 0; iy < nPartsY; iy++){
        for (int ix = 0; ix < nPartsX; ix++){
            // Calculate the index of the bottom-left corner node
            eleEdges[count][0] = ix + iy * (2*nPartsX+1);
            // Calculate the indices of the other three corner nodes
            eleEdges[count][3] = eleEdges[count][0] + nPartsX;
            eleEdges[count][2] = eleEdges[count][0] + 2*nPartsX + 1;
            eleEdges[count][1] = eleEdges[count][3] + 1;
            count++;
        }
    }
    return eleEdges;
}

// Distructor
QuadMesh::~QuadMesh()
{
};
