//
//  main.cpp
//  QuadrilateralFEM
//
//  Created by Champike A on 8/16/24.
//

#include <iostream>
#include <vector>
#include "QuadMesh.hpp"
#include "Display.hpp"



int main(int argc, const char * argv[]) {
    
    int nXparts = 10; // Number of partitions on X axis
    int nYparts = 15; // Number of partitions on Y axis
    int order = 2; // degree linear, quadratic cubic, . . .
    double leftX {-1.0}, rightX{1.0}, bottomY{-1.0}, topY{1.0}; // rectangular domain boundary
    
    // deform rectanular element into quadrilateral element
    // If the variable "defFlag == 1" the  deformation of rectangular elements is done in a random scale
    // If the variable "defFlag == 0" the  deformation of rectangular elements is done in the scale
    // of alphaX horizontally and alphaY vertically
    int defFlag = 0;
    //0 < alphaX,alphaY < 1
    double alphaX = 0.5; //deform scale alphaX/nXpart horizontally
    double alphaY = 0.3;  //deform scale alphaY/nYpart vertically
   // int boundary[4] = {0,1,1,0};
    
    double meshParamsDims[] = {leftX, rightX, bottomY, topY};
    int meshParamsParts[] = {nXparts, nYparts};
    double alphas[] = {alphaX,alphaY};
    
    QuadMesh mesh(meshParamsParts,meshParamsDims,order,alphas,defFlag);
   
    // Visualize  vectors,  matrices (only in int & double), and the mesh defined in QadMesh
    //displayMatrix(mesh.globalDOFs());
    
    /*
     * The variable parameter Controls the display of node indices:
     *              - 0: Display the mesh without any indices.
     *              - 1: Display indices at the corners of each element.
     *              - 2: Display indices of all nodes in the mesh.
     */
    int parameter = 1;
    displayMesh(mesh, nXparts, nYparts, order, parameter);
    
   
    return 0;
}
