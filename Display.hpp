//
//  Display.hpp
//  QuadrilateralFEM
//
//  Created by Champike Attanayake on 12/25/24.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "QuadMesh.hpp"

#ifndef Display_hpp
#define Display_hpp

void displayMatrix(const std::vector <std::vector <int>> mtx);
void displayMatrix(const std::vector <std::vector <double>> mtx);
void displayMatrix(const std::vector <int> vec);
void displayMatrix(const std::vector <double> vec);
void displayMesh(QuadMesh mesh, int nX, int nY, int deg, int param);

#endif /* Display_hpp */
