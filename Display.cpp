//
//  Display.cpp
//  QuadrilateralFEM
//
//  Created by Champike Attanayake on 12/25/24.
//
#include <iostream> // only for testing
#include "Display.hpp"

/**
 * @brief Displays a matrix with integers.
 *
 * This function displays the contents of a 2D vector of integers
 * in a visually readable format.
 *
 * @param mtx The 2D vector of integers to be displayed.
 *
 */
void displayMatrix(const std::vector <std::vector <int>> mtx){
    // Iterate through the rows
    for (const auto &row : mtx) {
        // Iterate through the elements in each row
        for (auto val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}

/**
 * @brief Displays a matrix with double.
 *
 * This function displays the contents of a 2D vector of integers
 * in a visually readable format.
 *
 * @param mtx The 2D vector of double to be displayed.
 *
 */
void displayMatrix(const std::vector <std::vector <double>> mtx){
    // Iterate through the rows
    for (const auto &row : mtx) {
        // Iterate through the elements in each row
        for (auto val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}

/**
 * @brief Displays a vector with integers.
 *
 * This function displays the contents of a  vector of integers
 * in a visually readable format.
 *
 * @param vec The vector of integers to be displayed.
 *
 */
void displayMatrix(const std::vector <int> vec){
    // Iterate through the vect
    for (auto val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

/**
 * @brief Displays a vector with double.
 *
 * This function displays the contents of a  vector of integers
 * in a visually readable format.
 *
 * @param vec The vector of integers to be displayed.
 *
 */
void displayMatrix(const std::vector <double> vec){
    // Iterate through the vect
    for (auto val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}



/**
 * @brief Displays mesh information.
 *
 * This function displays key information about the mesh,
 * including the number of partitions, polynomial degree, and optionally
 * node indices depending on the `param` value.
 *
 * @param mesh The mesh data.
 * @param nX Number of partitions in the X-direction.
 * @param nY Number of partitions in the Y-direction.
 * @param deg Degree of the polynomial basis functions.
 * @param param Controls the display of node indices:
 *              - 0: Display the mesh without any indices.
 *              - 1: Display indices at the corners of each element.
 *              - 2: Display indices of all nodes in the mesh.
 */
void displayMesh(QuadMesh mesh, int nX, int nY, int deg, int param){
    const std::vector <std::vector <int>> eleNodes = mesh.ptrEleNodes();
    std::vector <std::vector <double>> nodalCoords = mesh.globalDOFs();
    int nEle = nX*nY;
    
    
    // Create an output file for the mesh
    std::ofstream meshfile("/Users/attanac/Downloads/mesh.txt");

    for (int i = 0; i < nEle; i++){
        for (int j = 0; j < 4; j++){
            meshfile << nodalCoords[eleNodes[i][j]][0] << " " << nodalCoords[eleNodes[i][j]][1] << std::endl;
        }
        // Put coordinates of the first node to the end. So it completes the quadrlateral
        meshfile << nodalCoords[eleNodes[i][0]][0] << " " << nodalCoords[eleNodes[i][0]][1] << std::endl;
        // Add an empty line to separate elements
        meshfile << std::endl;
    }
    meshfile.close();
    
    // Create a Gnuplot script
    
    std::vector<std::string> Gnucommands = {
        "unset key",
        "unset border", "unset xtics", "unset ytics",
        "set title \"Quadrilateral Mesh\"",
        "set xlabel \"X\"",
        "set ylabel \"Y\""
    };
    
    // display nodal indices at the corners of the elements
    if (param == 1){
        std::vector<std::string> GnuText((nX+1)*(nY+1));
        int nodes{0};
        int count = 0;
        // loop over each row
        for (int i = 0; i <=  nY; i++){
            for (int j = 0; j <= nX; j++){
                nodes= j*deg + i*deg*(deg*nX+1);
                // insert the label
                GnuText[count] = "set label '" + std::to_string(nodes) + "' at " + std::to_string(nodalCoords[nodes][0]) + "," + std::to_string(nodalCoords[nodes][1]-0.01);
                count++;
            }
        }
        // combine vector with labels
        Gnucommands.insert(Gnucommands.end(), GnuText.begin(), GnuText.end());
        Gnucommands.push_back("plot '/Users/attanac/Downloads/mesh.txt' pt 1 lt 5 lc 7 w lp");
    }
    // - - - - - -
    else if (param == 2){
        std::vector<std::string> GnuText(nodalCoords.size());
        std::ofstream nodefile("/Users/attanac/Downloads/node.txt");

        for (int i = 0; i < nodalCoords.size(); i++ ){
            GnuText[i] = "set label '" + std::to_string(i) + "' at " + std::to_string(nodalCoords[i][0]) + "," + std::to_string(nodalCoords[i][1]-0.01);
            nodefile << nodalCoords[i][0]<< " " << nodalCoords[i][1]<< std::endl;
        }
        // combine vector with labels
        Gnucommands.insert(Gnucommands.end(), GnuText.begin(), GnuText.end());
        nodefile.close();
        Gnucommands.push_back("plot '/Users/attanac/Downloads/mesh.txt' pt 1 lt 5 lc 7 w lp, '/Users/attanac/Downloads/node.txt' pt 28  lc 5 w p");
    }
    // - - - - - -
    else {
        Gnucommands.push_back("plot '/Users/attanac/Downloads/mesh.txt' pt 1 lt 5 lc 7 w lp");
    }
    


    FILE* gnupipe = popen("gnuplot -persistent", "w");
    if (gnupipe == nullptr) {
        std::cerr << "Error: Could not open pipe to gnuplot." << std::endl;
      //  return 1;
    }
    
    for (const std::string& command : Gnucommands) {
        fprintf(gnupipe, "%s\n", command.c_str());
    }
    
    pclose(gnupipe);
}
