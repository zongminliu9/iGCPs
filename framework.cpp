#include "framework.h"

CubeElement::CubeElement(int precision) {
    degree = precision;
    order = precision + 1;
    numBasisFuncs = order * order * order;
    connectivity.resize(8);
    points.resize(8);
    for (int i = 0; i < 8; i++) {
        connectivity[i] = 0;
        points[i][0] = 0.;
        points[i][1] = 0.;
        points[i][2] = 0.;
    }
}

Point3D::Point3D() {
    coords[0] = 0.;
    coords[1] = 0.;
    coords[2] = 0.;
    marker = 0;
}


void CubeElement::ComputeBasis(double u, double v, double w, vector<double>& shapeFuncT, vector<array<double, 3>>& gradShapeFuncT) const {
    vector<double> shapeFuncU, shapeFuncV, shapeFuncW, gradShapeFuncU, gradShapeFuncV, gradShapeFuncW;
    BezierPoly(u, shapeFuncU, gradShapeFuncU);
    BezierPoly(v, shapeFuncV, gradShapeFuncV);
    BezierPoly(w, shapeFuncW, gradShapeFuncW);
    shapeFuncT.resize(numBasisFuncs);
    gradShapeFuncT.resize(numBasisFuncs);
    int i, j, k, loc = 0;
    for (k = 0; k < order; k++) {
        for (j = 0; j < order; j++) {
            for (i = 0; i < order; i++) {
                shapeFuncT[loc] = shapeFuncU[i] * shapeFuncV[j] * shapeFuncW[k];
                gradShapeFuncT[loc][0] = gradShapeFuncU[i] * shapeFuncV[j] * shapeFuncW[k];
                gradShapeFuncT[loc][1] = shapeFuncU[i] * gradShapeFuncV[j] * shapeFuncW[k];
                gradShapeFuncT[loc][2] = shapeFuncU[i] * shapeFuncV[j] * gradShapeFuncW[k];
                loc++;
            }
        }
    }
}


void CubeElement::BezierPoly(double paramU, vector<double>& shapeFuncU, vector<double>& gradShapeFuncU) const {
    if (degree == 3) {
        double shapeFunc0[4] = {(1. - paramU) * (1. - paramU) * (1. - paramU), 
                               3. * (1. - paramU) * (1. - paramU) * paramU, 
                               3. * (1. - paramU) * paramU * paramU, 
                               paramU * paramU * paramU};
        double gradShapeFunc0[4] = {-3. * (1. - paramU) * (1. - paramU), 
                                    3. - 12. * paramU + 9. * paramU * paramU, 
                                    3. * (2. - 3. * paramU) * paramU, 
                                    3. * paramU * paramU};
        shapeFuncU.resize(order);
        gradShapeFuncU.resize(order);
        for (int i = 0; i < order; i++) {
            shapeFuncU[i] = shapeFunc0[i];
            gradShapeFuncU[i] = gradShapeFunc0[i];
        }
    } else if (degree == 4) {
        double shapeFunc0[5] = {(1. - paramU) * (1. - paramU) * (1. - paramU) * (1. - paramU), 
                               4. * (1. - paramU) * (1. - paramU) * (1. - paramU) * paramU, 
                               6. * (1. - paramU) * (1. - paramU) * paramU * paramU, 
                               4. * (1. - paramU) * paramU * paramU * paramU, 
                               paramU * paramU * paramU * paramU};
        double gradShapeFunc0[5] = {-4. * (1. - paramU) * (1. - paramU) * (1. - paramU), 
                                    4. * (1. - paramU) * (1. - paramU) * (1. - 4. * paramU), 
                                    12. * paramU * (1. - 3. * paramU + 2. * paramU * paramU), 
                                    4. * (3. - 4. * paramU) * paramU * paramU, 
                                    4. * paramU * paramU * paramU};
        shapeFuncU.resize(order);
        gradShapeFuncU.resize(order);
        for (int i = 0; i < order; i++) {
            shapeFuncU[i] = shapeFunc0[i];
            gradShapeFuncU[i] = gradShapeFunc0[i];
        }
    }
}




void CubeElement::ParametricToPhysical(double u, double v, double w, double physPt[3]) const {
    vector<double> shapeFuncT;
    vector<array<double, 3>> gradShapeFuncT;
    ComputeBasis(u, v, w, shapeFuncT, gradShapeFuncT);
    physPt[0] = 0.;
    physPt[1] = 0.;
    physPt[2] = 0.;
    for (int i = 0; i < numBasisFuncs; i++) {
        physPt[0] += points[i][0] * shapeFuncT[i];
        physPt[1] += points[i][1] * shapeFuncT[i];
        physPt[2] += points[i][2] * shapeFuncT[i];
    }
}

void CubeElement::ParametricToPhysical(double u, double v, double w, double physPt[3]) const {
    vector<double> shapeFuncT;
    vector<array<double, 3>> gradShapeFuncT;
    ComputeBasis(u, v, w, shapeFuncT, gradShapeFuncT);
    physPt[0] = 0.;
    physPt[1] = 0.;
    physPt[2] = 0.;
    for (int i = 0; i < numBasisFuncs; i++) {
        physPt[0] += points[i][0] * shapeFuncT[i];
        physPt[1] += points[i][1] * shapeFuncT[i];
        physPt[2] += points[i][2] * shapeFuncT[i];
    }
}

void CubeElement::ComputeBasis(double u, double v, double w, vector<double>& shapeFuncT, vector<array<double, 3>>& gradShapeFuncT) const {
    vector<double> shapeFuncU, shapeFuncV, shapeFuncW, gradShapeFuncU, gradShapeFuncV, gradShapeFuncW;
    BezierPoly(u, shapeFuncU, gradShapeFuncU);
    BezierPoly(v, shapeFuncV, gradShapeFuncV);
    BezierPoly(w, shapeFuncW, gradShapeFuncW);
    shapeFuncT.resize(numBasisFuncs);
    gradShapeFuncT.resize(numBasisFuncs);
    int i, j, k, loc = 0;
    for (k = 0; k < order; k++) {
        for (j = 0; j < order; j++) {
            for (i = 0; i < order; i++) {
                shapeFuncT[loc] = shapeFuncU[i] * shapeFuncV[j] * shapeFuncW[k];
                gradShapeFuncT[loc][0] = gradShapeFuncU[i] * shapeFuncV[j] * shapeFuncW[k];
                gradShapeFuncT[loc][1] = shapeFuncU[i] * gradShapeFuncV[j] * shapeFuncW[k];
                gradShapeFuncT[loc][2] = shapeFuncU[i] * shapeFuncV[j] * gradShapeFuncW[k];
                loc++;
            }
        }
    }
}


void ConvertRawToVtk(string filename) {
    unsigned int numPoints, numElements;
    vector<array<double, 3>> points;
    vector<array<int, 8>> connectivity;
    double tempValue;
    string rawFile = filename + ".raw";
    ifstream inputFile;
    inputFile.open(rawFile);
    
    if (inputFile.is_open()) {
        inputFile >> numPoints >> numElements;
        points.resize(numPoints);
        connectivity.resize(numElements);
        for (unsigned int i = 0; i < numPoints; i++) {
            inputFile >> points[i][0] >> points[i][1] >> points[i][2] >> tempValue;
            inputFile >> tempValue;
        }
        for (unsigned int i = 0; i < numElements; i++) {
            for (int j = 0; j < 8; j++) inputFile >> connectivity[i][j];
            inputFile >> tempValue;
        }
        inputFile.close();
    } else {
        cerr << "Cannot open " << rawFile << "!\n";
    }
    
    string vtkFile = filename + ".vtk";
    ofstream outputFile;
    outputFile.open(vtkFile);
    
    if (outputFile.is_open()) {
        outputFile << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        outputFile << "POINTS " << points.size() << " float\n";
        for (unsigned int i = 0; i < points.size(); i++) {
            outputFile << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
        }
        outputFile << "\nCELLS " << connectivity.size() << " " << 9 * connectivity.size() << '\n';
        for (unsigned int i = 0; i < connectivity.size(); i++) {
            outputFile << "8 ";
            for (int j = 0; j < 8; j++) {
                outputFile << connectivity[i][j] << ' ';
            }
            outputFile << '\n';
        }
        outputFile << "\nCELL_TYPES " << connectivity.size() << '\n';
        for (unsigned int i = 0; i < connectivity.size(); i++) {
            outputFile << "12\n";
        }
        outputFile.close();
    } else {
        cerr << "Cannot open " << vtkFile << "!\n";
    }
}

void ConvertRawnToVtk(string filename) {
    unsigned int numPoints, numElements;
    vector<array<double, 3>> points;
    vector<array<int, 8>> connectivity;
    double tempValue;
    string rawnFile = filename + ".rawn";
    ifstream inputFile;
    inputFile.open(rawnFile);
    
    if (inputFile.is_open()) {
        inputFile >> numPoints >> numElements;
        points.resize(numPoints);
        connectivity.resize(numElements);
        for (unsigned int i = 0; i < numPoints; i++) {
            inputFile >> points[i][0] >> points[i][1] >> points[i][2] >> tempValue >> tempValue >> tempValue >> tempValue;
        }
        for (unsigned int i = 0; i < numElements; i++) {
            for (int j = 0; j < 8; j++) inputFile >> connectivity[i][j];
        }
        inputFile.close();
    } else {
        cerr << "Cannot open " << rawnFile << "!\n";
    }
    
    string vtkFile = filename + ".vtk";
    ofstream outputFile;
    outputFile.open(vtkFile);
    
    if (outputFile.is_open()) {
        outputFile << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        outputFile << "POINTS " << points.size() << " float\n";
        for (unsigned int i = 0; i < points.size(); i++) {
            outputFile << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
        }
        outputFile << "\nCELLS " << connectivity.size() << " " << 9 * connectivity.size() << '\n';
        for (unsigned int i = 0; i < connectivity.size(); i++) {
            outputFile << "8 ";
            for (int j = 0; j < 8; j++) {
                outputFile << connectivity[i][j] << ' ';
            }
            outputFile << '\n';
        }
        outputFile << "\nCELL_TYPES " << connectivity.size() << '\n';
        for (unsigned int i = 0; i < connectivity.size(); i++) {
            outputFile << "12\n";
        }
        outputFile.close();
    } else {
        cerr << "Cannot open " << vtkFile << "!\n";
    }
}
