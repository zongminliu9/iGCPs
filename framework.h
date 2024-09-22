#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <petsc.h>
#include <petscksp.h>

#include "petscsys.h"
#include "petscmat.h"

using namespace std;

/////////////////////////////////
class CoordPoint3D {
public:
    double coords[3];
    int marker; // ID for inlet and outlet
    CoordPoint3D();
};

class CubeElement {
public:
    int degree;
    int order;
    int numBasisFuncs;
    int elementType; // 0 for interior, 1 for boundary (for visualization)
    int isBezier; // 0 for spline, 1 for Bezier element
    int marker;

    vector<int> connectivity;
    vector<array<double, 3>> points; // Temporary points
    vector<array<double, 64>> controlMatrix;
    vector<int> boundaryConditionIDs;
    double velocity[3];

    CubeElement(int precision = 3);
    void BezierPoly(double u, vector<double>& shapeFuncU, vector<double>& gradShapeFuncU) const;
    void ComputeBasis(double u, double v, double w, vector<double>& shapeFuncT, vector<array<double, 3>>& gradShapeFuncT) const;
    void ParametricToPhysical(double u, double v, double w, double physPt[3]) const;
};

// Mesh format conversion
void ConvertRawToVtk(string filename);
void ConvertRawnToVtk(string filename);

#endif
