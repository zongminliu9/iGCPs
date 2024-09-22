#ifndef SETTINGS_H
#define SETTINGS_H

#include <iomanip>
#include "framework.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>


using namespace std;

// Problem configuration settings
class Settings {
public:
    Settings();
    void LoadVariables(string filename, vector<double>& variables);
    void InitializeCondition(int dof, vector<double>& initialVelocity, vector<double>& initialPressure, vector<CoordPoint3D>& points, const vector<array<double, 3>>& velocityNodes, const int dimension);
    void LoadMesh(string filename, vector<CoordPoint3D>& points, vector<CubeElement>& mesh);
    void LoadVelocityField(string filename, int numPoints, vector<array<double, 3>>& velocity);
    void AssignProcessor(string filename, int& numElements, vector<vector<int>>& elementProcessorMap);
};

#endif
