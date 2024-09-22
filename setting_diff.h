#ifndef SETTING_DIFF_H
#define SETTING_DIFF_H
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#endif
#ifdef __linux__ 
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "framework.h"
#endif
using namespace std;

class SettingDiff {
public:
    SettingDiff();
    

    void SetVariables(const string& paramFile, vector<double>& variables, string& outputPath);
    

    void SetInitialCondition(vector<Vertex3D>& points, vector<double>& variables, vector<double>& N0_ini, vector<double>& Nplus_ini, vector<double>& Nminus_ini);


    void ReadMesh(const string& meshFile, vector<Vertex3D>& points, vector<Element3D>& elements);
    

    void ReadVelocityField(const string& velocityFile, int pointCount, vector<array<double, 3>>& velocity);
    

    void AssignProcessor(const string& partitionFile, int& elementCount, vector<vector<int>>& elementProcess);


    void InitializeBoundaryConditions(const string& boundaryFile, vector<Vertex3D>& points, vector<double>& NplusBoundary, vector<double>& NminusBoundary);


    void ExportMesh(const string& outputFile, const vector<Vertex3D>& points, const vector<Element3D>& elements);


    void ProcessInputData(const string& filename, vector<Vertex3D>& points, vector<Element3D>& elements);
};

#endif
