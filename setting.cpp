#include "settings.h"

Settings::Settings() {}

void Settings::LoadVariables(string filename, vector<double>& variables) {
    variables.resize(12);
    string tempString;
    stringstream stream;
    ifstream file;
    file.open(filename);

    if (file.is_open()) {
        for (int i = 0; i < 12; i++) {
            file >> tempString >> variables[i];
        }
        file.close();
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}

void Settings::LoadMesh(string filename, vector<CoordPoint3D>& points, vector<CubeElement>& mesh) {
    string tempString;
    int numPoints, numElements, tempInt;
    ifstream file;
    file.open(filename);

    if (file.is_open()) {
        for (int i = 0; i < 4; i++) getline(file, tempString); // Skip lines
        file >> tempString >> numPoints >> tempString;
        points.resize(numPoints);
        for (int i = 0; i < numPoints; i++) {
            file >> points[i].coords[0] >> points[i].coords[1] >> points[i].coords[2];
        }
        getline(file, tempString);
        file >> tempString >> numElements >> tempInt;
        mesh.resize(numElements);
        for (int i = 0; i < numElements; i++) {
            file >> tempInt >> mesh[i].connectivity[0] >> mesh[i].connectivity[1] >> mesh[i].connectivity[2] >> mesh[i].connectivity[3] >>
                mesh[i].connectivity[4] >> mesh[i].connectivity[5] >> mesh[i].connectivity[6] >> mesh[i].connectivity[7];
            for (int j = 0; j < 8; j++) {
                mesh[i].points[j][0] = points[mesh[i].connectivity[j]].coords[0];
                mesh[i].points[j][1] = points[mesh[i].connectivity[j]].coords[1];
                mesh[i].points[j][2] = points[mesh[i].connectivity[j]].coords[2];
            }
        }
        for (int i = 0; i < numElements + 5; i++) getline(file, tempString); // Skip lines
        for (int i = 0; i < numPoints; i++) file >> points[i].marker;
        file.close();
        PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}


void Settings::InitializeCondition(int dof, vector<double>& initialVelocity, vector<double>& initialPressure, vector<CoordPoint3D>& points, const vector<array<double, 3>>& velocityNodes, const int dimension) {
    double initialValues[2] = {0.0, 0.0}; // Initial Velocity and Pressure
    initialVelocity.clear();
    initialPressure.clear();
    initialVelocity.resize(dof * dimension, initialValues[0]);
    initialPressure.resize(dof, initialValues[1]);
}
void Settings::LoadVelocityField(string filename, int numPoints, vector<array<double, 3>>& velocity) {
    ifstream file;
    file.open(filename);
    velocity.resize(numPoints);
    
    if (file.is_open()) {
        for (int i = 0; i < numPoints; i++) {
            file >> velocity[i][0] >> velocity[i][1] >> velocity[i][2];
        }
        file.close();
        PetscPrintf(PETSC_COMM_WORLD, "Velocity Field Loaded!\n");
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}

void Settings::DistributeProcessors(string filename, int& numElements, vector<vector<int>>& elementProcessorMap) {
    int elementID = 0;
    string tempString;
    ifstream file;
    file.open(filename, ios::in);
    
    if (file.is_open()) {
        while (!file.eof() && file.peek() != EOF) {
            int temp;
            file >> temp;
            elementProcessorMap[temp].push_back(elementID);
            elementID++;
            file.get();
        }
        numElements = elementID;
        PetscPrintf(PETSC_COMM_WORLD, "Mesh partition finished!\n");
        file.close();
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}
