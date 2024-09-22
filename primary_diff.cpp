#include <iostream>
#include <vector>
#include <array>
#include "primarydiff.h"  // 原来是 Transport.h
#include "settingdiff.h"  // 原来是 UserSetting.h
#include <sstream>
#include <iomanip>
#include "time.h"

using namespace std;

static char description[] = "Solve customized transport equation\n";

int main(int argc, char **argv)
{
    if (argc == 3)
    {
        stringstream inputStream, processStream;
        string inputPath;
        inputStream << argv[1];
        inputStream >> inputPath;
        processStream << argv[2];
        int totalProcesses = atoi(argv[2]);

        int rank, nProcs;
        PetscErrorCode err;

        err = PetscInitialize(&argc, &argv, nullptr, description);
        if (err) return err;

        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

        int bzMeshCount;
        vector<double> variables;
        vector<Vertex3D> controlPoints;
        vector<array<double, 3>> velocityField;
        vector<vector<int>> processElements;
        vector<Element3D> elements;
        vector<double> N0_init, Nplus_init, Nminus_init;
        processElements.resize(nProcs);

        string controlMeshFile(inputPath + "controlmesh.vtk");
        string bzMeshFile(inputPath + "bzmeshinfo.txt.epart." + processStream.str());
        string velocityFieldFile(inputPath + "velocityfield.txt");        
        string paramFile(inputPath + "simulation_parameters.txt");
        string outputPath(inputPath);

        SettingDiff *userConfig = new SettingDiff;
        userConfig->LoadMesh(controlMeshFile, controlPoints, elements);
        userConfig->LoadVelocityField(velocityFieldFile, controlPoints.size(), velocityField);
        userConfig->AssignProcesses(bzMeshFile, bzMeshCount, processElements);
        userConfig->InitializeVariables(paramFile, variables, outputPath);
        userConfig->SetInitialConditions(controlPoints, variables, N0_init, Nplus_init, Nminus_init);

        PrimaryDiff *solver = new PrimaryDiff;
        solver->InitializeSimulation(bzMeshCount, velocityField, N0_init, Nplus_init, Nminus_init, variables);
        solver->DistributeElements(processElements);
        solver->RunSimulation(controlPoints, velocityField, elements, inputPath, outputPath);

        PetscPrintf(PETSC_COMM_WORLD, "Simulation complete!\n");
        err = PetscFinalize();
        CHKERRQ(err);
    }
    else if (argc > 3)
    {
        cout << "Too many arguments.\n";
    }
    else
    {
        cout << "Two arguments expected.\n";
    }
    return 0;
}
