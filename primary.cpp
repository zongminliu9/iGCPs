
#include "FlowStationary.h"   
#include "UserConfig.h"       
#include "petscmat.h"
#include <sstream>

#include <vector>
#include <array>
#include "framework.h"        
#include <iostream>
#include <iomanip>
#include "time.h"

using namespace std;

static char help[] = "Solve 3D steady-state Navier-Stokes Equation\n";

int main(int argc, char **argv)
{
    if (argc == 3) 
    {
        stringstream pathStream, numProcsStream;
        string filePath;
        pathStream << argv[1];
        pathStream >> filePath;
        numProcsStream << argv[2];
        int numProcesses = atoi(argv[2]);

        int rank, numProcs;
        PetscErrorCode petscErr;
        /// Initialize PETSc
        petscErr = PetscInitialize(&argc, &argv, (char*)0, help); if (petscErr) return petscErr;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

        int spaceDim(3), numElements;
        vector<double> simParams;		
        vector<array<double, 3>> velocityField;
        vector<FrameworkVertex> controlPoints;
        vector<FrameworkElement> tetrahedralMesh;
        vector<vector<int>> elementDistribution;
        vector<double> initialVelocities, initialPressures;
        time_t startTime, endTime;
        elementDistribution.resize(numProcs);

        /// Define file paths for simulation parameters and mesh
        string controlMeshFile(filePath + "controlmesh.vtk");
        string bezMeshFile(filePath + "bzmeshinfo.txt.epart." + numProcsStream.str());
        string velocityFieldFile(filePath + "initial_velocityfield.txt");		
        string simParamsFile(filePath + "simulation_parameter.txt");

        /// Configure user settings
        UserConfig *userConfig = new UserConfig;
        userConfig->loadSimulationParameters(simParamsFile, simParams);
        userConfig->loadMesh(controlMeshFile, controlPoints, tetrahedralMesh);
        userConfig->loadVelocityField(velocityFieldFile, controlPoints.size(), velocityField);
        userConfig->distributeElements(bezMeshFile, numElements, elementDistribution);
        userConfig->initializeConditions(controlPoints.size(), initialVelocities, initialPressures, controlPoints, velocityField, spaceDim);

        /// Solve 3D steady-state Navier-Stokes problem
        FlowStationary* navierStokesSolver = new FlowStationary;
        navierStokesSolver->initializeProblem(controlPoints.size(), numElements, initialVelocities, initialPressures, simParams);
        navierStokesSolver->assignElementProcessors(elementDistribution);
        navierStokesSolver->runSolver(controlPoints, tetrahedralMesh, velocityField, filePath);
        
        PetscPrintf(PETSC_COMM_WORLD, "Simulation completed successfully!\n");
        petscErr = PetscFinalize(); CHKERRQ(petscErr);
    }
    else if (argc > 3) {
        cout << "Too many arguments provided.\n";
    }
    else {
        cout << "Expected two arguments: <file_path> <num_processes>.\n";
    }
    return 0;
}
