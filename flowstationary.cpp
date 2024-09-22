#include "flowstationary.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;

double Determinant(double transformation[3][3])
{
    double det = transformation[0][0] * transformation[1][1] * transformation[2][2] + 
                 transformation[0][1] * transformation[1][2] * transformation[2][0] + 
                 transformation[0][2] * transformation[2][1] * transformation[1][0] - 
                 (transformation[0][2] * transformation[1][1] * transformation[2][0] + 
                  transformation[0][0] * transformation[1][2] * transformation[2][1] + 
                  transformation[1][0] * transformation[0][1] * transformation[2][2]);

    return det;
}

void InverseMatrix(double transformation[3][3], double inverse[3][3])
{
    double det = Determinant(transformation);
    inverse[0][0] = 1 / det * (transformation[1][1] * transformation[2][2] - transformation[1][2] * transformation[2][1]);
    inverse[0][1] = 1 / det * (transformation[2][1] * transformation[0][2] - transformation[0][1] * transformation[2][2]);
    inverse[0][2] = 1 / det * (transformation[0][1] * transformation[1][2] - transformation[1][1] * transformation[0][2]);
    inverse[1][0] = 1 / det * (transformation[2][0] * transformation[1][2] - transformation[1][0] * transformation[2][2]);
    inverse[1][1] = 1 / det * (transformation[0][0] * transformation[2][2] - transformation[0][2] * transformation[2][0]);
    inverse[1][2] = 1 / det * (transformation[1][0] * transformation[0][2] - transformation[0][0] * transformation[1][2]);
    inverse[2][0] = 1 / det * (transformation[1][0] * transformation[2][1] - transformation[1][1] * transformation[2][0]);
    inverse[2][1] = 1 / det * (transformation[0][1] * transformation[2][0] - transformation[0][0] * transformation[2][1]);
    inverse[2][2] = 1 / det * (transformation[0][0] * transformation[1][1] - transformation[0][1] * transformation[1][0]);
}

FlowStationary::FlowStationary()
{
    comm = MPI_COMM_WORLD;
    mpiError = MPI_Comm_rank(comm, &comRank);
    mpiError = MPI_Comm_size(comm, &comSize);
    processCount = comSize;	
}


void FlowStationary::ApplyForces(double x, double y, double z, double& forceX, double& forceY, double& forceZ)
{
    forceX = 0;
    forceY = 0;
    forceZ = 0;
}

void FlowStationary::ReadElementData(string filePath)
{
    string elementData, tmp;
    int points, elements, nFunctions, tmpValue;

    ifstream cmatFile(filePath + "matrix.txt");
    if (!cmatFile.is_open()) {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open matrix file!\n");
        return;
    }

    cmatFile >> elements;
    mesh.resize(elementProcessing.size());

    for (int i = 0; i < elements; ++i) {
        if (i == elementProcessing[elementCounter]) {
            cmatFile >> tmpValue >> nFunctions >> mesh[elementCounter].type;
            mesh[elementCounter].matrixCoefficients.resize(nFunctions);
            mesh[elementCounter].functionIndices.resize(nFunctions);

            for (int j = 0; j < nFunctions; ++j) {
                cmatFile >> mesh[elementCounter].functionIndices[j];
            }

            for (int j = 0; j < nFunctions; ++j) {
                for (int k = 0; k < 64; ++k) {
                    cmatFile >> mesh[elementCounter].matrixCoefficients[j][k];
                }
            }
            ++elementCounter;
        } else {
            // Skipping non-processed elements
            for (int j = 0; j < nFunctions; ++j) {
                cmatFile >> tmp;
            }
            for (int j = 0; j < nFunctions; ++j) {
                for (int k = 0; k < 64; ++k) {
                    cmatFile >> tmp;
                }
            }
        }
    }
    cmatFile.close();
    PetscPrintf(PETSC_COMM_WORLD, "Element matrix data loaded!\n");
}

void FlowStationary::InitializeGaussian(int numPoints)
{
    gaussPoints.clear();
    gaussWeights.clear();

    switch (numPoints)
    {
    case 2:
        gaussPoints = {0.2113248654051871, 0.7886751345948129};
        gaussWeights = {1., 1.};
        break;
    case 3:
        gaussPoints = {0.1127016653792583, 0.5, 0.8872983346207417};
        gaussWeights = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
        break;
    case 4:
        gaussPoints = {0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262};
        gaussWeights = {0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539};
        break;
    case 5:
        gaussPoints = {0.046910077030668, 0.2307653449471585, 0.5, 0.7692346550528415, 0.953089922969332};
        gaussWeights = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};
        break;
    default:
        gaussPoints = {0.2113248654051871, 0.7886751345948129};
        gaussWeights = {1., 1.};
    }
}


void FlowStationary::InitializeSystem(int degreesOfFreedom, int numElements, const vector<double>& velocityInitial, const vector<double>& pressureInitial, const vector<double>& parameters)
{
    MPI_Barrier(comm);
    PetscPrintf(PETSC_COMM_WORLD, "Initializing flow system...\n");

    // Initializing Gaussian quadrature with 4 points
    InitializeGaussian(4);
    
    numElementsInMesh = numElements;
    velocity.resize(degreesOfFreedom * 3);
    pressure.resize(degreesOfFreedom);
    
    velocity = velocityInitial;
    pressure = pressureInitial;

    if (!parameters.empty())
    {
        viscosity = 0;
        density = 0.5;
        alphaM = 0.5 * (3 - density) / (1 + density);
        alphaF = 1 / (1 + density);
        gammaValue = 0.5 + alphaM - alphaF;
        maxVelocity = parameters[1];
    }
    else
    {
        cerr << "No parameters found!\n";
        getchar();
    }

    forceX = 0;
    forceY = 0;
    forceZ = 0;

    // Initialize PETSc vector and matrix
    PetscInt matrixSize = degreesOfFreedom * 4;
    ierr = MatCreate(PETSC_COMM_WORLD, &globalMatrix);
    ierr = MatSetSizes(globalMatrix, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize);
    ierr = MatSetType(globalMatrix, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(globalMatrix, 500, NULL, 500, NULL);
    MatGetOwnershipRange(globalMatrix, &rowStart, &rowEnd);
    ierr = MatSetOption(globalMatrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(globalMatrix);

    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, &globalResidual);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrixSize, &tempSolution);
}




void FlowStationary::CalculateBasisFunction(double paramU, double paramV, double paramW, int numFunctions, const vector<array<double, 3>>& points, const vector<array<double, 64>>& coefficients, vector<double>& basis, vector<array<double, 3>>& basisDerivatives, vector<array<array<double, 3>, 3>>& secondDerivatives, double gradientMap[3][3], double& detJacobian)
{
    // Local variables
    double shapeFunctions[4] = { (1 - paramU) * (1 - paramU) * (1 - paramU), 3 * (1 - paramU) * (1 - paramU) * paramU, 3 * (1 - paramU) * paramU * paramU, paramU * paramU * paramU };
    double derivU[4] = { -3 * (1 - paramU) * (1 - paramU), 3 - 12 * paramU + 9 * paramU * paramU, 3 * (2 - 3 * paramU) * paramU, 3 * paramU * paramU };

    // Initialize basis and derivative arrays
    basis.clear();
    basisDerivatives.clear();
    secondDerivatives.clear();
    basis.resize(numFunctions);
    basisDerivatives.resize(numFunctions, { 0.0 });
    secondDerivatives.resize(numFunctions, { {0.0} });

    // Compute basis functions and derivatives
    int index = 0;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                basis[index] = shapeFunctions[k] * shapeFunctions[j] * shapeFunctions[i];
                basisDerivatives[index][0] = derivU[k] * shapeFunctions[j] * shapeFunctions[i];
                basisDerivatives[index][1] = shapeFunctions[k] * derivU[j] * shapeFunctions[i];
                basisDerivatives[index][2] = shapeFunctions[k] * shapeFunctions[j] * derivU[i];
                index++;
            }
        }
    }

    // Compute the gradient map and its determinant
    double jacobian[3][3] = { {0} };
    for (int i = 0; i < numFunctions; ++i)
    {
        for (int a = 0; a < 3; ++a)
        {
            for (int b = 0; b < 3; ++b)
            {
                jacobian[a][b] += points[i][a] * basisDerivatives[i][b];
            }
        }
    }

    Matrix3DInverse(jacobian, gradientMap);
    detJacobian = MatrixDet(jacobian);
}

