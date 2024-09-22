#ifndef FLOW_STATIONARY_H
#define FLOW_STATIONARY_H

#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include "framework.h"
#include "time.h"

using namespace std;

const int flow_degree = 3;
const int space_dim = 3;
const int control_points = 64;

class FlowStationary
{
private:
	vector<double> gaussPoints;
	vector<double> weights;
	const double PI = 4 * atan(1.0);

	PetscErrorCode petscErr;
	MPI_Comm mpiComm;
	int mpiError;
	int rank;
	int size;
	int numProcesses;

	int rowStart, rowEnd;
	int bezMeshSize;
	vector<int> elementProcessors;
	vector<FrameworkElement> bezMeshProcessed;

	KSP kspSolver;
	PC preconditioner;
	Mat systemMatrix;
	Vec residualVector;
	Vec solutionVector;
	
	double timeStep;
	double maxVelocity;
	double viscosity;
	double density;
	double alphaMass;
	double alphaForce;
	double gammaValue;
	double forceX, forceY, forceZ;

	vector<double> velocities;
	vector<double> pressures;
	vector<double> boundaryConditions;
	vector<double> parameters; 

public:
	FlowStationary();

private:
	void computeBodyForce(double x, double y, double z, double &Fx, double &Fy, double &Fz);
	void loadBezierElements(string fileName);

	/* Analysis */
	void initializeGaussInfo(int numGaussPoints);	
	void computeBasisFunctions(double u, double v, double w, int numElements, const vector<array<double, 3>>& controlPoints, const vector<array<double, 64>>& controlMatrix, vector<double>& shapeFunctions, vector<array<double, 3>>& shapeGradients, vector<array<array<double, 3>, 3>>& shapeHessians, double jacobian[3][3], double& determinantJ);
	void computePointValue(vector<double>& shapeFunctions, const vector<array<double, 4>>& elementValues, double result[4]);
	void computePointGradient(vector<array<double, 3>>& shapeGradients, const vector<array<double, 4>>& elementValues, double result[4][3]);
	void computePointHessian(vector<array<array<double, 3>, 3>>& shapeHessians, const vector<array<double, 4>>& elementValues, double result[4][3][3]);
	void computeStabilizationParams(double jacobian[3][3], double velocity[4], double &tauM, double &tauC);
	void computeFineScale(double tauM, double tauC, double velocity[3], double velocityX[3], double velocityY[3], double velocityZ[3], double velocityXX[3], double velocityYY[3], double velocityZZ[3], double pressure, double pressureX, double pressureY, double pressureZ, double fineScaleVelocity[3], double &fineScalePressure);
	void computeResidual(vector<double>& shapeFunctions, vector<array<double, 3>>& shapeGradients, vector<array<array<double, 3>, 3>>& shapeHessians, double jacobian[3][3], const double detJ, const vector<array<double, 4>>& elementValues, vector<array<double, 4>>& residual);
	void computeTangentMatrix(vector<double>& shapeFunctions, vector<array<double, 3>>& shapeGradients, double jacobian[3][3], const double detJ, const vector<array<double, 4>>& elementValues, vector<array<vector<array<double, 4>>, 4>>& tangentMatrix);
	void buildSystemForProcessor(const vector<FrameworkVertex>& controlPoints, const vector<array<double, 3>>& velocityBC, const vector<double> velocityNodes, const vector<double> pressureNodes);
	void applyBoundaryCondition(const double bcValue, int pointIdx, int variableIdx, vector<array<vector<array<double, 4>>, 4>>& tangentMatrix, vector<array<double, 4>>& residual);
	void assembleMatrix(vector<array<vector<array<double, 4>>, 4>>& tangentMatrix, const vector<int>& connectivity, Mat& systemMatrix);
	void assembleResidual(vector<array<double, 4>>& residual, const vector<int>& connectivity, Vec& residualVec);

	/* Postprocessing */
	void computeResultForBezier(double u, double v, double w, const FrameworkElement& element, double point[3], double result[4], double jacobian[3], double& detJ);
	void visualizeControlMesh(const vector<FrameworkVertex>& controlPoints, const vector<FrameworkElement>& mesh, int step, string filename);
	void visualizePhysicalDomain(int step, string filename);
	void writeToVTK(const vector<array<double, 3>>& points, const vector<double> results, const vector<array<int, 8>>& elements, int step, string filename);

public:
	/* Preprocessing */
	void initializeProblem(const int numDOF, const int numElements, const vector<double>& initialVelocities, const vector<double>& initialPressures, const vector<double>& parameters);
	void assignProcessorTasks(vector<vector<int>>& elementDistribution);
	void runSolver(const vector<FrameworkVertex>& controlPoints, const vector<FrameworkElement>& mesh, const vector<array<double, 3>>& velocityBC, string filename);
};

#endif
