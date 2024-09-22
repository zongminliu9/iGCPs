#include "ConvectionDiffusion.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

const double PI = 3.1415926;

double CalculateDeterminant(double matrix[3][3]) // 改名为 CalculateDeterminant
{
    double determinant_value = matrix[0][0] * matrix[1][1] * matrix[2][2] 
                             + matrix[0][1] * matrix[1][2] * matrix[2][0] 
                             + matrix[0][2] * matrix[2][1] * matrix[1][0] 
                             - (matrix[0][2] * matrix[1][1] * matrix[2][0] 
                             + matrix[0][0] * matrix[1][2] * matrix[2][1] 
                             + matrix[1][0] * matrix[0][1] * matrix[2][2]);
    return determinant_value;
}

void Inverse3x3Matrix(double input_matrix[3][3], double output_matrix[3][3]) // 改名为 Inverse3x3Matrix
{
    double determinant = CalculateDeterminant(input_matrix);
    output_matrix[0][0] = 1 / determinant * (input_matrix[1][1] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][1]);
    output_matrix[0][1] = 1 / determinant * (input_matrix[2][1] * input_matrix[0][2] - input_matrix[0][1] * input_matrix[2][2]);
    output_matrix[0][2] = 1 / determinant * (input_matrix[0][1] * input_matrix[1][2] - input_matrix[1][1] * input_matrix[0][2]);
    output_matrix[1][0] = 1 / determinant * (input_matrix[2][0] * input_matrix[1][2] - input_matrix[1][0] * input_matrix[2][2]);
    output_matrix[1][1] = 1 / determinant * (input_matrix[0][0] * input_matrix[2][2] - input_matrix[0][2] * input_matrix[2][0]);
    output_matrix[1][2] = 1 / determinant * (input_matrix[1][0] * input_matrix[0][2] - input_matrix[0][0] * input_matrix[1][2]);
    output_matrix[2][0] = 1 / determinant * (input_matrix[1][0] * input_matrix[2][1] - input_matrix[1][1] * input_matrix[2][0]);
    output_matrix[2][1] = 1 / determinant * (input_matrix[0][1] * input_matrix[2][0] - input_matrix[0][0] * input_matrix[2][1]);
    output_matrix[2][2] = 1 / determinant * (input_matrix[0][0] * input_matrix[1][1] - input_matrix[0][1] * input_matrix[1][0]);
}

ConvectionDiffusionSolver::ConvectionDiffusionSolver()
{
    communicator = MPI_COMM_WORLD;
    mpi_error = MPI_Comm_rank(communicator, &rank);
    mpi_error = MPI_Comm_size(communicator, &process_count);
    total_processes = process_count;
    initialized = 0;
}

void ConvectionDiffusionSolver::LoadBezierElementData(string filepath) // 改名为 LoadBezierElementData
{
    string temp_string;
    int num_points, num_elements, num_functions, temp_int, element_count = 0;

    string matrix_filename = filepath + "matrix_data.txt";
    ifstream matrix_file;
    matrix_file.open(matrix_filename);

    if (matrix_file.is_open()) {
        matrix_file >> num_elements;
        elements_in_process.resize(process_elements.size());

        for (int i = 0; i < num_elements; i++) {
            if (i == process_elements[element_count]) {
                matrix_file >> temp_int >> num_functions >> elements_in_process[element_count].element_type;
                elements_in_process[element_count].control_matrix.resize(num_functions);
                elements_in_process[element_count].node_indices.resize(num_functions);

                for (int j = 0; j < num_functions; j++)
                    matrix_file >> elements_in_process[element_count].node_indices[j];

                for (int j = 0; j < num_functions; j++) {
                    for (int k = 0; k < 64; k++) {
                        matrix_file >> elements_in_process[element_count].control_matrix[j][k];
                    }
                }
                element_count++;
            } else {
                matrix_file >> temp_string >> num_functions >> temp_int;
                for (int j = 0; j < num_functions; j++)
                    matrix_file >> temp_string;
                for (int j = 0; j < num_functions; j++)
                    for (int k = 0; k < 64; k++)
                        matrix_file >> temp_string;
            }
        }
        matrix_file.close();
        PetscPrintf(PETSC_COMM_WORLD, "Control Matrices Loaded!\n");
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", matrix_filename.c_str());
    }

    string point_filename = filepath + "bezier_points.txt";
    ifstream point_file;
    point_file.open(point_filename);
    element_count = 0;

    if (point_file.is_open()) {
        point_file >> num_points;
        getline(point_file, temp_string);

        for (int e = 0; e < num_elements; e++) {
            if (e == process_elements[element_count]) {
                elements_in_process[element_count].control_points.resize(control_point_count);

                for (int i = 0; i < control_point_count; i++) {
                    point_file >> elements_in_process[element_count].control_points[i][0] 
                               >> elements_in_process[element_count].control_points[i][1] 
                               >> elements_in_process[element_count].control_points[i][2];
                }
                element_count++;
            } else {
                for (int i = 0; i < control_point_count; i++)
                    point_file >> temp_string >> temp_string >> temp_string;
            }
        }
        point_file.close();
        PetscPrintf(PETSC_COMM_WORLD, "Control Points Loaded!\n");
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", point_filename.c_str());
    }
}
void ConvectionDiffusionSolver::SetGaussPoints(int num_points) // 改名为 SetGaussPoints
{
    gauss_points.clear();
    gauss_weights.clear();

    switch (num_points)
    {
    case 2:
        gauss_points = {0.211324865, 0.788675135};
        gauss_weights = {1.0, 1.0};
        break;
    case 3:
        gauss_points = {0.112701665, 0.5, 0.887298335};
        gauss_weights = {0.55555556, 0.88888889, 0.55555556};
        break;
    case 4:
        gauss_points = {0.069431844, 0.330009478, 0.669990522, 0.930568156};
        gauss_weights = {0.347854845, 0.652145155, 0.652145155, 0.347854845};
        break;
    case 5:
        gauss_points = {0.046910077, 0.230765345, 0.5, 0.769234655, 0.953089923};
        gauss_weights = {0.236926885, 0.478628671, 0.568888889, 0.478628671, 0.236926885};
        break;
    default:
        gauss_points = {0.211324865, 0.788675135};
        gauss_weights = {1.0, 1.0};
        break;
    }
}

void ConvectionDiffusionSolver::InitializeSystem(const int num_elements, vector<array<double, 3>> &velocity_field, const vector<double>& initial_n0, const vector<double>& initial_nplus, const vector<double>& initial_nminus, const vector<double>& parameters)
{
    MPI_Barrier(communicator);

    PetscPrintf(PETSC_COMM_WORLD, "Starting Initialization...\n");

    SetGaussPoints(4);
    num_mesh_elements = num_elements;

    // 初始化各个变量
    n0.resize(initial_n0.size());
    nplus.resize(initial_nplus.size());
    nminus.resize(initial_nminus.size());
    vplus.resize(initial_n0.size());
    vminus.resize(initial_n0.size());

    if (!parameters.empty())
        param_values = parameters;
    else
        cerr << "No parameter values provided!\n";

    time_step = parameters[7];
    num_steps = parameters[8];

    n0 = initial_n0;
    nplus = initial_nplus;
    nminus = initial_nminus;

    PetscInt matrix_size = initial_n0.size() * 2;

    ierr = MatCreate(PETSC_COMM_WORLD, &global_matrix);
    ierr = MatSetSizes(global_matrix, PETSC_DECIDE, PETSC_DECIDE, matrix_size, matrix_size);
    ierr = MatSetType(global_matrix, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(global_matrix, 250, NULL, 250, NULL);
    ierr = MatSetUp(global_matrix);
    MatGetOwnershipRange(global_matrix, &row_start, &row_end);

    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrix_size, &global_rhs);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, matrix_size, &solution);
}

void ConvectionDiffusionSolver::ComputeSUPG(double velocity[3], vector<array<double, 3>>& gradient_N, double &tau_SUPG)
{
    double tau1(0.), tau2(0.);

    for (size_t i = 0; i < gradient_N.size(); i++)
        tau1 += fabs(velocity[0] * gradient_N[i][0] + velocity[1] * gradient_N[i][1] + velocity[2] * gradient_N[i][2]);

    tau1 = 1.0 / tau1;
    tau2 = time_step / 2.0;
    tau_SUPG = 1.0 / sqrt(1.0 / (tau1 * tau1) + 1.0 / (tau2 * tau2));
}

void ConvectionDiffusionSolver::ComputeBasisFunctions(double u, double v, double w, const int num_nodes, const vector<Vertex3D> &control_points, const vector<array<double, 64>> &cmat, vector<double> &basis_N, vector<array<double, 3>> &gradient_N, double dudx[3][3], double &jacobian_det)
{
    double Nu[4] = {(1 - u) * (1 - u) * (1 - u), 3 * (1 - u) * (1 - u) * u, 3 * (1 - u) * u * u, u * u * u};
    double Nv[4] = {(1 - v) * (1 - v) * (1 - v), 3 * (1 - v) * (1 - v) * v, 3 * (1 - v) * v * v, v * v * v};
    double Nw[4] = {(1 - w) * (1 - w) * (1 - w), 3 * (1 - w) * (1 - w) * w, 3 * (1 - w) * w * w, w * w * w};
    double dNdu[4] = {-3 * (1 - u) * (1 - u), 3 - 12 * u + 9 * u * u, 3 * (2 - 3 * u) * u, 3 * u * u};
    double dNdv[4] = {-3 * (1 - v) * (1 - v), 3 - 12 * v + 9 * v * v, 3 * (2 - 3 * v) * v, 3 * v * v};
    double dNdw[4] = {-3 * (1 - w) * (1 - w), 3 - 12 * w + 9 * w * w, 3 * (2 - 3 * w) * w, 3 * w * w};

    double Nt_bz[bzpt_num];
    double dNdt_bz[bzpt_num][3];

    basis_N.clear();
    gradient_N.clear();
    basis_N.resize(num_nodes, 0);
    gradient_N.resize(num_nodes, {0});

    int loc = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                Nt_bz[loc] = Nu[k] * Nv[j] * Nw[i];
                dNdt_bz[loc][0] = dNdu[k] * Nv[j] * Nw[i];
                dNdt_bz[loc][1] = Nu[k] * dNdv[j] * Nw[i];
                dNdt_bz[loc][2] = Nu[k] * Nv[j] * dNdw[i];
                loc++;
            }
        }
    }

    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < bzpt_num; j++) {
            basis_N[i] += cmat[i][j] * Nt_bz[j];
            for (int m = 0; m < 3; m++) {
                gradient_N[i][m] += cmat[i][j] * dNdt_bz[j][m];
            }
        }
    }

    double dxdt[3][3] = {{0}};
    for (int loc = 0; loc < num_nodes; loc++) {
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                dxdt[a][b] += control_points[loc].coor[a] * gradient_N[loc][b];
            }
        }
    }

    double dtdx[3][3] = {{0}};
    ComputeInverseMatrix(dxdt, dtdx);

    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < 3; j++) {
            gradient_N[i][j] = dNdt_bz[i][j] * dtdx[0][0] + dNdt_bz[i][1] * dtdx[1][0] + dNdt_bz[i][2] * dtdx[2][0];
        }
    }

    jacobian_det = ComputeDeterminant(dxdt);
    jacobian_det = 0.125 * jacobian_det;
}

void ConvectionDiffusionSolver::AssembleResidual(const int num_nodes, const double N0_value, const double Nplus_value, const vector<double> &basis_N, const vector<double> &weighted_basis, const double detJ, vector<double> &residual_vector)
{
    for (int i = 0; i < num_nodes; i++) {
        residual_vector[i] += N0_value * basis_N[i] * detJ;
        residual_vector[i + num_nodes] += Nplus_value * weighted_basis[i] * detJ;
    }
}

void ConvectionDiffusionSolver::AssembleMatrix(const int num_nodes, vector<double>& basis_N, vector<double>& weighted_basis, vector<array<double, 3>>& gradient_N, double velocity_plus[3], double velocity_minus[3], double detJ, vector<vector<double>>& element_matrix)
{
    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < num_nodes; j++) {
            element_matrix[i][j] += (basis_N[i] * basis_N[j] + time_step * diffusion_coefficient * (gradient_N[i][0] * gradient_N[j][0] + gradient_N[i][1] * gradient_N[j][1] + gradient_N[i][2] * gradient_N[j][2])) * detJ;
            element_matrix[i][j + num_nodes] += -reaction_rate * time_step * basis_N[i] * basis_N[j] * detJ;
        }
    }
}

void ConvectionDiffusionSolver::ApplyBoundaryCondition(const double boundary_value, int point_index, int variable_index, vector<vector<double>>& element_matrix, vector<double>& residual_vector)
{
    int total_nodes = residual_vector.size() / 2;

    for (int i = 0; i < total_nodes * 2; i++) {
        residual_vector[i] -= boundary_value * element_matrix[i][point_index + variable_index * total_nodes];
    }

    for (int i = 0; i < total_nodes * 2; i++) {
        element_matrix[i][point_index + variable_index * total_nodes] = 0.0;
        element_matrix[point_index + variable_index * total_nodes][i] = 0.0;
    }

    element_matrix[point_index + variable_index * total_nodes][point_index + variable_index * total_nodes] = 1.0;
    residual_vector[point_index + variable_index * total_nodes] = boundary_value;
}

void ConvectionDiffusionSolver::MatrixAssembly(vector<vector<double>>& element_matrix, const vector<int>& connectivity, Mat& global_matrix)
{
    int num_nodes = connectivity.size();
    PetscInt *node_indices = new PetscInt[num_nodes * 2];
    PetscReal *local_matrix_values = new PetscReal[num_nodes * num_nodes * 4];

    int index = 0;
    for (int i = 0; i < num_nodes; i++) {
        int global_index = connectivity[i];
        for (int j = 0; j < 2; j++) {
            node_indices[2 * i + j] = 2 * global_index + j;
            for (int k = 0; k < num_nodes; k++) {
                for (int l = 0; l < 2; l++) {
                    local_matrix_values[index] = element_matrix[i + j * num_nodes][k + l * num_nodes];
                    index++;
                }
            }
        }
    }
    MatSetValues(global_matrix, num_nodes * 2, node_indices, num_nodes * 2, node_indices, local_matrix_values, ADD_VALUES);
    delete[] node_indices;
    delete[] local_matrix_values;
}

void ConvectionDiffusionSolver::ResidualAssembly(vector<double>& element_residual, const vector<int>& connectivity, Vec& global_residual)
{
    int num_nodes = connectivity.size();
    PetscInt *node_indices = new PetscInt[num_nodes * 2];
    PetscReal *local_residual_values = new PetscReal[num_nodes * 2];

    int index = 0;
    for (int i = 0; i < num_nodes; i++) {
        int global_index = connectivity[i];
        node_indices[2 * i] = 2 * global_index;
        node_indices[2 * i + 1] = 2 * global_index + 1;
        local_residual_values[index] = element_residual[i];
        local_residual_values[index + 1] = element_residual[i + num_nodes];
        index += 2;
    }

    VecSetValues(global_residual, num_nodes * 2, node_indices, local_residual_values, ADD_VALUES);
    delete[] node_indices;
    delete[] local_residual_values;
}



void ConvectionDiff::CalculateWeighting(const double flow_velocity[3], const double& scalar, const double& tau_param, const vector<double>& shape_function, const vector<array<double, 3>>& grad_shape_function, vector<double>& weight_function)
{
    /* Compute weighting function based on SUPG stabilization parameter Tau */
    double velocity_magnitude = sqrt(flow_velocity[0] * flow_velocity[0] + flow_velocity[1] * flow_velocity[1] + flow_velocity[2] * flow_velocity[2]);
    for (int idx = 0; idx < shape_function.size(); idx++)
    {
        weight_function[idx] = shape_function[idx] + tau_param * (flow_velocity[0] * grad_shape_function[idx][0] +
                        flow_velocity[1] * grad_shape_function[idx][1] + flow_velocity[2] * grad_shape_function[idx][2]);
    }
}

void ConvectionDiff::CalculateElementValue(const vector<double>& shape_function, const vector<double>& node_values, double& interpolated_value)
{
    interpolated_value = 0.0;
    for (int idx = 0; idx < shape_function.size(); idx++)
    {
        interpolated_value += node_values[idx] * shape_function[idx];
    }
}

void ConvectionDiff::CalculateElementVelocity(const vector<double>& shape_function, const vector<array<double, 3>>& node_velocities, double interpolated_velocity[3])
{
    for (int j = 0; j < dim; j++)
    {
        interpolated_velocity[j] = 0;
    }
    for (int i = 0; i < shape_function.size(); i++)
    {
        for (int j = 0; j < dim; j++)
        {
            interpolated_velocity[j] += node_velocities[i][j] * shape_function[i];
        }
    }
}

void ConvectionDiff::ComputeTangentMatrix(const int num_elements, vector<double>& shape_function, vector<double>& supg_plus, vector<double>& supg_minus, vector<array<double, 3>>& grad_shape_function, double pos_velocity[3], double neg_velocity[3], double jacobian_det, vector<vector<double>>& stiffness_matrix)
{
    int i, j;
    for (i = 0; i < num_elements; i++)
    {
        for (j = 0; j < num_elements; j++)
        {
            stiffness_matrix[i][j] += ((1 + dt * (params[3] + params[4])) * shape_function[i] * shape_function[j] +
                                       dt * params[0] * (grad_shape_function[i][0] * grad_shape_function[j][0] +
                                                         grad_shape_function[i][1] * grad_shape_function[j][1] +
                                                         grad_shape_function[i][2] * grad_shape_function[j][2])) * jacobian_det;
            stiffness_matrix[i][j + num_elements] += -params[5] * dt * shape_function[i] * shape_function[j] * jacobian_det;
            stiffness_matrix[i + num_elements][j] += -params[3] * dt * supg_plus[i] * shape_function[j] * jacobian_det;
            stiffness_matrix[i + num_elements][j + num_elements] += (1 + dt * params[5]) * supg_plus[i] * shape_function[j] * jacobian_det +
                                                                     dt * supg_plus[i] * (pos_velocity[0] * grad_shape_function[j][0] +
                                                                                          pos_velocity[1] * grad_shape_function[j][1] +
                                                                                          pos_velocity[2] * grad_shape_function[j][2]) * jacobian_det;
        }
    }
}

void ConvectionDiff::ComputeResidual(const int num_elements, const double N0, const double N_plus, const double N_minus, const vector<double>& shape_function, const vector<double>& supg_plus, const vector<double>& supg_minus, const double jacobian_det, vector<double>& residual_vector)
{
    for (int i = 0; i < num_elements; i++)
    {
        residual_vector[i] += N0 * shape_function[i] * jacobian_det;
        residual_vector[i + num_elements] += N_plus * supg_plus[i] * jacobian_det;
    }
}

void ConvectionDiff::ApplyBoundaryConditions(const double boundary_value, int element_idx, int variable_idx, vector<vector<double>>& stiffness_matrix, vector<double>& residual_vector)
{
    int num_elements = residual_vector.size() / 2;
    for (int j = 0; j < num_elements * 2; j++)
    {
        residual_vector[j] -= boundary_value * stiffness_matrix[j][element_idx + variable_idx * num_elements];
    }
    for (int j = 0; j < num_elements * 2; j++)
    {
        stiffness_matrix[j][element_idx + variable_idx * num_elements] = 0.0;
        stiffness_matrix[element_idx + variable_idx * num_elements][j] = 0.0;
    }
    stiffness_matrix[element_idx + variable_idx * num_elements][element_idx + variable_idx * num_elements] = 1.0;
    residual_vector[element_idx + variable_idx * num_elements] = boundary_value;
}

void ConvectionDiff::AssembleMatrix(vector<vector<double>>& stiffness_matrix, const vector<int>& element_indices, Mat& global_matrix)
{
    int num_elements = element_indices.size();
    PetscInt* index_list = new PetscInt[num_elements * 2];
    PetscReal* local_matrix = new PetscReal[num_elements * 2 * num_elements * 2];
    int add = 0;
    for (int i = 0; i < element_indices.size(); i++)
    {
        int A = element_indices[i];
        for (int j = 0; j < 2; j++)
        {
            index_list[2 * i + j] = 2 * A + j;
            for (int k = 0; k < element_indices.size(); k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    local_matrix[add] = stiffness_matrix[i + j * num_elements][k + l * num_elements];
                    add++;
                }
            }
        }
    }
    MatSetValues(global_matrix, num_elements * 2, index_list, num_elements * 2, index_list, local_matrix, ADD_VALUES);
    delete[] index_list;
    delete[] local_matrix;
}

void ConvectionDiff::AssembleResidual(vector<double>& residual_vector, const vector<int>& element_indices, Vec& global_vector)
{
    int num_elements = element_indices.size();
    PetscInt* index_list = new PetscInt[num_elements * 2];
    PetscReal* local_residual = new PetscReal[num_elements * 2];
    int add = 0;
    for (int i = 0; i < element_indices.size(); i++)
    {
        int A = element_indices[i];
        index_list[2 * i + 0] = A * 2 + 0;
        index_list[2 * i + 1] = A * 2 + 1;
        local_residual[add] = residual_vector[i];
        local_residual[add + 1] = residual_vector[i + num_elements];
        add += 2;
    }
    VecSetValues(global_vector, num_elements * 2, index_list, local_residual, ADD_VALUES);
    delete[] index_list;
    delete[] local_residual;
}

void ConvectionDiff::MTPlusCalculation(double x, double y, double z, double& result, double result_grad[3])
{
    result = 1.0;
    result_grad[0] = 0.0;
    result_grad[1] = 0.0;
    result_grad[2] = 0.0;
}



void Solver::ComputeMatrix(const int nodeCount, vector<double> &shapeFunc, vector<double> &shapeFuncP, vector<double> &shapeFuncM, vector<array<double, dimensions>> &gradShapeFunc, double velPos[dimensions], double velNeg[dimensions], double linParam, double gradParam[dimensions], double detJacobian, vector<vector<double>> &globalMatrix)
{
    int row, col, idxA, idxB;
    for (row = 0; row < nodeCount; row++) {
        for (col = 0; col < nodeCount; col++) {
            globalMatrix[row + 0][col + 0] += ((1 + timeStep * (parameters[3] + parameters[4])) * shapeFunc[row] * shapeFunc[col] + timeStep * parameters[0] * (gradShapeFunc[row][0] * gradShapeFunc[col][0] + gradShapeFunc[row][1] * gradShapeFunc[col][1] + gradShapeFunc[row][2] * gradShapeFunc[col][2])) * detJacobian;
            globalMatrix[row + 0][col + nodeCount] += -parameters[5] * timeStep * linParam * shapeFunc[row] * shapeFunc[col] * detJacobian;
            globalMatrix[row + nodeCount][col + 0] += -parameters[3] * timeStep * shapeFuncP[row] * shapeFunc[col] * detJacobian;
            globalMatrix[row + nodeCount][col + nodeCount] += (1 + timeStep * parameters[5]) * linParam * shapeFuncP[row] * shapeFunc[col] * detJacobian + timeStep * shapeFuncP[row] * (velPos[0] * (linParam * gradShapeFunc[col][0] + shapeFunc[col] * gradParam[0]) + velPos[1] * (linParam * gradShapeFunc[col][1] + shapeFunc[col] * gradParam[1]) + velPos[2] * (linParam * gradShapeFunc[col][2] + shapeFunc[col] * gradParam[2])) * detJacobian;
        }
    }
}

void Solver::ComputeResidual(const int nodeCount, const double sol0, const double solP, const double solM, const vector<double> &shapeFunc, const vector<double> &shapeFuncP, const vector<double> &shapeFuncM, double linParam, double gradParam[dimensions], const double detJacobian, vector<double> &globalResidual)
{
    int row;
    for (row = 0; row < nodeCount; row++) {
        globalResidual[row] += sol0 * shapeFunc[row] * detJacobian;
        globalResidual[row + nodeCount] += linParam * solP * shapeFuncP[row] * detJacobian;
    }
}

void Solver::CalculateConcentration(double u, double v, double w, const Element &currentElement, double point[3], double &displacement, double gradU[3], double &detJacobian)
{
    double jacobianMat[3][3];
    vector<double> shapeFunc(currentElement.IEN.size());
    vector<array<double, 3>> gradShapeFunc(currentElement.IEN.size());
    vector<Vertex> controlPoints;
    controlPoints.resize(currentElement.IEN.size());

    for (int i = 0; i < currentElement.IEN.size(); i++)
        controlPoints[i] = controlPointsList[currentElement.IEN[i]];

    ShapeFunctions(u, v, w, currentElement.IEN.size(), controlPoints, currentElement.shapeMatrix, shapeFunc, gradShapeFunc, jacobianMat, detJacobian);
    
    displacement = 0.0;
    point[0] = point[1] = point[2] = 0.0;
    gradU[0] = gradU[1] = gradU[2] = 0.0;

    for (uint i = 0; i < currentElement.IEN.size(); i++) {
        point[0] += shapeFunc[i] * controlPoints[i].coords[0];
        point[1] += shapeFunc[i] * controlPoints[i].coords[1];
        point[2] += shapeFunc[i] * controlPoints[i].coords[2];

        displacement += shapeFunc[i] * (sol0[currentElement.IEN[i]] + solPlus[currentElement.IEN[i]]);
        gradU[0] += gradShapeFunc[i][0] * (sol0[currentElement.IEN[i]] + solPlus[currentElement.IEN[i]]);
        gradU[1] += gradShapeFunc[i][1] * (sol0[currentElement.IEN[i]] + solPlus[currentElement.IEN[i]]);
        gradU[2] += gradShapeFunc[i][2] * (sol0[currentElement.IEN[i]] + solPlus[currentElement.IEN[i]]);
    }
}

void Solver::ExportVTK_ControlMesh(const vector<Vertex> &vertexList, const vector<Element> &elementList, int step, string fileName)
{
    ofstream outFile;
    stringstream stepStr;
    stepStr << step;
    fileName += "/control_mesh_" + stepStr.str() + ".vtk";
    outFile.open(fileName.c_str());
    if (outFile.is_open()) {
        outFile << "# vtk DataFile Version 2.0\nMesh Data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        outFile << "POINTS " << vertexList.size() << " float\n";
        for (unsigned int i = 0; i < vertexList.size(); i++) {
            outFile << vertexList[i].coords[0] << " " << vertexList[i].coords[1] << " " << vertexList[i].coords[2] << "\n";
        }
        outFile << "\nCELLS " << elementList.size() << " " << 9 * elementList.size() << '\n';
        for (unsigned int i = 0; i < elementList.size(); i++) {
            outFile << "8 " << elementList[i].IEN[0] << " " << elementList[i].IEN[1] << " " << elementList[i].IEN[2] << " " << elementList[i].IEN[3]
                    << " " << elementList[i].IEN[4] << " " << elementList[i].IEN[5] << " " << elementList[i].IEN[6] << " " << elementList[i].IEN[7] << '\n';
        }
        outFile << "\nCELL_TYPES " << elementList.size() << '\n';
        for (unsigned int i = 0; i < elementList.size(); i++) {
            outFile << "12\n";
        }
        outFile << "POINT_DATA " << solPlus.size() << "\nSCALARS Displacement float 1\nLOOKUP_TABLE default\n";
        for (uint i = 0; i < solPlus.size(); i++) {
            outFile << sol0[i] + solPlus[i] + solMinus[i] << "\n";
        }
        outFile.close();
    } else {
        cout << "Unable to open file: " << fileName << endl;
    }
}


void Solver::ExportVTK_Domain(int step, string outputDir)
{
    vector<array<double, 3>> samplePoints;
    vector<double> sampleResults;
    vector<array<int, 8>> sampleElements;
    double jacobianDet;
    int localElementCount = elementProcessList.size();
    double pointProcessData[localElementCount * 24];
    double resultProcessData[localElementCount * 8];
    int elementProcessData[localElementCount * 8];

    for (unsigned int e = 0; e < localElementCount; e++)
    {
        int resolution = 2;
        vector<double> subDomain(resolution);
        for (int i = 0; i < resolution; i++)
        {
            subDomain[i] = double(i) / (double(resolution) - 1.);
        }

        int pointIdx = 0;
        for (int a = 0; a < resolution; a++)
        {
            for (int b = 0; b < resolution; b++)
            {
                for (int c = 0; c < resolution; c++)
                {
                    double point[3], gradU[3];
                    double result;
                    CalculateConcentration(subDomain[c], subDomain[b], subDomain[a], elementProcessList[e], point, result, gradU, jacobianDet);
                    pointProcessData[24 * e + pointIdx * 3 + 0] = point[0];
                    pointProcessData[24 * e + pointIdx * 3 + 1] = point[1];
                    pointProcessData[24 * e + pointIdx * 3 + 2] = point[2];
                    resultProcessData[8 * e + pointIdx] = result;
                    pointIdx++;
                }
            }
        }

        int elementNodes[2] = {resolution * resolution * resolution, resolution * resolution};
        for (int a = 0; a < resolution - 1; a++)
        {
            for (int b = 0; b < resolution - 1; b++)
            {
                for (int c = 0; c < resolution - 1; c++)
                {
                    elementProcessData[8 * e + 0] = 8 * e + a * elementNodes[1] + b * resolution + c;
                    elementProcessData[8 * e + 1] = 8 * e + a * elementNodes[1] + b * resolution + c + 1;
                    elementProcessData[8 * e + 2] = 8 * e + a * elementNodes[1] + (b + 1) * resolution + c + 1;
                    elementProcessData[8 * e + 3] = 8 * e + a * elementNodes[1] + (b + 1) * resolution + c;
                    elementProcessData[8 * e + 4] = 8 * e + (a + 1) * elementNodes[1] + b * resolution + c;
                    elementProcessData[8 * e + 5] = 8 * e + (a + 1) * elementNodes[1] + b * resolution + c + 1;
                    elementProcessData[8 * e + 6] = 8 * e + (a + 1) * elementNodes[1] + (b + 1) * resolution + c + 1;
                    elementProcessData[8 * e + 7] = 8 * e + (a + 1) * elementNodes[1] + (b + 1) * resolution + c;
                }
            }
        }
    }

    double *gatheredPoints = NULL;
    double *gatheredResults = NULL;
    int *gatheredElements = NULL;
    int *displacementsPoints = NULL;
    int *displacementsResults = NULL;
    int *displacementsElements = NULL;
    int *elementCounts = NULL;
    int *gatherCountsPoints = NULL;
    int *gatherCountsResults = NULL;
    int *gatherCountsElements = NULL;

    if (comRank == 0)
    {
        elementCounts = (int*)malloc(sizeof(int)*nProcess);
        gatherCountsPoints = (int*)malloc(sizeof(int)*nProcess);
        gatherCountsResults = (int*)malloc(sizeof(int)*nProcess);
        gatherCountsElements = (int*)malloc(sizeof(int)*nProcess);
    }

    MPI_Gather(&localElementCount, 1, MPI_INT, elementCounts, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Barrier(comm);

    if (comRank == 0)
    {
        gatheredPoints = (double*)malloc(sizeof(double) * 24 * totalElements);
        gatheredResults = (double*)malloc(sizeof(double) * 8 * totalElements);
        gatheredElements = (int*)malloc(sizeof(int) * 8 * totalElements);

        displacementsPoints = (int*)malloc(nProcess * sizeof(int));
        displacementsResults = (int*)malloc(nProcess * sizeof(int));
        displacementsElements = (int*)malloc(nProcess * sizeof(int));

        displacementsPoints[0] = 0;
        displacementsResults[0] = 0;
        displacementsElements[0] = 0;

        for (int i = 1; i < nProcess; i++)
        {
            displacementsPoints[i] = displacementsPoints[i - 1] + elementCounts[i - 1] * 24;
            displacementsResults[i] = displacementsResults[i - 1] + elementCounts[i - 1] * 8;
            displacementsElements[i] = displacementsElements[i - 1] + elementCounts[i - 1] * 8;
        }

        for (int i = 0; i < nProcess; i++)
        {
            gatherCountsPoints[i] = elementCounts[i] * 24;
            gatherCountsResults[i] = elementCounts[i] * 8;
            gatherCountsElements[i] = elementCounts[i] * 8;
        }
    }

    MPI_Gatherv(pointProcessData, localElementCount * 8 * 3, MPI_DOUBLE, gatheredPoints, gatherCountsPoints, displacementsPoints, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Gatherv(resultProcessData, localElementCount * 8, MPI_DOUBLE, gatheredResults, gatherCountsResults, displacementsResults, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Gatherv(elementProcessData, localElementCount * 8, MPI_INT, gatheredElements, gatherCountsElements, displacementsElements, MPI_INT, 0, PETSC_COMM_WORLD);

    if (comRank == 0)
    {
        for (int i = 0; i < totalElements; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                array<double, 3> point = {gatheredPoints[i * 24 + j * 3 + 0], gatheredPoints[i * 24 + j * 3 + 1], gatheredPoints[i * 24 + j * 3 + 2]};
                samplePoints.push_back(point);
                sampleResults.push_back(gatheredResults[i * 8 + j]);
            }
        }

        int sumElements = 0;
        int pointStart = 0;
        for (int i = 0; i < nProcess; i++)
        {
            for (int e = 0; e < elementCounts[i]; e++)
            {
                array<int, 8> element;
                element[0] = pointStart + gatheredElements[8 * sumElements + 0];
                element[1] = pointStart + gatheredElements[8 * sumElements + 1];
                element[2] = pointStart + gatheredElements[8 * sumElements + 2];
                element[3] = pointStart + gatheredElements[8 * sumElements + 3];
                element[4] = pointStart + gatheredElements[8 * sumElements + 4];
                element[5] = pointStart + gatheredElements[8 * sumElements + 5];
                element[6] = pointStart + gatheredElements[8 * sumElements + 6];
                element[7] = pointStart + gatheredElements[8 * sumElements + 7];
                sampleElements.push_back(element);
                sumElements++;
            }
            pointStart += elementCounts[i] * 8;
        }

        cout << "Writing VTK file for physical domain..." << endl;
        WriteVTKFile(samplePoints, sampleResults, sampleElements, step, outputDir);
    }
}

void Solver::WriteVTKFile(const vector<array<double, 3>> &points, const vector<double> &results, const vector<array<int, 8>> &elements, int step, const string &outputDir)
{
    stringstream fileName;
    fileName << step;
    string fullPath = outputDir + "/physical_mesh_" + fileName.str() + ".vtk";
    ofstream outFile(fullPath.c_str());

    if (outFile.is_open())
    {
        outFile << "# vtk DataFile Version 2.0\nHexahedron Mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        outFile << "POINTS " << points.size() << " float\n";
        for (size_t i = 0; i < points.size(); i++)
        {
            outFile << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
        }

        outFile << "\nCELLS " << elements.size() << " " << 9 * elements.size() << '\n';
        for (size_t i = 0; i < elements.size(); i++)
        {
            outFile << "8 " << elements[i][0] << " " << elements[i][1] << " " << elements[i][2] << " " << elements[i][3] << " "
                    << elements[i][4] << " " << elements[i][5] << " " << elements[i][6] << " " << elements[i][7] << '\n';
        }

        outFile << "\nCELL_TYPES " << elements.size() << '\n';
        for (size_t i = 0; i < elements.size(); i++)
        {
            outFile << "12\n";  // Hexahedron cell type
        }

        outFile << "POINT_DATA " << results.size() << "\nSCALARS scalarField float 1\nLOOKUP_TABLE default\n";
        for (size_t i = 0; i < results.size(); i++)
        {
            outFile << results[i] << "\n";
        }

        outFile.close();
    }
    else
    {
        cout << "Unable to open file " << fullPath << " for writing.\n";
    }
}
