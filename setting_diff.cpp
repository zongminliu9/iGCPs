#include "setting_diff.h"

SettingDiff::SettingDiff() 
{
    // 构造函数
}

void SettingDiff::SetParameters(const std::string& filename, std::vector<double>& parameters, std::string& outputPath) 
{
    parameters.resize(12);
    std::ifstream file(filename);
    if (file.is_open()) 
    {
        std::stringstream ss;
        std::string tmp;
        for (int i = 0; i < 12; i++) 
        {
            file >> tmp >> parameters[i];
            if (i < 11) {
                ss << tmp << parameters[i] << "_";
            } else {
                ss << tmp << parameters[i];
            }
        }
        outputPath += ss.str();
        file.close();
        PetscPrintf(PETSC_COMM_WORLD, "Parameters Loaded!\n");

        if (access(outputPath.c_str(), 0) == -1) 
        {
            CreateOutputDirectory(outputPath);
        }
    } 
    else 
    {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}

void SettingDiff::CreateOutputDirectory(const std::string& outputPath) 
{
    PetscPrintf(PETSC_COMM_WORLD, "%s is not existing\n", outputPath.c_str());
    
#ifdef _WIN32
    int flag = mkdir(outputPath.c_str());
#else
    int flag = mkdir(outputPath.c_str(), 0777);
#endif

    if (flag == 0) {
        PetscPrintf(PETSC_COMM_WORLD, "%s created successfully\n", outputPath.c_str());
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "%s creation failed\n", outputPath.c_str());
    }
}

void SettingDiff::InitializeConditions(const std::vector<Vertex3D>& vertices, const std::vector<double>& parameters, std::vector<double>& N0, std::vector<double>& Nplus, std::vector<double>& Nminus) 
{
    int numVertices = vertices.size();
    N0.resize(numVertices, 0.0);
    Nplus.resize(numVertices, 0.0);
    Nminus.resize(numVertices, 0.0);

    for (int i = 0; i < numVertices; i++) 
    {
        if (vertices[i].label == 1) 
        {
            N0[i] = parameters[9];
            Nplus[i] = parameters[10];
            Nminus[i] = parameters[11];
        }
    }
}

void SettingDiff::LoadMesh(const std::string& filename, std::vector<Vertex3D>& vertices, std::vector<Element3D>& elements) 
{
    std::ifstream file(filename);
    if (file.is_open()) 
    {
        int numVertices, numElements, tmp;
        std::string line;
        
        // 跳过前4行
        for (int i = 0; i < 4; i++) getline(file, line);

        file >> line >> numVertices >> line;
        vertices.resize(numVertices);
        for (int i = 0; i < numVertices; i++) 
        {
            file >> vertices[i].coordinates[0] >> vertices[i].coordinates[1] >> vertices[i].coordinates[2];
        }

        file >> line >> numElements >> tmp;
        elements.resize(numElements);
        for (int i = 0; i < numElements; i++) 
        {
            file >> tmp;
            for (int j = 0; j < 8; j++) 
            {
                file >> elements[i].IEN[j];
                elements[i].points[j][0] = vertices[elements[i].IEN[j]].coordinates[0];
                elements[i].points[j][1] = vertices[elements[i].IEN[j]].coordinates[1];
                elements[i].points[j][2] = vertices[elements[i].IEN[j]].coordinates[2];
            }
        }

        // 读取顶点标签
        for (int i = 0; i < numVertices; i++) 
        {
            file >> vertices[i].label;
        }
        file.close();
        PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
    } 
    else 
    {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}

void SettingDiff::LoadVelocityField(const std::string& filename, int numVertices, std::vector<std::array<double, 3>>& velocityField) 
{
    std::ifstream file(filename);
    velocityField.resize(numVertices);
    if (file.is_open()) 
    {
        for (int i = 0; i < numVertices; i++) 
        {
            file >> velocityField[i][0] >> velocityField[i][1] >> velocityField[i][2];
        }
        file.close();
        PetscPrintf(PETSC_COMM_WORLD, "Velocity Field Loaded!\n");
    } 
    else 
    {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}

void SettingDiff::AssignProcessorToElements(const std::string& filename, int& numElements, std::vector<std::vector<int>>& elementProcessorMap) 
{
    std::ifstream file(filename);
    int index = 0, tmp;
    if (file.is_open()) 
    {
        while (!file.eof() && file.peek() != EOF) 
        {
            file >> tmp;
            elementProcessorMap[tmp].push_back(index++);
            file.get(); // 读取换行符
        }
        numElements = index;
        PetscPrintf(PETSC_COMM_WORLD, "Mesh partitioning completed!\n");
        file.close();
    } 
    else 
    {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
    }
}
void SettingDiff::InitializeBoundaryConditions(const std::string& filename, std::vector<Vertex3D>& vertices, std::vector<double>& NplusBoundary, std::vector<double>& NminusBoundary)
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", filename.c_str());
        return;
    }
    
    int numVertices = vertices.size();
    NplusBoundary.resize(numVertices, 0.0);
    NminusBoundary.resize(numVertices, 0.0);

    for (int i = 0; i < numVertices; i++) 
    {
        // 根据需要加载边界条件
        if (vertices[i].label == 1) 
        {
            file >> NplusBoundary[i] >> NminusBoundary[i];
        }
    }
    PetscPrintf(PETSC_COMM_WORLD, "Boundary Conditions Loaded!\n");
    file.close();
}

void SettingDiff::ExportMesh(const std::string& outputPath, const std::vector<Vertex3D>& vertices, const std::vector<Element3D>& elements)
{
    std::ofstream outFile(outputPath);
    if (!outFile.is_open()) 
    {
        PetscPrintf(PETSC_COMM_WORLD, "Cannot open output file %s!\n", outputPath.c_str());
        return;
    }
    
    outFile << "# vtk DataFile Version 2.0\n";
    outFile << "Generated Mesh\n";
    outFile << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
    
    outFile << "POINTS " << vertices.size() << " float\n";
    for (const auto& vertex : vertices) 
    {
        outFile << vertex.coordinates[0] << " " << vertex.coordinates[1] << " " << vertex.coordinates[2] << "\n";
    }

    outFile << "\nCELLS " << elements.size() << " " << elements.size() * 9 << "\n";
    for (const auto& element : elements) 
    {
        outFile << "8 " << element.IEN[0] << " " << element.IEN[1] << " " << element.IEN[2] << " " << element.IEN[3]
                << " " << element.IEN[4] << " " << element.IEN[5] << " " << element.IEN[6] << " " << element.IEN[7] << "\n";
    }
    
    outFile << "\nCELL_TYPES " << elements.size() << "\n";
    for (size_t i = 0; i < elements.size(); ++i) 
    {
        outFile << "12\n";  // Assuming hexahedron elements
    }

    outFile.close();
    PetscPrintf(PETSC_COMM_WORLD, "Mesh Exported to %s!\n", outputPath.c_str());
}

void SettingDiff::ProcessInputData(const std::string& filename, std::vector<Vertex3D>& vertices, std::vector<Element3D>& elements)
{
    // Example to encapsulate multiple steps: reading mesh, loading conditions, etc.
    ReadMesh(filename + "mesh.vtk", vertices, elements);
    PetscPrintf(PETSC_COMM_WORLD, "Data Processing Complete!\n");
}
