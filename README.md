# iGCPs: Generalized Neurocytoskeletal-PDEs Model

We introduce the iGCPs (Generalized Neurocytoskeletal-PDEs) model, incorporating regulatory protein interaction networks for quantitative spatiotemporal simulations in neurocytoskeleton dynamics.

## Key Components and Scripts:

1. **2D_to_initialswc.py**  
   Used for processing 2D segmentation results to generate spatial parameters for key nodes in axons, which are essential for constructing the 3D neurocytoskeleton space.

2. **NeurocytoManifold.m**  
   Integrates the discretized spatial parameters into a continuous manifold, forming the structure of the neurocytoskeleton.

3. **generate_mesh.m**  
   Generates the `MeshStructure.vtk` 3D spatial file, used for initializing and executing the generalized neurocytoskeletal-PDEs solver with the provided `velocity_field.txt`.

## Building Executables:

- **make thbspline.mk**  
  Generates the executable `thespline` to process B-splines.

- **METIS Grid Partitioning**  
  Use METIS for grid partitioning to enable parallel computation, optimizing runtime performance.

- **make velocity.mk**  
  Generates the executable for calculating the steady-state velocity and pressure fields.

- **make convectiondiff.mk**  
  Generates the executable for solving the concentration field in the 3D neurocytoskeletal space.

## Parameters Utilized in iGCPs:

B         1     # BPAG1 density  
N       250     # NF-L density  
rho     100     # Microtubule density  
D       0.6     # Diffusion coefficient  
vplus   0.3     # Initial anterograde velocity  
vminus  1.2     # Initial retrograde velocity  
alphap  1.0     # Efficiency of kinesin in mediating mitochondria transport along MTs  
alpham  0.0     # Efficiency of dynein in mediating mitochondria transport along MTs  
k'plus  0.1     # Release rate of mitochondria from kinesin during transport on MTs  
k'minus 0.0     # Release rate of mitochondria from dynein during transport on MTs  
dt      0.1     # Time step  
nstep   300     # Total step number  
betap   1.0     # Density of mitochondria at neuron body  
betam   1.0     # Density of mitochondria at axon distal end  


# iGCPs: Generalized Neurocytoskeletal-PDEs Simulation Framework

iGCPs (Generalized Cytoskeletal PDEs) performs fully 3D spatiotemporal simulations of axonal mitochondrial transport by converting 2D neurocytoskeleton images into realistic 3D axonal meshes, solving steady-state Navier–Stokes for intracellular flow, and then simulating advection–diffusion–reaction PDEs that include BPAG1/NF-L regulatory effects.

## Dependencies
- MATLAB R2020a+ (with TREES Toolbox)
- Eigen (C++ linear algebra)
- METIS 5.x (mesh partitioning)
- PETSc 3.18.5 + MPI (PDE solvers)
- ParaView 5.x (visualization)

## Workflow
1. 2D_to_initialswc.py  
   Converts segmented 2D neurocytoskeleton images into `initial.swc` (axon centerlines + swelling nodes).  
   Run:
   ```bash
   python 2D_to_initialswc.py --input ./images/ --output initial.swc
   ```

2. NeurocytoManifold.m  
   Builds a smooth 3D manifold (`neurocyto_geometry.mat`) from `initial.swc`.  
   Run:
   ```matlab
   NeurocytoManifold('initial.swc')
   ```

3. generate_mesh.m  
   Generates an IGA control mesh `MeshStructure.vtk` from `neurocyto_geometry.mat`.  
   Run:
   ```matlab
   generate_mesh('neurocyto_geometry.mat')
   ```

4. THB-Spline Extraction  
   Compile and run Bezier extraction:
   ```bash
   make -f thbspline.mk
   ./thespline mesh/MeshStructure.vtk
   ```

5. Mesh Partitioning  
   Partition mesh for MPI:
   ```bash
   mpmetis bzmeshinfo.txt 28
   ```

6. nsvms (Navier–Stokes Solver)  
   Compile and run:
   ```bash
   make -f velocity.mk
   mpiexec -np 28 ./nsvms mesh/ 28
   ```

7. transport (PDE Solver)  
   Edit `UserSetting.cpp` for initial/boundary conditions, configure `simulation_parameter.txt`, then:
   ```bash
   make -f convectiondiff.mk
   mpiexec -np 28 ./transport mesh/ 28
   ```

## File I/O
- initial.swc → neurocyto_geometry.mat → MeshStructure.vtk → bzpt.txt, cmat.txt, bzmeshinfo.txt → partitioned files → velocityfield.txt → controlmesh_allparticle_*.vtk

## simulation_parameter.txt Example
```
B            1
N          250
rho        100
D          0.6
vplus      0.3
vminus     1.2
alphap     1.0
alpham     0.0
kprimeplus 0.1
kprimeminus0.0
dt         0.1
nstep      300
betap      1.0
betam      1.0
```

## Visualization
Open `.vtk` files in ParaView. Use volume rendering for mitochondrial concentration, glyphs for velocity vectors, and time animation to visualize transport front evolution.

## Citation
If you use iGCPs, please cite:  
Liu Z. et al., “iGCPs: A PDE-based Simulation Framework for Neurocytoskeletal Transport in Neurons,” 2025.

