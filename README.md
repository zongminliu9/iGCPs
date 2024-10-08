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
