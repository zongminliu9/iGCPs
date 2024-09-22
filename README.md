# iGCPs
We introduce an generalized neurocytoskeletal-PDEs (iGCPs) model, incorporating regulatory protein interaction networks for quantitative spatiotemporal simulations in Neurocytoskeleton Dynamics..

2D_to_initialswc.py is used for 生成关键节点的轴突空间参数 for construction of 3D neurocytoskeleton 空间
NeurocytoManifold.m用于将离散化的空间参数整合为空间manifold。
generate_mesh.m用于生成MeshStructure.vtk和velocity_field.txt

make thbspline.mk来生成可执行文件thespline，b样条
使用metis来划分网格，用来并行计算，减少计算时间
make velocity.mk来生成稳态速度场和压力场
make convectiondiff.mk来生成浓度场

B 1 BPAG1 density
N 250 NF-L density
rho 100 Microtubule density
D 0.6 Diffusion coefficient
vplus 0.3 Initial anterograde velocity
vminus 1.2 Initial retrograde velocity
alphap 1.0 Efficiency of Kinsein in mediating mitochondria to MTs 1.0
alpham 0.0 Efficiency of Dynein in mediating mitochondria to MTs 1.0
k'plus 0.1 Release rate of mitochondria from Kinsein on MTs during transport 
k'minus 0.0 Release rate of mitochondria from Dynein on MTs during transport 
dt 0.1 Time step
nstep 300 Total step number
betap 1.0 Density of mitochondria at neuron body
betam 1.0 Density of mitochondria at axon distal end
