# Include common files and import packages

include("../commons.jl")
include("../meshtools.jl")
include("../assemb2D.jl")

using Plots

using Images

##########################################################
###  Parameters and Mesh
##########################################################

MeshFile = "../mesh/carre-extra-fin.mat"  # test in the unit square
PenalConst = 1.e5

rho = 0.1 # Radius of the immersed circle
thetastep = 0.05# immersed boundary

# Define Lamé coefficients
Young = 2.5 # Modulus of rigidity
nu = 0.25 # Lamé coefficient


#########################################################
###   Reading the mesh and defining the immersed boundary
#########################################################

println("Reading the mesh")
TheMesh = volcmesh(MeshFile,1)

println("Définition de la frontière immergée")

global Xbound = [0.;0.]
for t in 0:thetastep:(2*pi-thetastep)
    global Xbound = hcat(Xbound,rho*[cos(t);sin(t)])
end
Xbound = copy(Xbound[:,2:end])
Xbound = Xbound .+ 0.5 # Centrage 

Gamma = reservoir2D(Xbound,TheMesh)
local_remesh_2Dline(Gamma)

println("the immersed boundary is created")

##########################################################
###  Functions for the RHS and associated exact solution
##########################################################

function Indicatrice(x::RealVec = [0.; 0.])::Real
    vecteur = x .- [0.5, 0.5]  # Element-wise subtraction

    if norm(vecteur) < rho
        return 0.0
    else
        return 1.0
    end
end

function f_test(x::RealVec = [0.; 0.])::RealVec
    return [1., 1.]
end

#########################################################
###   Building the system
#########################################################


println("Stiffness assembly for elasticity ")
K = asm_elast_stiff_2D_P1(TheMesh, Young, nu)

# We have Dirichlet on the boundary of the square and 
# we have a Neumann on the immersed boundary

println("Mass assembly for the boundary")
Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 0)
K = K  + PenalConst*Mb
# Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 1)
# K = K  + PenalConst*Mb
# Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 2)
# K = K  + PenalConst*Mb
# Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 3)
# K = K  + PenalConst*Mb
# Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 4)
# K = K  + PenalConst*Mb

# Assemble cross mass matrix
B_grad = asm_cross_mass_tensor_2D_P1P0(TheMesh, Gamma, Young, nu)

println("Mass and right-hand side assembly")
M_gamma = asm_bnd_mass_matrix_on_crack_2DOF_P0(Gamma)
F_gamma =  M_gamma*(midlinevector(Gamma, f_test, 2))

M = asm_mass_vect_2D_P1(TheMesh)
F = M*(meshvector(TheMesh,f_test, 2))

# Build big matrix & RHS
nddl = length(F)
nmult = size(B_grad,2)
ntotal = nddl + nmult
A = spzeros(ntotal,ntotal)

# Insert stiffness matrix with penalized Dirichlet condition
A[1:nddl,1:nddl] = K
A[1:nddl,(nddl+1):end] = -B_grad
A[(nddl+1):end,1:nddl] = -copy(transpose(B_grad))

# Build the RHS 
RHS = zeros(ntotal,1)
RHS[1:nddl] = F 
RHS[nddl+1: ntotal] = -10*F_gamma
#########################################################
###   Solving the system
#########################################################

Xsol = (A + 1e-16*I) \ RHS
U = Xsol[1:nddl]

#########################################################
###   Solving the system
#########################################################

y = TheMesh.nodes[2,:]
x = TheMesh.nodes[1,:]
U1 = U[1:2:end]
U2 = U[2:2:end]
quiver(x,y,quiver=(U1,U2), legend=false, arrowsize = 0.15)
savefig("../results/Test_imm.png")

plot2Dlinewithnormals(Gamma)


U_outside = zeros(nddl) # Initialize a vector of zeros
for i = 1:length(x)
    U_outside[2*i-1: 2*i] = Indicatrice([x[i]; y[i]]) * U[2*i-1: 2i]  # Correct Julia indexing with []
end

U_outside_1 = U_outside[1:2:end]
U_outside_2 = U_outside[2:2:end]

NormU = sqrt.((U_outside_1.^2) + (U_outside_2.^2))
print(maximum(NormU))

lam = B_grad * Xsol[nddl+1: ntotal]

# surface(x, y, lam(1:2:end), camera = (0,90))
p = surface(x, y, NormU, camera = (0,90))
savefig(p, "../results/surface_plot_gomme_immersed.png")

# plot2Dmesh!(zaza) # superimpose the mesh plot on the quiver graph
# immersed_sol =  move2Dmesh(TheMesh,U,1.)
# plot2Dmeshsave(immersed_sol)




