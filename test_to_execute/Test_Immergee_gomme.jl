# Include common files and import packages

include("../commons.jl")
include("../meshtools.jl")
include("../assemb2D.jl")
include("../inputoutput.jl")

using Plots

using Images

##########################################################
###  Parameters and Mesh
##########################################################

MeshFile = "../mesh/Mesh_N_320.mat"  # test in the unit square
MeshFileTrou = "../mesh/carretrou-fin.mat"
PenalConst = 1.e5

rho = 0.1 # Radius of the immersed circle
thetastep = 0.0005# immersed boundary

# Define Lamé coefficients
Young = 2.5 # Modulus of rigidity
nu = 0.25 # Lamé coefficient

lambda = 1.0;
mu = 1.0;


#########################################################
###   Reading the mesh and defining the immersed boundary
#########################################################

println("Reading the mesh")
TheMesh = volcmesh(MeshFile,1)
TheMeshTrou = volcmesh(MeshFileTrou, 1)

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

function f(x::RealVec = [0.; 0.])::RealVec
    vecteur = x .- [0.5, 0.5]  # Element-wise subtraction

    if norm(vecteur) < rho
        return ([0.0, 0.0])
    else
        return ([0.0, 0.0])
    end
end

function g(x::RealVec = [0.; 0.])::RealVec
    vecteur = (x - [0.5;0.5])
    vecteur /= norm(vecteur)

    return(vecteur)
end


function K_epsilon(x::RealVec = [0., 0.])::RealNum
    vecteur = x .- [0.5, 0.5]  # Element-wise subtraction

    if norm(vecteur) < rho
        return (0.00000001)
    else
        return (1.0)
    end
end

#########################################################
function uexact(x::RealVec = [0., 0.])::RealVec
    return [(x[1] * (x[1] - 1) * x[2] * (x[2] - 1)); (x[1] * (x[1] - 1) * x[2] * (x[2] - 1))]
end

function grad(x::RealVec = [0., 0.])::RealMat
    return (
            [(2 * x[1] - 1) * x[2] * (x[2] - 1)  (2 * x[2] - 1) * x[1] * (x[1] - 1);
             (2 * x[1] - 1) * x[2] * (x[2] - 1)  (2 * x[2] - 1) * x[1] * (x[1] - 1)]
           )
end

function strain_tensor(x::RealVec = [0., 0.])::RealMat

    # Calculate the gradient 
    grad_f = grad(x)

    # Create the strain tensor
    epsilon = zeros(2, 2)
    epsilon[1, 1] = grad_f[1, 1]  # Normal strain in x-direction
    epsilon[2, 2] = grad_f[2, 2]  # Normal strain in y-direction
    
    # Calculate shear strains
    epsilon[1, 2] = 0.5 * (grad_f[1,2] + grad_f[2,1] )
    epsilon[2, 1] =  epsilon[1, 2]  

    return epsilon
end

function stress_tensor(x::RealVec = [0., 0.])::RealMat

    epsilon = strain_tensor(x)
    trace_epsilon = epsilon[1, 1] + epsilon[2, 2]

    sigma = zeros(2, 2)
    sigma[1, 1] = 2 * mu * epsilon[1, 1] + lambda * trace_epsilon
    sigma[2, 2] = 2 * mu * epsilon[2, 2] + lambda * trace_epsilon
    sigma[1, 2] = 2 * mu * epsilon[1, 2]
    sigma[2, 1] = sigma[1, 2]

    return sigma
end

# Function to compute the divergence of the stress tensor
function divergence_of_stress(x::RealVec = [0., 0.])
    
    # Calculate divergence of the stress tensor
    div_sigma = zeros(2)  # Initialize divergence vector

    # Divergence components
    div_sigma[1] =  -partialsigma11(x) - partialsigma12(x)
    div_sigma[2] =  -partialsigma21(x) - partialsigma22(x)

    return div_sigma
end

# Partial derivatives of stress tensor components
function partialsigma11(x::RealVec = [0., 0.])
    return (4 * mu * x[2] * (x[2] - 1)) + lambda * ((2 * x[2] * (x[2] - 1)) +  (2 * x[1] - 1) * (2 * x[2] - 1))
end

function partialsigma12(x::RealVec = [0., 0.])
    return  mu * ((2 * x[1] * (x[1] - 1)) + (2 * x[2] - 1) * (2 * x[1] - 1))
end

function partialsigma21(x::RealVec = [0., 0.])
    return  mu * ((2 * x[2] - 1) * (2 * x[1] - 1) + (2 * x[2] * (x[2]-1)))
end

function partialsigma22(x::RealVec = [0., 0.])
    return (4 * mu * x[1] * (x[1]  -1)) + lambda * ((2 * x[1] * (x[1] - 1)) +  (2 * x[1] - 1) * (2 * x[2] - 1))
end

function f_manufactured(x::RealVec = [0., 0.])::RealVec
    return(K_epsilon(x).*divergence_of_stress(x))
end

function g_manufactured(x::RealVec)::RealVec
    vecteur = (x - [0.5;0.5])
    vecteur /= norm(vecteur)

    return(-stress_tensor(x)*vecteur) 
end
#########################################################
###   Building the system
#########################################################


println("Stiffness assembly for elasticity ")
K = asm_elast_stiff_2D_P1(TheMesh, Young, nu, K_epsilon)

# We have Dirichlet on the boundary of the square and 
# we have a Neumann on the immersed boundary

println("Mass assembly for the boundary")
Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 0)
K = K  + PenalConst*Mb


M = asm_mass_vect_2D_P1(TheMesh)
F = M*(meshvector(TheMesh, f_manufactured, 2))


println("Mass and right-hand side assembly")
M_gamma = asm_immersed_bnd_mass_2D_P1_2DOF(TheMesh, Gamma)
F_gamma =  M_gamma*(1.0 .* meshvector(TheMesh, g_manufactured, 2))

F = F + F_gamma
#########################################################
###   Solving the system
#########################################################
nddl = length(F)

Xsol = (K + 1e-16*I) \ F
U = Xsol[1:nddl]

Uex = meshvector(TheMesh, uexact, 2)
ErrorU = (U - Uex)

y = TheMesh.nodes[2,:]
x = TheMesh.nodes[1,:]

U_outside = zeros(nddl) 
for i = 1:length(x)
    U_outside[2*i-1: 2*i] = Indicatrice([x[i]; y[i]]) * (ErrorU[2*i-1: 2i]) 
end
println("norm Infinie:::", "  ",  norm(U_outside,Inf))

U_outside_1 = U_outside[1:2:end]
U_outside_2 = U_outside[2:2:end]

NormU_E2 = sqrt.((U_outside_1.^2) + (U_outside_2.^2))
println("norm eucli:::", "  ",maximum(NormU_E2))

p = surface(x, y, NormU_E2, camera = (0,90))
savefig(p, "NormU_E2.png")

NormU_L2 = sqrt(U_outside' * M * U_outside)
println("norm L2:::", "  ",NormU_L2)

K_H1 =  asm_elast_stiff_2D_P1_lambda_mu(TheMesh, 0.0, 1.0, K_epsilon)
NormU_H1 = sqrt(U_outside' * K_H1 * U_outside) + sqrt(U_outside' * M * U_outside)
println("norm H1:::", "  ",NormU_H1)
#########################################################
###   Solving the system
#########################################################

# y = TheMesh.nodes[2,:]
# x = TheMesh.nodes[1,:]
# U1 = U[1:2:end]
# U2 = U[2:2:end]
# quiver(x,y,quiver=(U1,U2), legend=false, arrowsize = 0.15)
# savefig("Test_imm_2.png")

# plot2Dlinewithnormals(Gamma)



# U_outside = zeros(nddl) # Initialize a vector of zeros
# for i = 1:length(x)
#     U_outside[2*i-1: 2*i] = Indicatrice([x[i]; y[i]]) * U[2*i-1: 2i]  # Correct Julia indexing with []
# end

# U_outside_1 = U_outside[1:2:end]
# U_outside_2 = U_outside[2:2:end]

# NormU = sqrt.((U_outside_1.^2) + (U_outside_2.^2))
# print(maximum(NormU))

# # surface(x, y, lam(1:2:end), camera = (0,90))
# p = surface(x, y, NormU, camera = (0,90))
# savefig(p, "surface_plot_gomme_immersed_2.png")

# # plot2Dmesh!(zaza) # superimpose the mesh plot on the quiver graph

# immersed_sol_2 =  move2Dmesh(TheMesh,U,1.)
# plot2Dmeshsave(immersed_sol_2)


# # Interpolate
# U_interpolate = Interpolate_in_another_mesh_2DOF(TheMeshTrou, TheMesh, U_outside)
# U_trou = read_vector_from_tex("trousolution_gomme.tex")

# println("U_interpolate" , size(U_interpolate))
# println("U_trou" , size(U_trou))


# U_diff = (U_interpolate - U_trou)

# U1_diff = U_diff[1:2:end]
# U2_diff = U_diff[2:2:end] 

# U1_trou = U_trou[1:2:end]
# U2_trou = U_trou[2:2:end] 

# U1_imm = U[1:2:end]
# U2_imm = U[2:2:end] 

# NormU = sqrt.((U1_diff.^2) + (U2_diff.^2))

# y_t = TheMeshTrou.nodes[2,:]
# x_t = TheMeshTrou.nodes[1,:]

# y_i = TheMesh.nodes[2,:]
# x_i = TheMesh.nodes[1,:]


# println(maximum(NormU))

# p = surface(x_t,y_t,NormU,title="Exact solution",camera = (0,90))
# p1 = surface(x_t, y_t, U1_trou, titre="U_1trou", camera = (0,90))
# p2 = surface(x_t, y_t, U2_trou, titre="U_2trou", camera = (0,90))

# p3 = surface(x_i, y_i, U1_imm, titre="U_1IMM", camera = (0,90))
# p4= surface(x_i, y_i, U2_imm, titre="U_2IMM", camera = (0,90))


# display(plot(p))
# savefig("Difference_norm_gomme.png")

# display(plot(p1,p3))
# savefig("Disp_gomme_U1.png")

# display(plot(p2,p4))
# savefig("Disp_gomme_U2.png")


# file = MAT.matopen("error.mat", "w")
# MAT.write(file, "x", x_t)
# MAT.write(file, "y", y_t)
# MAT.write(file, "u1", U1_diff)
# MAT.write(file,"u2", U2_diff)
# MAT.write(file, "Norme", NormU)
# MAT.close(file)






