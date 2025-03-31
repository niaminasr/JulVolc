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

MeshFile = "../mesh/carre_sans_trou_01.mat"# test in the unit square
# MeshFile = "carre-fin.mat"
PenalConst = 1.e5

rho = 0.1 # Radius of the immersed circle
# rho = 0.0025 # Radius of the immersed circle
thetastep = 0.0005# immersed boundary

# Define Lamé coefficients
Young = 2.5 # Modulus of rigidity
nu = 0.25 # Lamé coefficient
# Young = 16.0 # Modulus of rigidity
# nu = 0.25 # Lamé coefficient

pressure = 0.02
#########################################################
###   Reading the mesh and defining the immersed boundary
#########################################################

println("Reading the mesh")
TheMesh = volcmesh(MeshFile,1)

# println("Définition de la frontière immergée")
# global Xbound = [0.;0.]
# for t in 0:thetastep:(2*pi-thetastep)
#     global Xbound = hcat(Xbound,rho*[cos(t);sin(t)])
# end
# Xbound = copy(Xbound[:,2:end])
# Xbound = Xbound .+ 0.5 # Centrage 

# Gamma = reservoir2D(Xbound,TheMesh)
# local_remesh_2Dline(Gamma)

println("the immersed boundary is created")
s = 0.01
a = 0.7
b = 0.8
n = [-1.0; 0.5]
global Xbound = [0.;0.]

# Générer les points pour x entre x1 et x2
for x in a:thetastep:b
    y_shiftplus = ((x  - a)/2.0 + b) + s*n[2] # Calculer la coordonnée y selon l'équation y = ax + b
    x_shiftplus = x + s*n[1]  
    global Xbound = hcat(Xbound, [x_shiftplus; y_shiftplus])
    
    y_mid_2 = (b - a)/2.0 + b
    x_mid_2 = b
    global Xbound = hcat(Xbound, [x_mid_2; y_mid_2])

    y_shiftmoin = (((x) - a)/2.0 + b) - s*n[2] # Calculer la coordonnée y selon l'équation y = ax + b
    x_shiftmoin = (x) - s*n[1]  
    global Xbound = hcat(Xbound, [x_shiftmoin; y_shiftmoin])

    y_mid_1 = b
    x_mid_1 = a
    global Xbound = hcat(Xbound, [x_mid_1; y_mid_1])
end

println(size(Xbound))
Xbound = copy(Xbound[:,2:end])

Gamma = reservoir2D(Xbound,TheMesh)
local_remesh_2Dline(Gamma)

# ##########################################################
# ###  Functions for the RHS and associated exact solution
# ##########################################################
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

# function g(x::RealVec = [0.; 0.])::RealVec
#     vecteur = (x .- [0.5 0.5])
#     vecteur /= norm(vecteur)

#     return(vecteur)
# end

function d1(x::RealVec = [0.; 0.])::RealNum
    return (x[2] + s/2 - a)/2 + 1/2 + s - x[1]
end

function d2(x::RealVec = [0.; 0.])::RealNum
    return (x[2] - s/2 - a)/2 + 1/2 - s - x[1]
end

function d3(x::RealVec = [0.; 0.])::RealNum
    return 1/2-2*(x[2]-a) - x[1]
end

function d4(x::RealVec = [0.; 0.])::RealNum
    return (b-a)/2+1/2-2*(x[2]-b) - x[1]
end

function g(x::RealVec = [0.; 0.])::RealVec

    vecteur = [1.0; 1.0]

    if( d1(x) == 0)
         vecteur = [1/(sqrt(5)/2); -0.5/(sqrt(5)/2)]
    elseif(d2(x) == 0)
        vecteur =  [-1/(sqrt(5)/2); 0.5/(sqrt(5)/2)]
    elseif(d3(x) == 0)
        vecteur = [-1/(sqrt(5)); -2/(sqrt(5))]
    elseif(d4(x) == 0)
        vecteur = [1/(sqrt(5)); 2/(sqrt(5))]
    end
  
    return vecteur
end

function K_epsilon(x::RealVec = [0., 0.])::RealNum
    vecteur = x .- [0.5, 0.5]  # Element-wise subtraction

    if norm(vecteur) < rho
        return (0.00000001)
    else
        return (1.0)
    end
end

function uexact(x::RealVec = [0., 0.])::RealNum
    
    mu = 1.0
    depth = 0.5
    # deltaV = (pi * pressure * (rho * rho))/ mu
    deltaV = (((2*pi*rho*rho*(1-nu*nu))/Young)*pressure)
    # println(deltaV)
    # deltaV = 7.85e-4
    xi = (x[1] - 0.5) / depth

    return ((2.0 * (1 - nu) * deltaV) / (pi * depth)) * (1 / (1 + (xi*xi)))
end
#########################################################
###   Building the system
#########################################################
println("Stiffness assembly for elasticity ")
K = asm_elast_stiff_2D_P1(TheMesh, Young, nu, K_epsilon)

# We have Dirichlet on the boundary of the square and 
# we have a Neumann on the immersed boundary
println("Mass assembly for the boundary")
Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 1)
K = K  + PenalConst*Mb
Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 2)
K = K  + PenalConst*Mb
Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 4)
K = K  + PenalConst*Mb
# Mb = asm_bnd_mass_vect_2D_P1(TheMesh, 3)
# K = K  + PenalConst*Mb

M = asm_mass_vect_2D_P1(TheMesh)
F = M*(meshvector(TheMesh,f, 2))

println("Mass and right-hand side assembly")
M_gamma = asm_immersed_bnd_mass_2D_P1_2DOF(TheMesh, Gamma)
F_gamma =  M_gamma*((+pressure) .* meshvector(TheMesh, g, 2))
F = F + F_gamma

#########################################################
###   Solving the system
#########################################################
nddl = length(F)

Xsol = (K + 1e-16*I) \ F
U = Xsol[1:nddl]

Ux = U[1:2:end]
Uy = U[2:2:end] 

#########################################################
###  Find the nodes on the Ground Boundary
#########################################################
BndSegments = findall(isequal(3) , TheMesh.edges[3,:])
N = length(BndSegments)
Node_number = Vector{Int}(undef, (N + 1))

for i in 1:N
    TheEdge = TheMesh.edges[:,BndSegments[i]]
    Node_number[i] = TheEdge[1]
end
TheEdge = TheMesh.edges[:,BndSegments[N]]
Node_number[N+1] = TheEdge[2]

println(Node_number)

Node_coordinate = zeros(2, (N + 1))
U_err = zeros((N +1))
for i in 1:(N + 1)
    Node_coordinate[:,i] = TheMesh.nodes[:, Node_number[i]]

    U_err[i] = Uy[Node_number[i]] - uexact(Node_coordinate[:,i])
    #  println(Uy[Node_number[i]], "  ", uexact(Node_coordinate[:,i]))
end

Norm_inf =  norm(U_err,Inf)
println("norm Infinie:::", "  ",  norm(U_err,Inf))
# println(U)
#########################################
###  Plots
#########################################
moved_mesh =  move2Dmesh(TheMesh,U,1.)
p = plot2Dmesh(moved_mesh, "../results/fracture.png")











