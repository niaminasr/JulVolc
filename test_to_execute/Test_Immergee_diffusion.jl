# Include common files and import packages

include("../commons.jl")
include("../meshtools.jl")
include("../assemb2D.jl")
include("../inputoutput.jl")

using Plots

##########################################################
###  Parameters and Mesh
##########################################################

MeshFile = "mesh/carre-fin.mat"  # test in the unit square
PenalConst = 1.e5

rho = 0.1 # Radius of the immersed circle
thetastep = 0.00005# immersed boundary

# Define Lamé coefficients
Young = 2.5 # Modulus of rigidity
nu = 0.25 # Lamé coefficient


#########################################################
###   Reading the mesh and defining the immersed boundary
#########################################################

println("Reading the mesh")
TheMesh = volcmesh(MeshFile,1)


MeshFileTrou = "../mesh/carretrou-fin.mat"
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


# Functions for the RHS and associated exact solution
function f(x::RealVec = [0., 0.])::RealNum
    vecteur = x .- [0.5, 0.5]  # Element-wise subtraction

    if norm(vecteur) < rho
        return (0.0)
    else
        return (1.0)
    end
end

function g(x::RealVec = [0., 0.])::RealNum
    return(3.0) 
end

function K(x::RealVec = [0., 0.])::RealNum
    vecteur = x .- [0.5, 0.5]  # Element-wise subtraction

    if norm(vecteur) < rho
        return (0.0)
    else
        return (1.0)
    end
end



#########################################################
###   Building the system
#########################################################

println("Assemblage rigidité")
A = asm_stiff_2D_P1(TheMesh, K)

# # println("Assemblage Masse frontière")
Mb = asm_bnd_mass_2D_P1(TheMesh, 0)
A = A + PenalConst*Mb


println("Assemblage Masse et second membre")

M = asm_mass_2D_P1(TheMesh)
F = M*(meshvector(TheMesh,f))

M_gamma = asm_immersed_bnd_mass_2D_P1(TheMesh, Gamma)
G = M_gamma*(meshvector(TheMesh,g))



F = F + G 

Xsol = (A + 1e-16*I) \ F
U = Xsol


println("Solve the system and computation of error")

y = TheMesh.nodes[2,:]
x = TheMesh.nodes[1,:]

U_outside = zeros(size(U)) # Initialize a vector of zeros
for i = 1:length(x)
    U_outside[i] = K([x[i]; y[i]]) * U[i]  # Correct Julia indexing with []
end



#########################################
# Plots
#########################################
# p1 = surface(x,y,U_outside,title="Approximate solution",camera = (0,90))
# savefig("Test_diffusion_immersed.png")
plot2Dmesh(TheMesh,"../results/mech_immersed.png")


U_interpolate = Interpolate_in_another_mesh(TheMeshTrou, TheMesh, U_outside)
U_trou = read_vector_from_tex("trousolution.tex")
U_diff = (U_interpolate - U_trou)
NormU = sqrt.((U_diff .^2))
println(maximum(NormU))

y_t = TheMeshTrou.nodes[2,:]
x_t = TheMeshTrou.nodes[1,:]
println(size(y_t), size(U_diff), size(U_trou), size(U_interpolate))
p = surface(x_t,y_t,NormU,title="Exact solution",camera = (0,90))
# display(plot(p))#,layout = (3,1))
# savefig("Difference_norm.png")




