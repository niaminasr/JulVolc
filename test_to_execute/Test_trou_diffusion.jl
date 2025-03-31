# Include common files and import packages

include("../commons.jl")
include("../meshtools.jl")
include("../assemb2D.jl")
include("../inputoutput.jl")

using Plots
#using GRUtils

MeshFile = "../mesh/carretrou-fin.mat" #unit square with a hole
PenalConst = 1.e14


# Functions for the RHS and associated exact solution
function f(x::RealVec = [0., 0.])::RealNum
    return(1.0)
end

function g(x::RealVec = [0., 0.])::RealNum
    return(3.0)
end

function K(x::RealVec = [0., 0.])::RealNum
    return(1.0)
end


function Indicatrice(x::RealVec = [0.; 0.])::Real
    vecteur = x .- [0.5, 0.5]  # Element-wise subtraction

    if norm(vecteur) < rho
        return 0.0
    else
        return 1.0
    end
end

println("Lecture du maillage")

TheMesh = volcmesh(MeshFile,1)

println("Assemblage rigidité")

A = asm_stiff_2D_P1(TheMesh, K)

println("Assemblage Masse frontière")

Mdiric = asm_bnd_mass_2D_P1(TheMesh,2) + asm_bnd_mass_2D_P1(TheMesh,4) +
        asm_bnd_mass_2D_P1(TheMesh,1) + asm_bnd_mass_2D_P1(TheMesh,3)
A = A + PenalConst*Mdiric


println("Assemblage Masse et second membre")
M = asm_mass_2D_P1(TheMesh)
F = M*meshvector(TheMesh,f)

Mneumann = asm_bnd_mass_2D_P1(TheMesh,5) +asm_bnd_mass_2D_P1(TheMesh,6) +
        asm_bnd_mass_2D_P1(TheMesh,7) +asm_bnd_mass_2D_P1(TheMesh,8)
G = Mneumann*meshvector(TheMesh, g)

F = F + G

println("Solve the system and computation of error")

println("Il y a ", length(F), " d.d.l")
U = \(A,F)


y = TheMesh.nodes[2,:]
x = TheMesh.nodes[1,:]


p1 = surface(x,y,U,title="Approximate solution",camera = (0,90))
display(plot(p1))

savefig("../results/Test_diffusion.png")
plot2Dmesh(TheMesh,"../results/mech.png")

print(maximum(U))

write_vector_to_tex(U, "../results/trousolution.tex")





