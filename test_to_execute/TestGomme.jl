# Load a simple mesh of a square
# the upper boundary is labeled 3, the lower 1
# the sides are 2 (right) and 4 (left)

include("../commons.jl")
include("../meshtools.jl")
include("../assemb2D.jl")

using Plots
#using GRUtils

MeshFile = "../mesh/carre-fin.mat" #"test.mat" #"carre.mat"  # test in the unit square
PenalConst = 1.e14


toto = volcmesh(MeshFile,1)  

# Force Function

function force(x::RealVec = [0., 0.])::RealVec
    if (x[1]>= 0.4)&&(x[1]<=0.6)
        return([0. , -1.])
    else
        return([0.,0.])
    end
end


println("Nous avons ", nbnodes(toto), "noeuds soit ", 2*nbnodes(toto), " d.d.l")


# Boundary mass for part 1 and 2

println("Assemblage masses frontières")


Mdiric = asm_bnd_mass_vect_2D_P1(toto,1)
Mneumann = asm_bnd_mass_vect_2D_P1(toto,3)



# Stiffness Matrix
println("Assemblage rigidité")


Young = 2.5
nu = 0.25  # yields lambda = mu = 1 in terms of Lamé coef

K = asm_elast_stiff_2D_P1(toto,Young,nu)


# Global Mass and global force

# M = asm_mass_2D_P1_2DOF(toto)
# F = M * meshvector(toto,force,2) # Body force

F = Mneumann * meshvector(toto,force,2) # ground force


# Solve
println("Résolution (en mode bourrinator)")


A = K + PenalConst*Mdiric
U = \(A,F)


# A few graphics


zaza = move2Dmesh(toto,U,1.) # Move the points of the mesh according to U

y = zaza.nodes[2,:]
x = zaza.nodes[1,:]

U1 = U[1:2:end]
U2 = U[2:2:end]

# quiver(x,y,quiver=(U1,U2), legend=false, arrowsize = 0.15)
# plot2Dmesh!(zaza) # superimpose the mesh plot on the quiver graph

# plot2Dmesh(zaza) # superimpose the mesh plot on the quiver graph