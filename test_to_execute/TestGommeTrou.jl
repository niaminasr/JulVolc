# Load a simple mesh of a square
# the upper boundary is labeled 3, the lower 1
# the sides are 2 (right) and 4 (left)

include("../meshtools.jl")
include("../assemb2D.jl")
include("../inputoutput.jl")

#using Plots
#using GRUtils

MeshFile = "../mesh/carre_trou_003_carre02.mat" #unit square with a hole
PenalConst = 1.e16


toto = volcmesh(MeshFile,1)  

# Force Function on the hole pushes inside or outside (accord. to sign)

function normalforce(x::RealVec = [0.; 0.])::RealVec
    vecteur = (x - [1.0, 1.0])
    vecteur /= norm(vecteur)

    return(0.02 .* vecteur)
end


println("Nous avons ", nbnodes(toto), "noeuds soit ", 2*nbnodes(toto), " d.d.l")


# Boundary mass for part 1 and 2

println("Assemblage masses frontières")

#cotés 1,2,3,4 = tour du carré (encastré)
# cotés 5,6,7,8 = bord du trou (poussé)


Mdiric = asm_bnd_mass_vect_2D_P1(toto,2) + asm_bnd_mass_vect_2D_P1(toto,4) +
        asm_bnd_mass_vect_2D_P1(toto,1) 
        

Mneumann = asm_bnd_mass_vect_2D_P1(toto,5) +asm_bnd_mass_vect_2D_P1(toto,6) +
        asm_bnd_mass_vect_2D_P1(toto,7) +asm_bnd_mass_vect_2D_P1(toto,8)

# Build RHS


#Version 1 : compute the boundary normals on the hole

# Force = zeros(2*nbnodes(toto))
# BndSegments = [ findall(isequal(5) , toto.edges[3,:]); findall(isequal(6) , toto.edges[3,:]) ;
#         findall(isequal(7) , toto.edges[3,:]) ; findall(isequal(8) , toto.edges[3,:]) ]
# NN = mk_edge_normals(toto)[:,BndSegments]



# NN = hcat(NN , NN[:,1]) # add first point at the end
# NN = (NN[:,1:end-1] + NN[:,2:end]) / 2.0

# println(dimension(toto))


# indices = toto.edges[1,BndSegments]


# for iloop in 1:length(indices)
#         jndex = indices[iloop]
#         Force[jndex:(jndex+1)] = NN[:,iloop]
# end

# Force = 5*Force

# Version 2 : use the function normalforce

Force = meshvector(toto,normalforce,2)



F = Mneumann * Force # force on the hole



# Stiffness Matrix

println("Assemblage rigidité")


Young = 2.5
nu = 0.25  # yields lambda = mu = 1 in terms of Lamé coef

K = asm_elast_stiff_2D_P1(toto,Young,nu)


# Global Mass and global force

#M = asm_mass_2D_P1_2DOF(toto)


# F = M * meshvector(toto,force,2) # Body force


# Solve

println("Résolution (en mode bourrinator)")


A = K + PenalConst*Mdiric

U = \(A,F)


zaza = move2Dmesh(toto,U,1.) # Move the points of the mesh according to U

plot2Dmesh(zaza, "../results/moved_mesh_trou.png")





y = toto.nodes[2,:]
x = toto.nodes[1,:]

U1 = U[1:2:end]
U2 = U[2:2:end] 

NormU = sqrt.((U1.^2) + (U2.^2))
print(maximum(NormU))

p = surface(x, y, NormU, camera = (0,90))

savefig(p, "../results/surface_plot_trou.png")





# # y = zaza.nodes[2,:]
# # x = zaza.nodes[1,:]
# # U1 = U[1:2:end]
# # U2 = U[2:2:end]
# # quiver(x,y,quiver=(U1,U2), legend=false, arrowsize = 0.15)
# # plot2Dmesh!(zaza) # superimpose the mesh plot on the quiver graph

# zaza = move2Dmesh(toto,U,1.)

plot2Dmesh(zaza) # superimpose the mesh plot on the quiver graph

write_vector_to_tex(U, "../results/trousolution_gomme.tex")

Ux = U[1:2:end]
Uy = U[2:2:end] 

#########################################################
###  Find the nodes on the Ground Boundary
#########################################################
BndSegments = findall(isequal(3) , toto.edges[3,:])
N = length(BndSegments)
Node_number = Vector{Int}(undef, (N + 1))

for i in 1:N
    TheEdge = toto.edges[:,BndSegments[i]]
    Node_number[i] = TheEdge[1]
end
TheEdge = toto.edges[:,BndSegments[N]]
Node_number[N+1] = TheEdge[2]

println(Node_number)

Node_coordinate = zeros(2, (N + 1))
# U_err = zeros((N +1))
for i in 1:(N + 1)
    Node_coordinate[:,i] = toto.nodes[:, Node_number[i]]

#     U_err[i] = Uy[Node_number[i]] - uexact(Node_coordinate[:,i])
#      println(Uy[Node_number[i]], "  ", uexact(Node_coordinate[:,i]))
end

# Norm_inf =  norm(U_err,Inf)
# println("norm Infinie:::", "  ",  norm(U_err,Inf))


#########################################
###  Plots
#########################################
moved_mesh =  move2Dmesh(toto,U,1.)
p = plot2Dmesh(moved_mesh, "moove.png")

x_t1 = zeros((N+1))
y_t1 = zeros((N+1))
u_num = zeros((N+1))
for i in 1:(N + 1)
    x_t1[i] = Node_coordinate[1,i]
    u_num[i] = Uy[Node_number[i]]
end

u_n = [x_t1, u_num]

println(u_n)

sorted_indices = sortperm(x_t1)  # Get sorted indices
x_t1_sorted = x_t1[sorted_indices]  # Sort x_t1
u_num_sorted = u_num[sorted_indices]


# u_ex_sorted = u_e[sortperm(u_e[:, 1]), :]
# u_num_sorted = u_n[sortperm(u_n[:, 1]), :]

p2 = plot(x_t1_sorted, u_num_sorted, title = "Num solution")

# p1 = plot(x_t1,u_ex ,title="Exact solution")
# p2 = plot(x_t1,u_num ,title="Num solution")
display(plot(p2))
savefig("../results/valerie_sol_2.png")








