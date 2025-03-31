# Julvolc
#
# Routines for 2D problems P1 Elements
#

### TODO : 
# [1] BOUNDARY mass matrix with more then one ddl at each node (elasticity)
# [2] stiffness matrix for operators in the form -div(A(x) Grad u(x)) where A is a MATRIX
# [3] evaluate basis functions at a given point in a given triangle
#   (not documented for users, only for developpers)

# O. Bodart 2022/10/26


# Include common files and import packages

include("commons.jl")

include("meshtools.jl")

import SparseArrays
using SparseArrays


######################################################
######################################################
######################################################
#     P1 elements
######################################################
######################################################
######################################################


###########################################################################
#################                                       ###################
#################        Specific Data Structures       ###################
#################                                       ###################
###########################################################################


mutable struct basisfunvals_2D_P1
    #
    # Data structure containing the triangle numbers in which a point
    # can be and, for each triangle (3 at most in the extreme case
    # where the point is the vertex of a triangle of the mesh)
    #
    elements::IndexVec  # List of the elements where th point lies
    funvalues::RealMat  # For each element values of the 3 basis functions
    isempty::Bool       # flag set to true if the structure contains no value
    #
    #
    # Default constructor (void struct)
    basisfunvals_2D_P1() = new(@EmptyIntVec, @EmptyRealMat,true)
    #
    # Constructor from a 2D point X and a mesh
    #
    #
    function basisfunvals_2D_P1(X::RealVec, TheMesh::volcmesh)
        (ElementZ, ValueZ, Flaggy) = basis_func_pt_val_2D_P1(X, TheMesh)
        return new(ElementZ,ValueZ,Flaggy)
    end


end # of basisfunvals struct definition




###########################################################################
#################                                       ###################
################# Basis functions evaluation routines   ###################
#################                                       ###################
###########################################################################


# First version : returning a tuple to be stored in a basisfunvals_2D_P1 structure

function basis_func_pt_val_2D_P1(X::RealVec, TheMesh::volcmesh)
    # This routine locates the 2D Point X in the mesh and
    # computes the values of the P1 basis functions in the triangles
    # where it belongs. It output a structure of 
    # Input variables:
    # - X : 2D point 
    # - TheMesh : mesh of class volcmesh
    # - BV : basisfunvals_2D_P1 structure, empty or not
    # Output: none
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #

    TheTriangle = zeros(Index, 3)

    # Detect the elements where the point lies.
    # Calling the intriangle function will compute the element surfaces
    # and store them in the mesh if necessary.
    # No extra check is needed here.

    ElemZzz = intriangle(X,TheMesh)

    if (ElemZzz != 0) # if the point X lies in the mesh then compute the values

        NumElemZzz = length(ElemZzz)
        ValueZzz = zeros(RealNum,(3,NumElemZzz))


        for index in 1:NumElemZzz

            # Set basic triangle parameters, number, surface and vertices

            TheTriangle = TheMesh.elements[:,ElemZzz[index]]
            The2Surf = 2*TheMesh.element_measures[ElemZzz[index]]

            x1 = TheMesh.nodes[:, TheTriangle[1]]
            x2 = TheMesh.nodes[:, TheTriangle[2]]
            x3 = TheMesh.nodes[:, TheTriangle[3]]

            # Compute the inverse of the Jacobian matrix of the change of variables
            # Notice that the map from the real space to the reference element is 
            # given by [xi, eta] = InvJac*[x-x1,y-y1]

            InvJac = zeros(RealNum, (2,2))
            InvJac[1,1] = x3[2] - x1[2]
            InvJac[1,2] = x1[1] - x3[1]
            InvJac[2,1] = x1[2] - x2[2]
            InvJac[2,2] = x2[1] - x1[1]

            InvJac /= The2Surf  # The determinant of the jacobian is 2 times the surface of
                                # the triangle under consideration

            # Apply the change of variables

            XiEta = InvJac*(X-x1)

            # Compute the values of the 3 basis functions and store them in the corresponding
            # column of the "funvalues" matrix in BV

            ValueZzz[:,index] = [1 - XiEta[1] - XiEta[2] , XiEta[1], XiEta[2]]
            
        end # of loop over relevant triangles

        return(ElemZzz,ValueZzz,false)

    else # if the point is outside the mesh return empty arrays and the empty flag on
        return(@EmptyIntVec,@EmptyRealMat,true)
        
    end

end # of basis_func_pt_val_2D_P1 v1
    

# Second version : modifies a basisfunvals_2D_P1 structure provided in the input parameters

function basis_func_pt_val_2D_P1(X::RealVec, TheMesh::volcmesh, BV::basisfunvals_2D_P1)
    # This routine locates the 2D Point X in the mesh and
    # computes the values of the P1 basis functions in the triangles
    # where it belongs. It output a structure of 
    # Input variables:
    # - X : 2D point 
    # - TheMesh : mesh of class volcmesh
    # - BV : basisfunvals_2D_P1 structure, empty or not
    # Output: none
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #

    TheTriangle = zeros(Index, 3)

    #
    # Empty the input structure BV
    #
    
    BV.elements = @EmptyIntVec
    BV.funvalues = @EmptyRealMat
    BV.isempty = true

    # Detect the elements where the point lies.
    # Calling the intriangle function will compute the element surfaces
    # and store them in the mesh if necessary.
    # No extra check is needed here.

    ElemZzz = intriangle(X,TheMesh)

    if (ElemZzz == 0) # if the point is not in the mesh return an empty result
        println("ERROR in basis_func_pt_val_2D_P1 - given point not in mesh")
        
    else # if OK then fill BV with the right values
        BV.elements = copy(ElemZzz) # List of the relevant elements
        NumElemZzz = length(ElemZzz)
        BV.funvalues = zeros(RealNum,(3,NumElemZzz))


        for index in 1:NumElemZzz

            # Set basic triangle parameters, number, surface and vertices

            TheTriangle = TheMesh.elements[:,ElemZzz[index]]
            The2Surf = 2*TheMesh.element_measures[ElemZzz[index]]

            x1 = TheMesh.nodes[:, TheTriangle[1]]
            x2 = TheMesh.nodes[:, TheTriangle[2]]
            x3 = TheMesh.nodes[:, TheTriangle[3]]

            # Compute the inverse of the Jacobian matrix of the change of variables
            # Notice that the map from the real space to the reference element is 
            # given by [xi, eta] = InvJac*[x-x1,y-y1]

            InvJac = zeros(RealNum, (2,2))
            InvJac[1,1] = x3[2] - x1[2]
            InvJac[1,2] = x1[1] - x3[1]
            InvJac[2,1] = x1[2] - x2[2]
            InvJac[2,2] = x2[1] - x1[1]

            InvJac /= The2Surf  # The determinant of the jacobian is 2 times the surface of
                                # the triangle under consideration

            # Apply the change of variables

            XiEta = InvJac*(X-x1)

            # Compute the values of the 3 basis functions and store them in the corresponding
            # column of the "funvalues" matrix in BV

            BV.funvalues[:,index] = [1 - XiEta[1] - XiEta[2] , XiEta[1], XiEta[2]]
            
        end # of loop over relevant triangles

        BV.isempty = false
        return(nothing)

    end

end # of basis_func_pt_val_2D_P1 v2


# Developper version #1: modifies a basisfunvals_2D_P1 structure 
# The triangle of the mesh where the basis function shuld be evaluated is provided
# as an input parameter

function basis_func_pt_val_2D_P1_DEV!(X::RealVec, TheMesh::volcmesh, ElemZzz::Index,
                                    BV::RealVec)
    # This routine locates the 2D Point X in the mesh and
    # computes the values of the P1 basis functions in the triangles
    # where it belongs. It output a structure of 
    # Input variables:
    # - X : 2D point 
    # - TheMesh : mesh of class volcmesh
    # - ElemZzz : triangle of the mesh in which the point lies
    # - BV : here, just an array of three values to be filled
    #
    # Output: none
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    #
    # WARNING : no verification is made, this function is for developpers
    # of new assembly routines

    TheTriangle = zeros(Index, 3)

    # Set basic triangle parameters, number, surface and vertices

    TheTriangle = TheMesh.elements[:,ElemZzz]
    The2Surf = 2*TheMesh.element_measures[ElemZzz]

    x1 = TheMesh.nodes[:, TheTriangle[1]]
    x2 = TheMesh.nodes[:, TheTriangle[2]]
    x3 = TheMesh.nodes[:, TheTriangle[3]]

    # Compute the inverse of the Jacobian matrix of the change of variables
    # Notice that the map from the real space to the reference element is 
    # given by [xi, eta] = InvJac*[x-x1,y-y1]

    InvJac = zeros(RealNum, (2,2))
    InvJac[1,1] = x3[2] - x1[2]
    InvJac[1,2] = x1[1] - x3[1]
    InvJac[2,1] = x1[2] - x2[2]
    InvJac[2,2] = x2[1] - x1[1]

    InvJac /= The2Surf  # The determinant of the jacobian is 2 times the surface of
                        # the triangle under consideration

    # Apply the change of variables

    XiEta = InvJac*(X-x1)

    # Compute the values of the 3 basis functions and store them in the corresponding
    # column of the "funvalues" matrix in BV

    BV[1]= 1 - XiEta[1] - XiEta[2] 
    BV[2] = XiEta[1]
    BV[3] = XiEta[2]

    # println("               BV = ", BV)
            
    return(nothing)

end # of basis_func_pt_val_2D_P1 v2


### Gradients of P1 basis functions
#
# The gradients are constant over each element, only de change of variables
# has to be applied.

# First version : returning a tuple to be stored in a basisfunvals_2D_P1 structure

function basis_func_grad_pt_val_2D_P1(X::RealVec, TheMesh::volcmesh)
    # This routine locates the 2D Point X in the mesh and
    # computes the values of the P1 basis functions in the triangles
    # where it belongs. It output a structure of 
    # Input variables:
    # - X : 2D point 
    # - TheMesh : mesh of class volcmesh
    # - BV : basisfunvals_2D_P1 structure, empty or not
    # Output: none
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    
    TheTriangle = zeros(Index, 3)

    # Detect the elements where the point lies.
    # Calling the intriangle function will compute the element surfaces
    # and store them in the mesh if necessary.
    # No extra check is needed here.

    ElemZzz = intriangle(X,TheMesh)

    # println(ElemZzz)

    if (ElemZzz != 0) # if the point X lies in the mesh then compute the values

        NumElemZzz = length(ElemZzz)
        ValueZzz = zeros(RealNum,(2,3,NumElemZzz)) #Here we have to store NumelemZzz matrices

        
        GradMatrix = zeros(RealNum, (2,3))


        for index in 1:NumElemZzz

            # Set basic triangle parameters, number, surface and vertices

            TheTriangle = TheMesh.elements[:,ElemZzz[index]]
            The2Surf = 2*TheMesh.element_measures[ElemZzz[index]]

            x1 = TheMesh.nodes[:, TheTriangle[1]]
            x2 = TheMesh.nodes[:, TheTriangle[2]]
            x3 = TheMesh.nodes[:, TheTriangle[3]]

            # Compute the explicit gradient of the basis functions
            # WARNING : the gradients are stored columnwise

            # First basis function
            GradMatrix[1,1] = x2[2] - x3[2]
            GradMatrix[2,1] = x3[1] - x2[1]
            # Second basis function
            GradMatrix[1,2] = x3[2] - x1[2]
            GradMatrix[2,2] = x1[1] - x3[1]
            # Third basis function
            GradMatrix[1,3] = x1[2] - x2[2]
            GradMatrix[2,3] = x2[1] - x1[1]

            GradMatrix /= The2Surf  # The determinant of the jacobian is 2 times the surface of
                                # the triangle under consideration


            # Store them in the corresponding "column" of the result

            ValueZzz[:,:,index] = GradMatrix
            
        end # of loop over relevant triangles

        return(ElemZzz,ValueZzz,false)

    else # if the point is outside the mesh return empty arrays and the empty flag on
        return(@EmptyIntVec,@EmptyRealMat,true)
        
    end

end # of basis_func_grad_pt_val_2D_P1 v1 (IS v2 Useful ?)

function basis_func_stress_tensor_pt_val_2D_P1(X::RealVec, TheMesh::volcmesh, mu::Real, lambda::Real)
    # This routine locates the 2D Point X in the mesh and
    # computes the stress tensor for the P1 basis functions at that point.

    # Detect the elements where the point lies.
    ElemZzz = intriangle(X, TheMesh)

    if (ElemZzz != 0) # If the point X lies in the mesh
        NumElemZzz = length(ElemZzz)
        StressTensor = zeros(RealNum, (2, 2, 6, NumElemZzz)) # 2D stress tensor for each triangle
        epsilon_phi = zeros(RealNum, (2,2,6))
        trace_epsilon = zeros(RealNum, 6)
        

        for index in 1:NumElemZzz
            # Set basic triangle parameters
            TheTriangle = TheMesh.elements[:, ElemZzz[index]]
            The2Surf = 2*TheMesh.element_measures[ElemZzz[index]]

            x1 = TheMesh.nodes[:, TheTriangle[1]]
            x2 = TheMesh.nodes[:, TheTriangle[2]]
            x3 = TheMesh.nodes[:, TheTriangle[3]]


            # ICI BUG Corrigé dans les gradients : on a
            # PHI1 = [phi1; 0], PHI2 = [0, phi1] --> Noeud 1
            # PHI3 = [phi2; 0], PHI4 = [0, phi2] --> Noeud 2
            # PHI5 = [phi3; 0], PHI6 = [0, phi3] --> Noeud 3
            #
            # avec phi1 etc... fonctions de base scalires associées à chaque noeud
            #
            #
            # Donc la 
            # Compute the gradients of the basis functions (same as before)
            GradMatrix = zeros(RealNum, (2, 2, 6))
            # First node

            GradMatrix[1, 1, 1] = x2[2] - x3[2]
            GradMatrix[2, 1, 1] = 0
            GradMatrix[1, 2, 1] = x3[1] - x2[1]
            GradMatrix[2, 2, 1] = 0

            GradMatrix[1, 1, 2] = 0
            GradMatrix[2, 1, 2] = x2[2] - x3[2]
            GradMatrix[1, 2, 2] = 0
            GradMatrix[2, 2, 2] = x3[1] - x2[1]

            # Second node

            GradMatrix[1, 1, 3] = x3[2] - x1[2]
            GradMatrix[2, 1, 3] = 0
            GradMatrix[1, 2, 3] = x1[1] - x3[1]
            GradMatrix[2, 2, 3] = 0

            GradMatrix[2, 1, 4] = x3[2] - x1[2]
            GradMatrix[1, 1, 4] = 0
            GradMatrix[2, 2, 4] = x1[1] - x3[1]
            GradMatrix[1, 2, 4] = 0

            # Third node
           
            GradMatrix[2, 1, 5] = 0
            GradMatrix[1, 1, 5] = x1[2] - x2[2]
            GradMatrix[2, 2, 5] = 0
            GradMatrix[1, 2, 5] = x2[1] - x1[1]
            

            
            GradMatrix[1, 1, 6] = 0
            GradMatrix[2, 1, 6] = x1[2] - x2[2]
            GradMatrix[1, 2, 6] = 0
            GradMatrix[2, 2, 6] = x2[1] - x1[1]
            

            GradMatrix /= The2Surf # The determinant of the jacobian is 2 times the surface of
                                   # the triangle under consideration

            # Compute strain tensor epsilon(phi) as GradMatrix
            for i in 1:6
                epsilon_phi[:, :, i] = 0.5 * (GradMatrix[:, :, i] + GradMatrix[:, :, i]')  # Symmetric part
            end

            # Compute the trace of the strain tensor
            trace_epsilon[:] = epsilon_phi[1, 1, :] + epsilon_phi[2, 2, :]

            # Compute the stress tensor
            for i in 1:6
                StressTensor[:, :, i, index] = 2 * mu * epsilon_phi[:, :, i] + lambda * trace_epsilon[i] * I(2)
            end

        end # of loop over relevant triangles

        return(ElemZzz, StressTensor)

    else # If the point is outside the mesh
        return(@EmptyIntVec, @EmptyRealMat)
    end
end


### Stress tensor evaluation of basis functions for 2D elasticity
#
# The tensor being a 2x2 matrix we have to store 3 matrices per element in
# which the given point lies.
#
function basis_func_elast_tensor_2D_P1(X::RealVec, Young::RealNum , nu::RealNum, TheMesh::volcmesh)
    # This routine locates the 2D Point X in the mesh and
    # computes the values of the P1 basis functions in the triangles
    # where it belongs. It output a structure of 
    # Input variables:
    # - X : 2D point 
    # - TheMesh : mesh of class volcmesh
    # - Young, Nu : (scalar) values for the (constant) Young Modulus and Poisson Coef nu
    # Output: none
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #

    TheTriangle = zeros(Index, 3)

    # Detect the elements where the point lies.
    # Calling the intriangle function will compute the element surfaces
    # and store them in the mesh if necessary.
    # No extra check is needed here.

    ElemZzz = intriangle(X,TheMesh)

    if (ElemZzz != 0) # if the point X lies in the mesh then compute the values

        NumElemZzz = length(ElemZzz)
        ValueZzz = zeros(RealNum,(2,2,3,NumElemZzz)) # 3*NumElemZzz 2x2 matrices

        GradMatrix = zeros(RealNum, (2,3))

        lambda = (nu*Young) / ((1-nu)*(1+nu)) # Lamé coefficients
        mu = Young / (2*(1+nu))

        # FOR TEST PURPOSES -- REMOVE LATER
        # lambda = 1.0
        # mu = 1.0


        for index in 1:NumElemZzz # Loop over elements where the point lies


            #### First compute the gradients of the shape functions


            # Set basic triangle parameters, number, surface and vertices

            TheTriangle = TheMesh.elements[:,ElemZzz[index]]
            The2Surf = 2*TheMesh.element_measures[ElemZzz[index]]

            x1 = TheMesh.nodes[:, TheTriangle[1]]
            x2 = TheMesh.nodes[:, TheTriangle[2]]
            x3 = TheMesh.nodes[:, TheTriangle[3]]

            # Compute the explicit gradient of the basis functions
            # WARNING : the gradients are stored columnwise

            # First basis function
            GradMatrix[1,1] = x2[2] - x3[2]
            GradMatrix[2,1] = x3[1] - x2[1]
            # Second basis function
            GradMatrix[1,2] = x3[2] - x1[2]
            GradMatrix[2,2] = x1[1] - x3[1]
            # Third basis function
            GradMatrix[1,3] = x1[2] - x2[2]
            GradMatrix[2,3] = x2[1] - x1[1]

            GradMatrix /= The2Surf  # The determinant of the jacobian is 2 times the surface of
                                # the triangle under consideration

            # Warning : GradMatrix is for scalar fields. We have to take the vector structure
            # into account.

            #### Second, for each node, compute the tensor :
            # build two times the symmetric part of the gradient (epsilon) and the identity
            # matrix times the trace of the symmetric part of the gradient.

            # Node 1 (column 1 of GradMatrix) #
            TwoTimesEpsilon = [2.0*GradMatrix[1,1] (GradMatrix[1,1]+GradMatrix[2,1]);
                        (GradMatrix[1,1]+GradMatrix[2,1]) 2.0*GradMatrix[2,1]]
            TraceEpsId = [(GradMatrix[1,1]+GradMatrix[2,1]) 0.0 ;
                            0.0 (GradMatrix[1,1]+GradMatrix[2,1])]
            sigma1 = mu*TwoTimesEpsilon + lambda*TraceEpsId

            # Node 2 (column 2 of GradMatrix) #
            TwoTimesEpsilon = [2.0*GradMatrix[1,2] (GradMatrix[1,2]+GradMatrix[2,2]);
                        (GradMatrix[1,2]+GradMatrix[2,2]) 2.0*GradMatrix[2,2]]
            TraceEpsId = [(GradMatrix[1,2]+GradMatrix[2,2]) 0.0 ;
                            0.0 (GradMatrix[1,2]+GradMatrix[2,2])]
            sigma2 = mu*TwoTimesEpsilon + lambda*TraceEpsId

            # Node 3 (column 3 of GradMatrix) #
            TwoTimesEpsilon = [2.0*GradMatrix[1,3] (GradMatrix[1,3]+GradMatrix[2,3]);
                        (GradMatrix[1,3]+GradMatrix[2,3]) 2.0*GradMatrix[2,3]]
            TraceEpsId = [(GradMatrix[1,3]+GradMatrix[2,3]) 0.0 ;
                            0.0 (GradMatrix[1,3]+GradMatrix[2,3])]
            sigma3 = mu*TwoTimesEpsilon + lambda*TraceEpsId

            # Store them in the corresponding "column" of the result

            ValueZzz[:,:,1,index] = sigma1
            ValueZzz[:,:,2,index] = sigma2
            ValueZzz[:,:,3,index] = sigma3
            
            
        end # of loop over relevant triangles

        return(ElemZzz,ValueZzz,false)

    else # if the point is outside the mesh return empty arrays and the empty flag on
        return(@EmptyIntVec,@EmptyRealMat,true)
        
    end

end # of basis_func_elast_tensor_2D_P1 (IS v2 Useful ?)


###########################################################################
#################                                       ###################
#################       Mass and stiffness matrices     ###################
#################                                       ###################
###########################################################################

########################
### Basic mass matrices
########################


# Basic mass matrix assembly (2D, P1 elements)
#
# M = [integral(Phi_i(x)Phi_j(x)dx)]_ij
#
# where {Phi_i} are the P1 shape functions on the mesh
#
# On each triangle of the mesh, the integrals are explicitely known :
# m_ii = |T| / 6, m_ij = |T| / 12
# where |T| is the surface of the triangle
function asm_mass_2D_P1(TheMesh::volcmesh)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - Mass : the mass matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    M = spzeros(nbdofs,nbdofs)
    #
    nbtriangles = nbelements(TheMesh)
    TheTriangle = zeros(Index, 3)
    #
    # Check if the element measures are present
    # if not, compute them and store them into the mesh object
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    #
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        #
        # Compute the surface of the triangle
        #
        TheMeasureDiv6 = TheMesh.element_measures[index] / 6.0  # This is |T|/6
        TheMeasureDiv12 = 0.5*TheMeasureDiv6    # This is |T|/12
        #
        # Assemble the elementary mass terms
        #
        M[TheTriangle[1],TheTriangle[1]] += TheMeasureDiv6
        M[TheTriangle[2],TheTriangle[2]] += TheMeasureDiv6
        M[TheTriangle[3],TheTriangle[3]] += TheMeasureDiv6
        M[TheTriangle[1],TheTriangle[2]] += TheMeasureDiv12
        M[TheTriangle[1],TheTriangle[3]] += TheMeasureDiv12
        M[TheTriangle[2],TheTriangle[1]] += TheMeasureDiv12
        M[TheTriangle[2],TheTriangle[3]] += TheMeasureDiv12
        M[TheTriangle[3],TheTriangle[1]] += TheMeasureDiv12
        M[TheTriangle[3],TheTriangle[2]] += TheMeasureDiv12
        #
    end

    return(M)
end # of assemb_mass_2D_P1

# Vector mass matrix assembly 
# i.e. mass matrix assembly with 2 dofs per node.
#
# M = [integral(Phi_i(x)Phi_j(x)dx)]_ij
#
# where {Phi_i} are the P1 shape functions on the mesh
#
# On each triangle of the mesh, the integrals are explicitely known :
# m_ii = |T| / 6, m_ij = |T| / 12
# where |T| is the surface of the triangle

function asm_mass_vect_2D_P1(TheMesh::volcmesh)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # with two dofs per node
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - Mass : the mass matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs::Index = 2*nbnodes(TheMesh) # We have 2 dofs at each node, one for each component of the vector field

    M = spzeros(nbdofs,nbdofs)
    #
    nbtriangles::Index = nbelements(TheMesh)
    TheTriangle = zeros(Index, (3,1))
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    #
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        Node1 = 2*TheTriangle[1] - 1
        Node2 = 2*TheTriangle[2] - 1
        Node3 = 2*TheTriangle[3] - 1
        Node1bis = 2*TheTriangle[1]
        Node2bis = 2*TheTriangle[2]
        Node3bis = 2*TheTriangle[3]
        #
        # Compute the surface of the triangle
        #
        TheMeasureDiv6 = TheMesh.element_measures[index] / 6.  # This is |T|/12
        TheMeasureDiv12 = 0.5*TheMeasureDiv6    # This is |T|/12
        #
        # Assemble the elementary mass terms
        # Here Node1 has 2 ddls : 2*Node1 and 2*Node1 + 1 etc...
        #
        # "Diagonal terms"
        M[Node1,Node1] += TheMeasureDiv6
        M[Node2,Node2] += TheMeasureDiv6
        M[Node3,Node3] += TheMeasureDiv6
        M[Node1bis,Node1bis] += TheMeasureDiv6
        M[Node2bis,Node2bis] += TheMeasureDiv6
        M[Node3bis,Node3bis] += TheMeasureDiv6
        # "Non-diagonal terms"
        M[Node1,Node2] += TheMeasureDiv12
        M[Node1,Node3] += TheMeasureDiv12
        M[Node2,Node1] += TheMeasureDiv12
        M[Node2,Node3] += TheMeasureDiv12
        M[Node3,Node1] += TheMeasureDiv12
        M[Node3,Node2] += TheMeasureDiv12
        M[Node1bis,Node2bis] += TheMeasureDiv12
        M[Node1bis,Node3bis] += TheMeasureDiv12
        M[Node2bis,Node1bis] += TheMeasureDiv12
        M[Node2bis,Node3bis] += TheMeasureDiv12
        M[Node3bis,Node1bis] += TheMeasureDiv12
        M[Node3bis,Node2bis] += TheMeasureDiv12
        #
    end # for index
    
    return(M)
end # of assemb_mass_vect


# Boundary mass matrix assembly
# This function can be used to penalize dirichlet boundary conditions
# or impose Neumann boundary conditions.

function asm_bnd_mass_simpson_2D_P1(TheMesh::volcmesh, BoundaryNumber::Index = 0)
    # Mass matrix approximation via the trapeze formula
    # on each segment of the boundary.
    # Input:
    # - TheMesh : a volcmesh object
    # - BoundaryNumber : number of the boundary to treat. If 0, all the boundary will be processed.
    # Output:
    # - Mb : boundary mass matrix
    #
    # The resulting matrix has the same dimension as the number of unknows but containes 0 for non boundary nodes.
    # First read the number of dofs and declare the matrix
    
    #
    nbdofs = nbnodes(TheMesh)
    Mb = spzeros(nbdofs,nbdofs)
    
    #
    # Check if the edge measures are present
    #
    if !TheMesh.edg_meas_flag
        mk_edge_measures(TheMesh)
    end
    #
    #
    #
    # Extract the indexes of the segments to scan
    # 
    if (BoundaryNumber != 0)
        #
        # Detect the edges labeled as part of BoundaryNumber
        #
        BndSegments = findall(isequal(BoundaryNumber) , TheMesh.edges[3,:])
    else
        #
        # In the case when BoundaryNumber is 0, we treat the whole boundary
        # 
        BndSegments = 1:nbedges(TheMesh)
    end # If Else BoundaryNumber != 0
    #
    # Init. various variables
    #
    TheEdge = zeros(Index, (3,1))   # this has nothing to do with U2 of course.
    #
    # Loop over edges and build the elementary matrix
    #
    for numedge in BndSegments
        #
        TheEdge = TheMesh.edges[:,numedge]
        Begin = TheEdge[1]
        End = TheEdge[2]

        x_beg = TheMesh.nodes[:, Begin]
        x_end = TheMesh.nodes[:, End]
        x_mid = (x_beg + x_end)/2

        _,basis_func_val_mid,_ = basis_func_pt_val_2D_P1(x_mid, TheMesh)

        #
        # Compute the length of the segment (no function call to gain speed) 
        #
        # HalfEdgeLength = 0.5 * TheMesh.edge_measures[numedge]
        length = (TheMesh.edge_measures[numedge] / 6.) 

        # Find the indices of non-zero and non-null values
        indices = findall(x -> !isnan(x) && x != 0, basis_func_val_mid)

        # Assuming you know there are exactly 2 non-null values, store them in variables
        i1, i2 = indices[1], indices[2]
        #
        # Assemble the elementary terms (the elem matrix is diagonal)
        #
        # Mb[Begin,Begin] += HalfEdgeLength
        # Mb[End,End] += HalfEdgeLength

        # Comment avoir la numerotations de la bonne fonction de base. 
        # est ce que c'est lie avec le numerio de l'edge
        Mb[Begin,Begin] += length*(1 + 4*basis_func_val_mid[i1]*basis_func_val_mid[i1])
        Mb[End,End] += length*(1 + 4*basis_func_val_mid[i2]*basis_func_val_mid[i2])
        Mb[Begin, End] += length*(4*basis_func_val_mid[i1]*basis_func_val_mid[i2])
        Mb[End, Begin] += length*(4*basis_func_val_mid[i1]*basis_func_val_mid[i2])
        #
    end # for numedge

    #
    return(Mb)
end # of assemb_boundary_mass


function asm_bnd_mass_simpson_2D_P1_vect(TheMesh::volcmesh, BoundaryNumber::Index = 0)
    # Mass matrix approximation via the trapeze formula
    # on each segment of the boundary.
    # Input:
    # - TheMesh : a volcmesh object
    # - BoundaryNumber : number of the boundary to treat. If 0, all the boundary will be processed.
    # Output:
    # - Mb : boundary mass matrix
    #
    # The resulting matrix has the same dimension as the number of unknows but containes 0 for non boundary nodes.
    # First read the number of dofs and declare the matrix
    
    #
    nbdofs = 2*nbnodes(TheMesh)
    Mb = spzeros(nbdofs,nbdofs)
    
    #
    # Check if the edge measures are present
    #
    if !TheMesh.edg_meas_flag
        mk_edge_measures(TheMesh)
    end
    #
    #
    #
    # Extract the indexes of the segments to scan
    # 
    if (BoundaryNumber != 0)
        #
        # Detect the edges labeled as part of BoundaryNumber
        #
        BndSegments = findall(isequal(BoundaryNumber) , TheMesh.edges[3,:])
    else
        #
        # In the case when BoundaryNumber is 0, we treat the whole boundary
        # 
        BndSegments = 1:nbedges(TheMesh)
    end # If Else BoundaryNumber != 0
    #
    # Init. various variables
    #
    TheEdge = zeros(Index, (3,1))   # this has nothing to do with U2 of course.
    #
    # Loop over edges and build the elementary matrix
    #
    for numedge in BndSegments
        #
        TheEdge = TheMesh.edges[:,numedge]
        Begin = TheEdge[1]
        End = TheEdge[2]

        x_beg = TheMesh.nodes[:, Begin]
        x_end = TheMesh.nodes[:, End]
        x_mid = (x_beg + x_end)/2


        _,basis_func_val_mid,_ = basis_func_pt_val_2D_P1(x_mid, TheMesh)

        #
        # Compute the length of the segment (no function call to gain speed) 
        #
        # HalfEdgeLength = 0.5 * TheMesh.edge_measures[numedge]
        length = (TheMesh.edge_measures[numedge] / 6.) 

        # Find the indices of non-zero and non-null values
        indices = findall(x -> !isnan(x) && x != 0, basis_func_val_mid)

        # Assuming you know there are exactly 2 non-null values, store them in variables
        i1, i2 = indices[1], indices[2]
        #
        # Assemble the elementary terms (the elem matrix is diagonal)
        #
        # Mb[Begin,Begin] += HalfEdgeLength
        # Mb[End,End] += HalfEdgeLength

        # Comment avoir la numerotations de la bonne fonction de base. 
        # est ce que c'est lie avec le numerio de l'edge
        Mb[2 * Begin, 2 * Begin] += length*(1 + 4*basis_func_val_mid[i1]*basis_func_val_mid[i1])
        Mb[2 * End, 2 * End] += length*(1 + 4*basis_func_val_mid[i2]*basis_func_val_mid[i2])
        Mb[2 * Begin, 2 * End] += length*(4*basis_func_val_mid[i1]*basis_func_val_mid[i2])
        Mb[2 * End, 2 * Begin] += length*(4*basis_func_val_mid[i1]*basis_func_val_mid[i2])
        
        Mb[2 * Begin - 1, 2 * Begin -1] += length*(1 + 4*basis_func_val_mid[i1]*basis_func_val_mid[i1])
        Mb[2 * End - 1, 2 * End - 1] += length*(1 + 4*basis_func_val_mid[i2]*basis_func_val_mid[i2])
        Mb[2 * Begin - 1, 2 * End -1] += length*(4*basis_func_val_mid[i1]*basis_func_val_mid[i2])
        Mb[2 * End - 1, 2 * Begin - 1] += length*(4*basis_func_val_mid[i1]*basis_func_val_mid[i2])
        #
    end # for numedge

    #
    return(Mb)
end # of assemb_boundary_mass

function asm_bnd_mass_2D_P1(TheMesh::volcmesh, BoundaryNumber::Index = 0)
    # Mass matrix approximation via the trapeze formula
    # on each segment of the boundary.
    # Input:
    # - TheMesh : a volcmesh object
    # - BoundaryNumber : number of the boundary to treat. If 0, all the boundary will be processed.
    # Output:
    # - Mb : boundary mass matrix
    #
    # The resulting matrix has the same dimension as the number of unknows but containes 0 for non boundary nodes.
    # First read the number of dofs and declare the matrix
    
    #
    nbdofs = nbnodes(TheMesh)
    Mb = spzeros(nbdofs,nbdofs)
    
    #
    # Check if the edge measures are present
    #
    if !TheMesh.edg_meas_flag
        mk_edge_measures(TheMesh)
    end
    #
    #
    #
    # Extract the indexes of the segments to scan
    # 
    if (BoundaryNumber != 0)
        #
        # Detect the edges labeled as part of BoundaryNumber
        #
        BndSegments = findall(isequal(BoundaryNumber) , TheMesh.edges[3,:])
    else
        #
        # In the case when BoundaryNumber is 0, we treat the whole boundary
        # 
        BndSegments = 1:nbedges(TheMesh)
    end # If Else BoundaryNumber != 0
    #
    # Init. various variables
    #
    TheEdge = zeros(Index, (3,1))   # this has nothing to do with U2 of course.
    #
    # Loop over edges and build the elementary matrix
    #
    for numedge in BndSegments
        #
        TheEdge = TheMesh.edges[:,numedge]
        Begin = TheEdge[1]
        End = TheEdge[2]

        #
        # Compute the length of the segment (no function call to gain speed) 
        #
        HalfEdgeLength = 0.5 * TheMesh.edge_measures[numedge]
        #
        # Assemble the elementary terms (the elem matrix is diagonal)
        #
        Mb[Begin,Begin] += HalfEdgeLength
        Mb[End,End] += HalfEdgeLength
        #
    end # for numedge

    #
    return(Mb)
end # of assemb_boundary_mass

# Boundary mass matrix assembly for 2 DOFs per node
# This function can be used to penalize Dirichlet boundary conditions
# or impose Neumann boundary conditions.

function asm_bnd_mass_vect_2D_P1(TheMesh::volcmesh, BoundaryNumber::Index = 0)
    # Mass matrix approximation via the trapeze formula
    # on each segment of the boundary.
    # Input:
    # - TheMesh : a volcmesh object
    # - BoundaryNumber : number of the boundary to treat. If 0, all the boundary will be processed.
    # Output:
    # - Mb : boundary mass matrix (with 2 DOFs per node)
    
    # First read the number of nodes and declare the matrix
    nbdofs = 2 * nbnodes(TheMesh)  # 2 DOFs per node
    Mb = spzeros(nbdofs, nbdofs)

    # Check if the edge measures are present
    if !TheMesh.edg_meas_flag
        mk_edge_measures(TheMesh)
    end

    # Extract the indexes of the segments to scan
    if (BoundaryNumber != 0)
        # Detect the edges labeled as part of BoundaryNumber
        BndSegments = findall(isequal(BoundaryNumber), TheMesh.edges[3, :])
    else
        # In the case when BoundaryNumber is 0, we treat the whole boundary
        BndSegments = 1:nbedges(TheMesh)
    end

    # Loop over edges and build the elementary matrix
    for numedge in BndSegments
        TheEdge = TheMesh.edges[:, numedge]
        Begin = TheEdge[1]
        End = TheEdge[2]
        
        # Compute the length of the segment (no function call to gain speed)
        HalfEdgeLength = 0.5 * TheMesh.edge_measures[numedge]
        
        # Assemble the elementary terms (the elem matrix is diagonal for 2 DOFs)
        Mb[2 * Begin - 1, 2 * Begin - 1] += HalfEdgeLength    # DOF 1 for Begin
        Mb[2 * Begin, 2 * Begin] += HalfEdgeLength            # DOF 2 for Begin
        Mb[2 * End - 1, 2 * End - 1] += HalfEdgeLength        # DOF 1 for End
        Mb[2 * End, 2 * End] += HalfEdgeLength                # DOF 2 for End
    end

    return Mb
end

function asm_immersed_bnd_mass_2D_P1(TheMesh::volcmesh,  TheLine::linefrac2D)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - Mass : the mass matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    M = spzeros(nbdofs,nbdofs)

    Nf = nbelements(TheLine) 
    #
    # Compute the intersection points of the mesh and the segments
    (HasIntersection, InterSectPoints, Seglengths) = frac_mesh_intersect(TheLine)

    # Compute the segment normals if necessary

    if !TheLine.elm_norm_flag
        normals = mk_segment_normals(TheLine) # normals is a boolean
    end
       

    # Loop over the segments of TheLine
    for jndex = 1:Nf
        CurrentElem = TheLine.elements[:, jndex]
        

        # Get the start and end points of the current segment
        p1 = TheLine.nodes[:, CurrentElem[1]] # First point
        p2 = TheLine.nodes[:, CurrentElem[2]] # Last point

        # Get the triangles of the mesh containing each point
        Tr1 = TheLine.embed_triangles[CurrentElem[1]]
        Tr2 = TheLine.embed_triangles[CurrentElem[2]]

        if (HasIntersection[jndex]) # If the segment intersects two triangles
            Y = InterSectPoints[:, jndex]
            lg = Seglengths[:, jndex]

            ############
            # Triangle 1
            ############

            TheTriangle = TheMesh.elements[:,Tr1]

            _,phi_i,_ = basis_func_pt_val_2D_P1(p1, TheMesh)
            _,phi_j,_ = basis_func_pt_val_2D_P1(Y, TheMesh)

    
            for i in 1:3
                for j in 1:3
                    # println(TheTriangle[i], "  " ,TheTriangle[j], "  ", phi_i[i], "   ", phi_j[j, Tr1])
                    M[TheTriangle[i], TheTriangle[j]] += (lg[1]/2.0) * (phi_i[i] * phi_i[j] + phi_j[i,1] * phi_j[j, 1])
                end
            end

            ############
            # Triangle 2
            ############

            TheTriangle = TheMesh.elements[:,Tr2]

            _,phi_i,_ = basis_func_pt_val_2D_P1(Y, TheMesh)
            _,phi_j,_ = basis_func_pt_val_2D_P1(p2, TheMesh)


            for i in 1:3
                for j in 1:3
                    M[TheTriangle[i], TheTriangle[j]] += (lg[2]/2.0) * (phi_j[i] * phi_j[j] + phi_i[i,2] * phi_i[j,2])
                end
            end

        else # Segment lies entirely in one triangle

            TheTriangle = TheMesh.elements[:, Tr1]

            _,phi_i,_ = basis_func_pt_val_2D_P1(p1, TheMesh)
            _,phi_j,_ = basis_func_pt_val_2D_P1(p2, TheMesh)

             for i in 1:3
                for j in 1:3
                    M[TheTriangle[i], TheTriangle[j]] += (TheLine.element_measures[jndex]/2.0) * (phi_i[i] * phi_i[j] + phi_j[i] * phi_j[j])
                end
            end

         end
    end

    return M
end # of assemb_mass_2D_P1


function asm_immersed_bnd_mass_2D_P1_2DOF(TheMesh::volcmesh,  TheLine::linefrac2D)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - Mass : the mass matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    M = spzeros(2*nbdofs,2*nbdofs)

    Nf = nbelements(TheLine) 
    #
    # Compute the intersection points of the mesh and the segments
    (HasIntersection, InterSectPoints, Seglengths) = frac_mesh_intersect(TheLine)

    # Compute the segment normals if necessary

    if !TheLine.elm_norm_flag
        normals = mk_segment_normals(TheLine) # normals is a boolean
    end
       

    # Loop over the segments of TheLine
    for jndex = 1:Nf
        CurrentElem = TheLine.elements[:, jndex]
        

        # Get the start and end points of the current segment
        p1 = TheLine.nodes[:, CurrentElem[1]] # First point
        p2 = TheLine.nodes[:, CurrentElem[2]] # Last point

        # Get the triangles of the mesh containing each point
        Tr1 = TheLine.embed_triangles[CurrentElem[1]]
        Tr2 = TheLine.embed_triangles[CurrentElem[2]]

        if (HasIntersection[jndex]) # If the segment intersects two triangles
            Y = InterSectPoints[:, jndex]
            lg = Seglengths[:, jndex]

            ############
            # Triangle 1
            ############

            TheTriangle = TheMesh.elements[:,Tr1]

            _,phi_i,_ = basis_func_pt_val_2D_P1(p1, TheMesh)
            _,phi_j,_ = basis_func_pt_val_2D_P1(Y, TheMesh)

    
            for i in 1:3
                for j in 1:3
                    M[2 * TheTriangle[i], 2 * TheTriangle[j]] += (lg[1]/2.0) * (phi_i[i] * phi_i[j] + phi_j[i,1] * phi_j[j, 1])
                    M[2 * TheTriangle[i] - 1, 2 * TheTriangle[j] - 1] += (lg[1]/2.0) * (phi_i[i] * phi_i[j] + phi_j[i,1] * phi_j[j, 1])
                
                end
            end

            ############
            # Triangle 2
            ############

            TheTriangle = TheMesh.elements[:,Tr2]

            _,phi_i,_ = basis_func_pt_val_2D_P1(Y, TheMesh)
            _,phi_j,_ = basis_func_pt_val_2D_P1(p2, TheMesh)


            for i in 1:3
                for j in 1:3
                    M[2 * TheTriangle[i], 2 * TheTriangle[j]] += (lg[2]/2.0) * (phi_j[i] * phi_j[j] + phi_i[i,2] * phi_i[j,2])
                    M[2 * TheTriangle[i] - 1, 2 * TheTriangle[j] - 1] += (lg[2]/2.0) * (phi_j[i] * phi_j[j] + phi_i[i,2] * phi_i[j,2])
                end
            end

        else # Segment lies entirely in one triangle

            TheTriangle = TheMesh.elements[:, Tr1]

            _,phi_i,_ = basis_func_pt_val_2D_P1(p1, TheMesh)
            _,phi_j,_ = basis_func_pt_val_2D_P1(p2, TheMesh)

             for i in 1:3
                for j in 1:3
                    M[2 * TheTriangle[i], 2 * TheTriangle[j]] += (TheLine.element_measures[jndex]/2.0) * (phi_i[i] * phi_i[j] + phi_j[i] * phi_j[j])
                    M[2 * TheTriangle[i] - 1, 2 * TheTriangle[j] - 1] += (TheLine.element_measures[jndex]/2.0) * (phi_i[i] * phi_i[j] + phi_j[i] * phi_j[j])
                end
            end

         end
    end

    return M
end # of assemb_mass_2D_P1


###############################################################
# Weighted mass matrices : a weight vector or function is given
###############################################################


function asm_mass_2D_P1(TheMesh::volcmesh, W::RealVec)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # - W : a vector with values of a function defined on the mesh nodes
    # Output:
    # - Mass : the mass matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    M = spzeros(nbdofs,nbdofs)
    #
    nbtriangles = nbelements(TheMesh)
    TheTriangle = zeros(Index, 3)
    #
    # Check if the element measures are present
    # if not, compute them and store them into the mesh object
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    #
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the node numbers of the element
        #
        TheTriangle = TheMesh.elements[:,index]
        #        
        # Interpolate the weight at the center of the triangle
        #
        Wbar = (W[TheTriangle[1]] + W[TheTriangle[2]] + W[TheTriangle[3]]) / 3.0
        #
        # Compute the surface of the triangle
        #
        TheMeasureDiv6 = (Wbar*TheMesh.element_measures[index]) / 6.0  # This is Wbar*|T|/6
        TheMeasureDiv12 = 0.5*TheMeasureDiv6    # This is Wbar*|T|/12
        #
        # Assemble the elementary mass terms
        #
        M[TheTriangle[1],TheTriangle[1]] += TheMeasureDiv6
        M[TheTriangle[2],TheTriangle[2]] += TheMeasureDiv6
        M[TheTriangle[3],TheTriangle[3]] += TheMeasureDiv6
        M[TheTriangle[1],TheTriangle[2]] += TheMeasureDiv12
        M[TheTriangle[1],TheTriangle[3]] += TheMeasureDiv12
        M[TheTriangle[2],TheTriangle[1]] += TheMeasureDiv12
        M[TheTriangle[2],TheTriangle[3]] += TheMeasureDiv12
        M[TheTriangle[3],TheTriangle[1]] += TheMeasureDiv12
        M[TheTriangle[3],TheTriangle[2]] += TheMeasureDiv12
        #
    end

    return(M)
end # of assemb_mass_2D_P1


function asm_mass_2D_P1(TheMesh::volcmesh, fonk::Function)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # - F : a function with values computable on the mesh nodes
    # Output:
    # - M : the mass matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    M = spzeros(nbdofs,nbdofs)
    #
    nbtriangles = nbelements(TheMesh)
    TheTriangle = zeros(Index, 3)
    #
    # Check if the element measures are present
    # if not, compute them and store them into the mesh object
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    #
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the node numbers of the element
        #
        TheTriangle = TheMesh.elements[:,index]
        #        
        # Interpolate the weight at the center of the triangle
        #
        W1 = fonk(TheMesh.nodes[:,TheTriangle[1]])
        W2 = fonk(TheMesh.nodes[:,TheTriangle[2]])
        W3 = fonk(TheMesh.nodes[:,TheTriangle[3]])
        Wbar = (W1 + W2 + W3) / 3.0
        #
        # Compute the surface of the triangle
        #
        TheMeasureDiv6 = (Wbar*TheMesh.element_measures[index]) / 6.0  # This is Wbar*|T|/6
        TheMeasureDiv12 = 0.5*TheMeasureDiv6    # This is Wbar*|T|/12
        #
        # Assemble the elementary mass terms
        #
        M[TheTriangle[1],TheTriangle[1]] += TheMeasureDiv6
        M[TheTriangle[2],TheTriangle[2]] += TheMeasureDiv6
        M[TheTriangle[3],TheTriangle[3]] += TheMeasureDiv6
        M[TheTriangle[1],TheTriangle[2]] += TheMeasureDiv12
        M[TheTriangle[1],TheTriangle[3]] += TheMeasureDiv12
        M[TheTriangle[2],TheTriangle[1]] += TheMeasureDiv12
        M[TheTriangle[2],TheTriangle[3]] += TheMeasureDiv12
        M[TheTriangle[3],TheTriangle[1]] += TheMeasureDiv12
        M[TheTriangle[3],TheTriangle[2]] += TheMeasureDiv12
        #
    end

    return(M)
end # of assemb_mass_2D_P1

function asm_mass_2D_P1_2DOF(TheMesh::volcmesh)
    # This routine builds the mass matrix for P1 elements with two degrees of freedom on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - M : the mass matrix on the mesh for P1 elements with two degrees of freedom
    
    # First read the number of dofs and declare the matrix
    nbdofs = 2 * nbnodes(TheMesh)  # Two degrees of freedom per node
    M = spzeros(nbdofs, nbdofs)
    
    # Initialize various variables
    nbtriangles::Index = nbelements(TheMesh)
    TheTriangle = zeros(Int, 3)
    
    # Check if the element measures are present
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    
    # Loop over triangles and build the elementary matrix
    for index in 1:nbtriangles
        # Get the node numbers of the element
        TheTriangle = TheMesh.elements[:, index]
        
        # Compute the surface of the triangle
        TheMeasureDiv6 = TheMesh.element_measures[index] / 6.0  # This is |T|/6
        TheMeasureDiv12 = 0.5 * TheMeasureDiv6                   # This is |T|/12
        
        # Assemble the elementary mass terms for two degrees of freedom
        for i in 1:3
            for j in 1:3
                if i == j
                    M[2 * TheTriangle[i] - 1, 2 * TheTriangle[j] - 1] += TheMeasureDiv6
                    M[2 * TheTriangle[i], 2 * TheTriangle[j]] += TheMeasureDiv6
                else
                    M[2 * TheTriangle[i] - 1, 2 * TheTriangle[j] - 1] += TheMeasureDiv12
                    M[2 * TheTriangle[i], 2 * TheTriangle[j]] += TheMeasureDiv12
                    M[2 * TheTriangle[i] - 1, 2 * TheTriangle[j]] += TheMeasureDiv12
                    M[2 * TheTriangle[i], 2 * TheTriangle[j] - 1] += TheMeasureDiv12
                end
            end
        end
    end

    return M
end  # of asm_mass_2D_P1_2dof



function asm_bnd_mass_2D_P1(TheMesh::volcmesh, W::RealVec, BoundaryNumber::Index = 0)
    # Mass matrix approximation via the trapeze formula
    # on each segment of the boundary.
    # Input:
    # - TheMesh : a volcmesh object
    # - W : a vector with values of a function defined on the mesh nodes
    # - BoundaryNumber : number of the boundary to treat. If 0, all the boundary will be processed.
    # Output:
    # - Mb : boundary mass matrix
    #
    # The resulting matrix has the same dimension as the number of unknows but containes 0 for non boundary nodes.
    # First read the number of dofs and declare the matrix
    
    #
    nbdofs = nbnodes(TheMesh)
    Mb = spzeros(nbdofs,nbdofs)
    
    #
    # Check if the edge measures are present
    #
    if !TheMesh.edg_meas_flag
        mk_edge_measures(TheMesh)
    end
    #
    #
    #
    # Extract the indexes of the segments to scan
    # 
    if (BoundaryNumber != 0)
        #
        # Detect the edges labeled as part of BoundaryNumber
        #
        BndSegments = findall(isequal(BoundaryNumber) , TheMesh.edges[3,:])
    else
        #
        # In the case when BoundaryNumber is 0, we treat the whole boundary
        # 
        BndSegments = 1:nbedges(TheMesh)
    end # If Else BoundaryNumber != 0
    #
    # Init. various variables
    #
    TheEdge = zeros(Index, (3,1))   # this has nothing to do with U2 of course.
    #
    # Loop over edges and build the elementary matrix
    #
    for numedge in BndSegments
        #
        TheEdge = TheMesh.edges[:,numedge]
        Begin = TheEdge[1]
        End = TheEdge[2]

    
        #
        # Compute the length of the segment (no function call to gain speed) 
        #
        HalfEdgeLength = 0.5 * TheMesh.edge_measures[numedge]
        #
        # Assemble the elementary terms (the elem matrix is diagonal)
        #
        Mb[Begin,Begin] += HalfEdgeLength*W[Begin]
        Mb[End,End] += HalfEdgeLength*W[End]
        #
    end # for numedge
    #
    return(Mb)
end # of assemb_boundary_mass_weight

function asm_bnd_mass_2D_P1(TheMesh::volcmesh, fonk::Function, BoundaryNumber::Index = 0)
    # Mass matrix approximation via the trapeze formula
    # on each segment of the boundary.
    # Input:
    # - TheMesh : a volcmesh object
    # - fonk : a function taking values on the mesh nodes
    # - BoundaryNumber : number of the boundary to treat. If 0, all the boundary will be processed.
    # Output:
    # - Mb : boundary mass matrix
    #
    # The resulting matrix has the same dimension as the number of unknows but containes 0 for non boundary nodes.
    # First read the number of dofs and declare the matrix
    
    #
    nbdofs = nbnodes(TheMesh)
    Mb = spzeros(nbdofs,nbdofs)
    
    #
    # Check if the edge measures are present
    #
    if !TheMesh.edg_meas_flag
        mk_edge_measures(TheMesh)
    end
    #
    #
    #
    # Extract the indexes of the segments to scan
    # 
    if (BoundaryNumber != 0)
        #
        # Detect the edges labeled as part of BoundaryNumber
        #
        BndSegments = findall(isequal(BoundaryNumber) , TheMesh.edges[3,:])
    else
        #
        # In the case when BoundaryNumber is 0, we treat the whole boundary
        # 
        BndSegments = 1:nbedges(TheMesh)
    end # If Else BoundaryNumber != 0
    #
    # Init. various variables
    #
    TheEdge = zeros(Index, (3,1))   # this has nothing to do with U2 of course.
    #
    # Loop over edges and build the elementary matrix
    #
    for numedge in BndSegments
        #
        TheEdge = TheMesh.edges[:,numedge]
        Begin = TheEdge[1]
        End = TheEdge[2]

        println(Begin, End)
        #
        # Compute the length of the segment (no function call to gain speed) 
        #
        HalfEdgeLength = 0.5 * TheMesh.edge_measures[numedge]
        #
        # Compute the values of the weight function at start and end nodes of the segment
        #
        Wbegin = fonk(TheMesh.nodes[:,Begin])
        Wend = fonk(TheMesh.nodes[:,End])
        #
        # Assemble the elementary terms (the elem matrix is diagonal)
        #
        Mb[Begin,Begin] += HalfEdgeLength*Wbegin
        Mb[End,End] += HalfEdgeLength*Wend
        #
    end # for numedge
    #
    return(Mb)
end # of assemb_boundary_mass_weight


function asm_bnd_mass_matrix_on_crack_2DOF_P0(TheLine::linefrac2D)::RealMat
    # Get the number of elements in the line
    Ne = nbelements(TheLine)
    
    # Initialize the mass matrix (size Ne x Ne)
    mass_matrix = zeros(RealNum, 2*Ne, 2*Ne)

    # Loop over each element to populate the mass matrix
    for i in 1:Ne
        # Length of the current segment
        length_i = TheLine.element_measures[i]

        
        # Each entry in the mass matrix for order p0
        # Mass matrix entry is length_i for p0 elements
        mass_matrix[2*i-1 , 2*i-1] = length_i  # DOF 1
        mass_matrix[2*i , 2*i] = length_i      # DOF 2
    end
    
    return mass_matrix
end

function asm_bnd_mass_matrix_on_crack_P0(TheLine::linefrac2D)::RealMat
    # Get the number of elements in the line
    Ne = nbelements(TheLine)
    
    # Initialize the mass matrix (size 2*Ne x 2µNe)
    mass_matrix = zeros(RealNum, Ne, Ne)

    # Loop over each element to populate the mass matrix
    for i in 1:Ne
        # Length of the current segment
        length_i = TheLine.element_measures[i]
        
        # Each entry in the mass matrix for order p0
        # Mass matrix entry is length_i for p0 elements
        mass_matrix[i, i] = length_i 
    end
    
    return mass_matrix
end


#####################################################################################
# Basic stiffness matrix assembly
#
# M = [integral(Grad Phi_i(x) \dot Grad Phi_j(x)dx)]_ij
#
# where {Phi_i} are the P1 shape functions on the mesh
#
# On each triangle of the mesh, the integrals are explicitely known :
# k_ij = c_ij / (4|T|)
# where |T| is the surface of the triangle
# and C is a 3x3 matrix given by C = G*G' (' = transposed)
# where the lines of G are
#       | (y1-y3)   (x3-x2) |
# G =   | (y3-y1)   (x1-x3) |
#       | (y1-y2)   (x2-x1) |
#
function asm_stiff_2D_P1(TheMesh::volcmesh)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - K : the stiffness matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    nbtriangles = nbelements(TheMesh)
    TheTriangle = zeros(Index, (3,1))

    K = spzeros(nbdofs,nbdofs)

    x1 = zeros(RealNum, (2,1))
    x2 = zeros(RealNum, (2,1))
    x3 = zeros(RealNum, (2,1)) # Nodes of a triangle
    XX = zeros(RealNum, (3,1))
    YY = zeros(RealNum, (3,1))
    CX = zeros(RealNum, (3,3))
    CY = zeros(RealNum, (3,3))
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        x1 = TheMesh.nodes[:, TheTriangle[1] ]
        x2 = TheMesh.nodes[:, TheTriangle[2] ]
        x3 = TheMesh.nodes[:, TheTriangle[3] ]
        #
        # Compute the surface of the triangle
        #
        TheCoefficient = 1/(4*TheMesh.element_measures[index])   # This 1 / (4 x |T|)
        #
        # Build the matrices C (not calling np.array to avoid memory reallocations at each element)
        #
        YY[1] = x2[2] - x3[2]
        YY[2] = x3[2] - x1[2]
        YY[3] = x1[2] - x2[2]
        XX[1] = x3[1] - x2[1]
        XX[2] = x1[1] - x3[1]
        XX[3] = x2[1] - x1[1]
        #
        CX = XX*XX'
        CY = YY*YY'
        #
        #
        # Assemble the elementary stiffness terms
        #
        for lnode in 1:3
            for cnode in 1:3
                K[TheTriangle[lnode],TheTriangle[cnode]] += TheCoefficient * (CX[lnode,cnode]+CY[lnode,cnode])
            end # for cnode
        end # for lnode
    end # for index

    return(K)
end  # of asm_stiff_2D_P1

# Stiffness matrix assembly for problems in the form
#
#   -div(a(x) Grad u(x)) = f(x)
#
# where a(x) is a SCALAR function
#
# In this version a is given by a meshvector
#
# M = [integral(Grad Phi_i(x) \dot Grad Phi_j(x)dx)]_ij
#
# where {Phi_i} are the P1 shape functions on the mesh
#
# On each triangle of the mesh, the integrals are explicitely known :
# k_ij = c_ij / (4|T|)
# where |T| is the surface of the triangle
# and C is a 3x3 matrix given by C = G*G' (' = transposed)
# where the lines of G are
#       | (y1-y3)   (x3-x2) |
# G =   | (y3-y1)   (x1-x3) |
#       | (y1-y2)   (x2-x1) |
#
function asm_stiff_2D_P1(TheMesh::volcmesh, a::RealVec)
    # """This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # - a : a "meshvector" (vector defining the values of a on the mesh nodes)
    # Output:
    # - K : the stiffness matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    nbtriangles = nbelements(TheMesh)
    TheTriangle = zeros(Index, (3,1))

    K = spzeros(nbdofs,nbdofs)

    x1 = zeros(RealNum, (2,1))
    x2 = zeros(RealNum, (2,1))
    x3 = zeros(RealNum, (2,1)) # Nodes of a triangle
    XX = zeros(RealNum, (3,1))
    YY = zeros(RealNum, (3,1))
    CX = zeros(RealNum, (3,3))
    CY = zeros(RealNum, (3,3))
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        x1 = TheMesh.nodes[:, TheTriangle[1] ]
        x2 = TheMesh.nodes[:, TheTriangle[2] ]
        x3 = TheMesh.nodes[:, TheTriangle[3] ]
        #
        # Compute the coefficient involving the weight and the surface of the triangle
        #
        abar = (a[TheTriangle[1]] + a[TheTriangle[2]] + a[TheTriangle[3]]) / 3.0    # approx value of a at the center of element
        TheCoefficient = abar/(4*TheMesh.element_measures[index])   # This abar / (4 x |T|)
        #
        # Build the matrices C (not calling np.array to avoid memory reallocations at each element)
        #
        YY[1] = x2[2] - x3[2]
        YY[2] = x3[2] - x1[2]
        YY[3] = x1[2] - x2[2]
        XX[1] = x3[1] - x2[1]
        XX[2] = x1[1] - x3[1]
        XX[3] = x2[1] - x1[1]
        #
        CX = XX*XX'
        CY = YY*YY'
        #
        #
        # Assemble the elementary stiffness terms
        #
        for lnode in 1:3
            for cnode in 1:3
                K[TheTriangle[lnode],TheTriangle[cnode]] += TheCoefficient * (CX[lnode,cnode]+CY[lnode,cnode])
            end # for cnode
        end # for lnode
    end # for index

    return(K)
end  # of assemb_stiffness_2D_P1 (with coefficients vector)

# Stiffness matrix assembly for problems in the form
#
#   -div(a(x) Grad u(x)) = f(x)
#
# where a(x) is a SCALAR function
#
# In this version a is given by a function taking values on the mesh nodes
#
# M = [integral(Grad Phi_i(x) \dot Grad Phi_j(x)dx)]_ij
#
# where {Phi_i} are the P1 shape functions on the mesh
#
# On each triangle of the mesh, the integrals are explicitely known :
# k_ij = c_ij / (4|T|)
# where |T| is the surface of the triangle
# and C is a 3x3 matrix given by C = G*G' (' = transposed)
# where the lines of G are
#       | (y1-y3)   (x3-x2) |
# G =   | (y3-y1)   (x1-x3) |
#       | (y1-y2)   (x2-x1) |
#
function asm_stiff_2D_P1(TheMesh::volcmesh, afonk::Function)
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # - a : a function defined on the mesh nodes and returning a scalar real value
    # Output:
    # - K : the stiffness matrix on the mesh for P1 elements
    #
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs = nbnodes(TheMesh)
    nbtriangles = nbelements(TheMesh)
    TheTriangle = zeros(Index, (3,1))

    K = spzeros(nbdofs,nbdofs)

    x1 = zeros(RealNum, (2,1))
    x2 = zeros(RealNum, (2,1))
    x3 = zeros(RealNum, (2,1)) # Nodes of a triangle
    XX = zeros(RealNum, (3,1))
    YY = zeros(RealNum, (3,1))
    CX = zeros(RealNum, (3,3))
    CY = zeros(RealNum, (3,3))
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        x1 = TheMesh.nodes[:, TheTriangle[1] ]
        x2 = TheMesh.nodes[:, TheTriangle[2] ]
        x3 = TheMesh.nodes[:, TheTriangle[3] ]
        #
        # Compute the coefficient involving the weight and the surface of the triangle
        #
        a1 = afonk(TheMesh.nodes[:,TheTriangle[1]])
        a2 = afonk(TheMesh.nodes[:,TheTriangle[2]])
        a3 = afonk(TheMesh.nodes[:,TheTriangle[3]])
        abar = (a1 + a2 + a3) / 3.0    # approx value of a at the center of element
        TheCoefficient = abar/(4*TheMesh.element_measures[index])   # This abar / (4 x |T|)
        #
        # Build the matrices C (not calling np.array to avoid memory reallocations at each element)
        #
        YY[1] = x2[2] - x3[2]
        YY[2] = x3[2] - x1[2]
        YY[3] = x1[2] - x2[2]
        XX[1] = x3[1] - x2[1]
        XX[2] = x1[1] - x3[1]
        XX[3] = x2[1] - x1[1]
        #
        CX = XX*XX'
        CY = YY*YY'
        #
        #
        # Assemble the elementary stiffness terms
        #
        for lnode in 1:3
            for cnode in 1:3
                K[TheTriangle[lnode],TheTriangle[cnode]] += TheCoefficient * (CX[lnode,cnode]+CY[lnode,cnode])
            end # for cnode
        end # for lnode
    end # for index

    return(K)
end  # of assemb_stiffness_2D_P1 (with coefficients vector)


###########################################################################
#################                                       ###################
#################          LINEAR ELASTICITY            ###################
#################                                       ###################
###########################################################################
function asm_elast_stiff_2D_P1(TheMesh::volcmesh, Young::RealNum , nu::RealNum)   
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - K : the elasticity matrix on the mesh for P1 elements
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs::Index = 2*nbnodes(TheMesh)    # In 2D elasticity problems there are 2 dofs at each node.
    nbtriangles::Index = nbelements(TheMesh)
   
    K = spzeros(nbdofs,nbdofs)
    
    
    TheTriangle = zeros(Index, (3,1))
    x1 = zeros(RealNum, (2,1))
    x2 = zeros(RealNum, (2,1))
    x3 = zeros(RealNum, (2,1)) # Nodes of a triangle
    XX = zeros(RealNum, (3,1))
    YY = zeros(RealNum, (3,1))
    CX = zeros(RealNum, (3,3))
    CY = zeros(RealNum, (3,3))
    CXY = zeros(RealNum, (3,3))
    Rvoigt = zeros(RealNum, (6,6))
    NumUnknowns = zeros(Index, 6)
    NumDof::Index = 1
    #
    # Lamé coefficients for the operator (lambda + 2mu)
    # lambda being reserved in some languages (e.g. python) , we use lamda (c) V. Cayol :-)
    #
    lamda = (nu*Young) / ((1-2*nu)*(1+nu))
    mu = Young / (2*(1+nu))
    LameCoef = lamda + 2*mu
    
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        x1 = TheMesh.nodes[:, TheTriangle[1]]
        x2 = TheMesh.nodes[:, TheTriangle[2]]
        x3 = TheMesh.nodes[:, TheTriangle[3]]
        #
        # Compute the surface of the triangle
        #
        TheCoefficient = 1 / (4*TheMesh.element_measures[index])   # This 1 / (4 x |T|)
        lamh = lamda*TheCoefficient
        muh = mu*TheCoefficient
        Lameh = LameCoef*TheCoefficient
        #
        # Build the matrices C (not calling np.array to avoid memory reallocations at each element)
        #
        YY[1] = x2[2] - x3[2]
        YY[2] = x3[2] - x1[2]
        YY[3] = x1[2] - x2[2]
        XX[1] = x3[1] - x2[1]
        XX[2] = x1[1] - x3[1]
        XX[3] = x2[1] - x1[1]
        #
        CX = XX*XX'
        CY = YY*YY'
        CXY = XX*YY'
        #
        #
        # Compute the elementary elasticity terms
        # Diagonal blocks (trying to save cputime)
        Rvoigt[1:3,1:3] = Lameh*CY + muh*CX
        Rvoigt[4:6,4:6] = Lameh*CX + muh*CY
        Rvoigt[1:3,4:6] = lamh*CXY' + muh*CXY
        Rvoigt[4:6,1:3] = lamh*CXY + muh*CXY'
        #
        # Assemble the term, reordering the unknowns in "standard order"
        #
        # First store the dofs number in "Voigt" order
        # Remind that the number of the triangles start from 1
        #
        # So at node T[nk] correspond 2 consecutive dofs in the system : 2*T[nk] - 2 and 2*T[Nk] - 1
        # The order of the unknowns in the Voigt formulation is 1 3 5 2 4 6
        for nk in 1:3
            NumDof = 2*TheTriangle[nk] - 1  
            NumUnknowns[nk] = NumDof
            NumUnknowns[nk+3] = NumDof + 1
        end
        # Now assemble the matrix terms
        # Scanning the matrix Rvoigt and using the NumUnknowns permutation table
        #
        for ldof in 1:6
            for cdof in 1:6
                K[NumUnknowns[ldof],NumUnknowns[cdof]] += Rvoigt[ldof,cdof]
            end # cnode
        end # for lnode
    end # for index

    return(K)
end # of elast_stiffness

function asm_elast_stiff_2D_P1(TheMesh::volcmesh, Young::RealNum , nu::RealNum, afonk::Function)   
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - K : the elasticity matrix on the mesh for P1 elements
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs::Index = 2*nbnodes(TheMesh)    # In 2D elasticity problems there are 2 dofs at each node.
    nbtriangles::Index = nbelements(TheMesh)
   
    K = spzeros(nbdofs,nbdofs)
    
    
    TheTriangle = zeros(Index, (3,1))
    x1 = zeros(RealNum, (2,1))
    x2 = zeros(RealNum, (2,1))
    x3 = zeros(RealNum, (2,1)) # Nodes of a triangle
    XX = zeros(RealNum, (3,1))
    YY = zeros(RealNum, (3,1))
    CX = zeros(RealNum, (3,3))
    CY = zeros(RealNum, (3,3))
    CXY = zeros(RealNum, (3,3))
    Rvoigt = zeros(RealNum, (6,6))
    NumUnknowns = zeros(Index, 6)
    NumDof::Index = 1
    #
    # Lamé coefficients for the operator (lambda + 2mu)
    # lambda being reserved in some languages (e.g. python) , we use lamda (c) V. Cayol :-)
    #
    lamda = (nu*Young) / ((1-2*nu)*(1+nu))
    mu = Young / (2*(1+nu))
    LameCoef = lamda + 2*mu
    
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        x1 = TheMesh.nodes[:, TheTriangle[1]]
        x2 = TheMesh.nodes[:, TheTriangle[2]]
        x3 = TheMesh.nodes[:, TheTriangle[3]]
        #
        # Compute the surface of the triangle
        #
        a1 = afonk(TheMesh.nodes[:,TheTriangle[1]])
        a2 = afonk(TheMesh.nodes[:,TheTriangle[2]])
        a3 = afonk(TheMesh.nodes[:,TheTriangle[3]])
        abar = (a1 + a2 + a3) / 3.0

        TheCoefficient = abar / (4*TheMesh.element_measures[index])   # This 1 / (4 x |T|)
        lamh = lamda*TheCoefficient
        muh = mu*TheCoefficient
        Lameh = LameCoef*TheCoefficient
        #
        # Build the matrices C (not calling np.array to avoid memory reallocations at each element)
        #
        YY[1] = x2[2] - x3[2]
        YY[2] = x3[2] - x1[2]
        YY[3] = x1[2] - x2[2]
        XX[1] = x3[1] - x2[1]
        XX[2] = x1[1] - x3[1]
        XX[3] = x2[1] - x1[1]
        #
        CX = XX*XX'
        CY = YY*YY'
        CXY = XX*YY'
        #
        #
        # Compute the elementary elasticity terms
        # Diagonal blocks (trying to save cputime)
        Rvoigt[1:3,1:3] = Lameh*CY + muh*CX
        Rvoigt[4:6,4:6] = Lameh*CX + muh*CY
        Rvoigt[1:3,4:6] = lamh*CXY' + muh*CXY
        Rvoigt[4:6,1:3] = lamh*CXY + muh*CXY'
        #
        # Assemble the term, reordering the unknowns in "standard order"
        #
        # First store the dofs number in "Voigt" order
        # Remind that the number of the triangles start from 1
        #
        # So at node T[nk] correspond 2 consecutive dofs in the system : 2*T[nk] - 2 and 2*T[Nk] - 1
        # The order of the unknowns in the Voigt formulation is 1 3 5 2 4 6
        for nk in 1:3
            NumDof = 2*TheTriangle[nk] - 1  
            NumUnknowns[nk] = NumDof
            NumUnknowns[nk+3] = NumDof + 1
        end
        # Now assemble the matrix terms
        # Scanning the matrix Rvoigt and using the NumUnknowns permutation table
        #
        for ldof in 1:6
            for cdof in 1:6
                K[NumUnknowns[ldof],NumUnknowns[cdof]] += Rvoigt[ldof,cdof]
            end # cnode
        end # for lnode
    end # for index

    return(K)
end # of elast_stiffness

function asm_elast_stiff_2D_P1_lambda_mu(TheMesh::volcmesh, lamda::RealNum , mu::RealNum, afonk::Function)   
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - K : the elasticity matrix on the mesh for P1 elements
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs::Index = 2*nbnodes(TheMesh)    # In 2D elasticity problems there are 2 dofs at each node.
    nbtriangles::Index = nbelements(TheMesh)
   
    K = spzeros(nbdofs,nbdofs)
    
    
    TheTriangle = zeros(Index, (3,1))
    x1 = zeros(RealNum, (2,1))
    x2 = zeros(RealNum, (2,1))
    x3 = zeros(RealNum, (2,1)) # Nodes of a triangle
    XX = zeros(RealNum, (3,1))
    YY = zeros(RealNum, (3,1))
    CX = zeros(RealNum, (3,3))
    CY = zeros(RealNum, (3,3))
    CXY = zeros(RealNum, (3,3))
    Rvoigt = zeros(RealNum, (6,6))
    NumUnknowns = zeros(Index, 6)
    NumDof::Index = 1
    #
    # Lamé coefficients for the operator (lambda + 2mu)
    # lambda being reserved in some languages (e.g. python) , we use lamda (c) V. Cayol :-)
    #
    LameCoef = lamda + 2*mu
    
    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        x1 = TheMesh.nodes[:, TheTriangle[1]]
        x2 = TheMesh.nodes[:, TheTriangle[2]]
        x3 = TheMesh.nodes[:, TheTriangle[3]]
        #
        # Compute the surface of the triangle
        #
        a1 = afonk(TheMesh.nodes[:,TheTriangle[1]])
        a2 = afonk(TheMesh.nodes[:,TheTriangle[2]])
        a3 = afonk(TheMesh.nodes[:,TheTriangle[3]])
        abar = (a1 + a2 + a3) / 3.0

        TheCoefficient = abar / (4*TheMesh.element_measures[index])   # This 1 / (4 x |T|)
        lamh = lamda*TheCoefficient
        muh = mu*TheCoefficient
        Lameh = 0.0
        #
        # Build the matrices C (not calling np.array to avoid memory reallocations at each element)
        #
        YY[1] = x2[2] - x3[2]
        YY[2] = x3[2] - x1[2]
        YY[3] = x1[2] - x2[2]
        XX[1] = x3[1] - x2[1]
        XX[2] = x1[1] - x3[1]
        XX[3] = x2[1] - x1[1]
        #
        CX = XX*XX'
        CY = YY*YY'
        CXY = XX*YY'
        #
        #
        # Compute the elementary elasticity terms
        # Diagonal blocks (trying to save cputime)
        Rvoigt[1:3,1:3] = Lameh*CY + muh*CX
        Rvoigt[4:6,4:6] = Lameh*CX + muh*CY
        Rvoigt[1:3,4:6] = lamh*CXY' + muh*CXY
        Rvoigt[4:6,1:3] = lamh*CXY + muh*CXY'
        #
        # Assemble the term, reordering the unknowns in "standard order"
        #
        # First store the dofs number in "Voigt" order
        # Remind that the number of the triangles start from 1
        #
        # So at node T[nk] correspond 2 consecutive dofs in the system : 2*T[nk] - 2 and 2*T[Nk] - 1
        # The order of the unknowns in the Voigt formulation is 1 3 5 2 4 6
        for nk in 1:3
            NumDof = 2*TheTriangle[nk] - 1  
            NumUnknowns[nk] = NumDof
            NumUnknowns[nk+3] = NumDof + 1
        end
        # Now assemble the matrix terms
        # Scanning the matrix Rvoigt and using the NumUnknowns permutation table
        #
        for ldof in 1:6
            for cdof in 1:6
                K[NumUnknowns[ldof],NumUnknowns[cdof]] += Rvoigt[ldof,cdof]
            end # cnode
        end # for lnode
    end # for index

    return(K)
end # of elast_stiffness


function asm_elast_stiff_2D_P1(TheMesh::volcmesh, Young::RealVec , nu::RealVec)   
    # This routine builds the mass matrix for P1 elements on a mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # - Young, nu : vectors giving the pointwise values of the Young modulus and Poisson coefficient on the nodes in the mesh
    # Output:
    # - K : the elasticity matrix on the mesh for P1 elements
    #
    # Check if the element measures are present
    #
    if !TheMesh.elm_meas_flag
        mk_element_measures(TheMesh)
    end
    # First read the number of dofs and declare the matrix
    # Init. various variables
    #
    nbdofs::Index = 2*nbnodes(TheMesh)    # In 2D elasticity problems there are 2 dofs at each node.
    nbtriangles::Index = nbelements(TheMesh)
   
    K = spzeros(nbdofs,nbdofs)
    
    
    TheTriangle = zeros(Index, (3,1))
    x1 = zeros(RealNum, (2,1))
    x2 = zeros(RealNum, (2,1))
    x3 = zeros(RealNum, (2,1)) # Nodes of a triangle
    XX = zeros(RealNum, (3,1))
    YY = zeros(RealNum, (3,1))
    CX = zeros(RealNum, (3,3))
    CY = zeros(RealNum, (3,3))
    CXY = zeros(RealNum, (3,3))
    Rvoigt = zeros(RealNum, (6,6))
    NumUnknowns = zeros(Index, 6)
    NumDof::Index = 1
    #
    # Build meshvectors for the Lamé coefficients
    #
    # Lamé coefficients for the operator (lambda + 2mu)
    # lambda being reserved in some languages (e.g. python) , we use lamda (c) V. Cayol :-)
    #
    # Using vector brodcasting of usual operations

    # Assuming nu and Young are vectors of the same length

    lamda = (nu .* Young) ./ ((1 .- 2 .* nu) .* (1 .+ nu))  # Use .- and .+ for element-wise operations
    mu = Young ./ (2 .* (1.0 .+ nu))  # Use .+ for element-wise addition
    LameCoef = lamda .+ 2 .* mu  # Use .+ for element-wise addition

    #
    # Loop over triangles and build the elementary matrix
    #
    for index in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,index]
        x1 = TheMesh.nodes[:, TheTriangle[1]]
        x2 = TheMesh.nodes[:, TheTriangle[2]]
        x3 = TheMesh.nodes[:, TheTriangle[3]]
        #
        # Compute the mean Lamé coefficient on the triangle (value at the center of the triangle)
        #
        lamdabar = (lamda[TheTriangle[1]] + lamda[TheTriangle[2]] + lamda[TheTriangle[3]]) / 3.0
        mubar = (mu[TheTriangle[1]] + mu[TheTriangle[2]] + mu[TheTriangle[3]]) / 3.0
        LameCoefbar = (LameCoef[TheTriangle[1]] + LameCoef[TheTriangle[2]] + LameCoef[TheTriangle[3]]) / 3.0
        #
        # Compute the coefficient involving surface of the triangle
        #
        TheCoefficient = 1 / (4*TheMesh.element_measures[index])   # This 1 / (4 x |T|)
        lamh = lamdabar * TheCoefficient
        muh = mubar * TheCoefficient
        Lameh = LameCoefbar * TheCoefficient
        #
        # Build the matrices C (not calling np.array to avoid memory reallocations at each element)
        #
        YY[1] = x2[2] - x3[2]
        YY[2] = x3[2] - x1[2]
        YY[3] = x1[2] - x2[2]
        XX[1] = x3[1] - x2[1]
        XX[2] = x1[1] - x3[1]
        XX[3] = x2[1] - x1[1]
        #
        CX = XX*XX'
        CY = YY*YY'
        CXY = XX*YY'
        #
        #
        # Compute the elementary elasticity terms
        # Diagonal blocks (trying to save cputime)
        Rvoigt[1:3,1:3] = Lameh*CY + muh*CX
        Rvoigt[4:6,4:6] = Lameh*CX + muh*CY
        Rvoigt[1:3,4:6] = lamh*CXY' + muh*CXY
        Rvoigt[4:6,1:3] = lamh*CXY + muh*CXY'
        #
        # Assemble the term, reordering the unknowns in "standard order"
        #
        # First store the dofs number in "Voigt" order
        # Remind that the number of the triangles start from 1
        #
        # So at node T[nk] correspond 2 consecutive dofs in the system : 2*T[nk] - 2 and 2*T[Nk] - 1
        # The order of the unknowns in the Voigt formulation is 1 3 5 2 4 6, or, starting from 0 :
        # 0 2 4 1 3 5
        for nk in 1:3
            NumDof = 2*TheTriangle[nk] - 1  #
            NumUnknowns[nk] = NumDof
            NumUnknowns[nk+3] = NumDof + 1
        end
        # Now assemble the matrix terms
        # Scanning the matrix Rvoigt and using the NumUnknowns permutation table
        #
        for ldof in 1:6
            for cdof in 1:6
                K[NumUnknowns[ldof],NumUnknowns[cdof]] += Rvoigt[ldof,cdof]
            end # cnode
        end # for lnode
    end # for index

    return(K)
end # of elast_stiffness


#################################################################################
##################                                                  #############
##################  Specific routines for cross matrices between    #############
##################  main meshes and fractures                       #############
##################                                                  #############
#################################################################################


function asm_mass_mesh_frac_2D_P1P0(TheMesh::volcmesh, TheLine::linefrac2D)
    #
    #
    # le P1 pour les mult et p0 pour les u ?

    #### Mettre tests conformité divers et exception throws au kazou
    ### faire une version avec les params d'intersection

    # Init. 

    Ne = nbnodes(TheMesh)   # Number of lines of the resulting matrix
    Nf = nbelements(TheLine) # Number of columns

    Bmat = spzeros(Ne,Nf)

    # Compute the intersection points of the mesh and the segments

    (HasIntersection,InterSectPoints,Seglengths) = frac_mesh_intersect(TheLine)

    # Small arrays used throughout the routine

    p1 = zeros(RealNum,(2,1)) 
    p2 = zeros(RealNum,(2,1))
    Y = zeros(RealNum,(2,1))
    lg = zeros(RealNum,(2,1))
    PhiX = zeros(RealNum,3)
    PhiY = zeros(RealNum,3)
    IntegralZzz = zeros(RealNum,3)

    #
    # Loop over the segments of TheLine
    #

    for jndex = 1:Nf

        CurrentElem = TheLine.elements[:,jndex]

        # Get the start and end points of the current segment
        p1 = TheLine.nodes[:,CurrentElem[1]] # First point
        p2 = TheLine.nodes[:,CurrentElem[2]] # Last point

        # Get the triangles of the mesh containing each point

        Tr1 = TheLine.embed_triangles[CurrentElem[1]]
        Tr2 = TheLine.embed_triangles[CurrentElem[2]]

        if (HasIntersection[jndex]) # If the segment is split between two triangles, 2 terms

            Y = InterSectPoints[:,jndex]
            lg = Seglengths[:,jndex]

            ############
            # Triangle 1
            ############
            
            TheTriangle = TheMesh.elements[:,Tr1]
            # x1 = TheMesh.nodes[:, TheTriangle[1]]
            # x2 = TheMesh.nodes[:, TheTriangle[2]]
            # x3 = TheMesh.nodes[:, TheTriangle[3]]
            
            # Evaluate basis function values at P1 and Y
            # Remember that Y lies in both triangles

            basis_func_pt_val_2D_P1_DEV!(p1,TheMesh,Tr1,PhiX)
            basis_func_pt_val_2D_P1_DEV!(Y,TheMesh,Tr1,PhiY)

            # Compute the triangle formula on the segment

            IntegralZzz = 0.5*lg[1]*(PhiX+PhiY)  # Warning P0 on the fracture so that the basis function on the seg is 1

            # Store the results in Bmat

            Bmat[TheTriangle[1],jndex] = Bmat[TheTriangle[1],jndex] + IntegralZzz[1]
            Bmat[TheTriangle[2],jndex] = Bmat[TheTriangle[2],jndex] + IntegralZzz[2]
            Bmat[TheTriangle[3],jndex] = Bmat[TheTriangle[3],jndex] + IntegralZzz[3]

            ############
            # Triangle 2
            ############
            
            TheTriangle = TheMesh.elements[:,Tr2]
            
            # Evaluate basis function values at P1 and Y
            # Remember that Y lies in both triangles

            basis_func_pt_val_2D_P1_DEV!(p2,TheMesh,Tr2,PhiX)
            basis_func_pt_val_2D_P1_DEV!(Y,TheMesh,Tr2,PhiY)

            # Compute the triangle formula on the segment

            IntegralZzz = 0.5*lg[2]*(PhiX+PhiY)

            # Store the results in Bmat

            Bmat[TheTriangle[1],jndex] = Bmat[TheTriangle[1],jndex] + IntegralZzz[1]
            Bmat[TheTriangle[2],jndex] = Bmat[TheTriangle[2],jndex] + IntegralZzz[2]
            Bmat[TheTriangle[3],jndex] = Bmat[TheTriangle[3],jndex] + IntegralZzz[3]

        else # Here the segment lies entirely in one triangle

            TheTriangle = TheMesh.elements[:,Tr1]
            
            # Evaluate basis function values at P1 and Y
            # Remember that Y lies in both triangles

            basis_func_pt_val_2D_P1_DEV!(p1,TheMesh,Tr1,PhiX)
            basis_func_pt_val_2D_P1_DEV!(p2,TheMesh,Tr1,PhiY)

            # Compute the triangle formula on the segment

            IntegralZzz = 0.5*TheLine.element_measures[jndex]*(PhiX+PhiY)

            # Store the results in Bmat

            Bmat[TheTriangle[1],jndex] = Bmat[TheTriangle[1],jndex] + IntegralZzz[1]
            Bmat[TheTriangle[2],jndex] = Bmat[TheTriangle[2],jndex] + IntegralZzz[2]
            Bmat[TheTriangle[3],jndex] = Bmat[TheTriangle[3],jndex] + IntegralZzz[3]
            # println("Triangle nodes : ", TheTriangle)
            # println("       p1 = ", p1, " --- Y = ", Y)
            # println("       PhiX = ", PhiX , " --- PhiY = ", PhiY," --- IntegralZzz = ", IntegralZzz)

        end
    end

    return(Bmat)

end

function asm_cross_mass_tensor_2D_P1P0(TheMesh::volcmesh, TheLine::linefrac2D, Young::Real, nu::Real)
    # Initialization
    Ne = nbnodes(TheMesh)   # Number of nodes
    Nf = nbelements(TheLine) # Number of segments

    Bmat = spzeros(2*Ne, 2*Nf)   # Sparse matrix for results (2D P1/P0)

    # Compute the intersection points of the mesh and the segments
    (HasIntersection, InterSectPoints, Seglengths) = frac_mesh_intersect(TheLine)

    lamda = (nu*Young) / ((1-2*nu)*(1+nu))
    mu = Young / (2*(1+nu))
    LameCoef = lamda + 2*mu

    # Compute the segment normals if necessary

    if !TheLine.elm_norm_flag
        normals = mk_segment_normals(TheLine) # normals is a boolean
    end
       
    IntegralZzz = zeros(RealNum, (6,2))
    NormStresses = zeros(2,6) # Normal stresses in ONE element (6 of them since P1 elements)
    NormalVec = zeros(2,1) # Normal vector on ONE segment of the linfe/rerservoir

    # Loop over the segments of TheLine
    for jndex = 1:Nf
        CurrentElem = TheLine.elements[:, jndex]

        # Get the start and end points of the current segment
        p1 = TheLine.nodes[:, CurrentElem[1]] # First point
        p2 = TheLine.nodes[:, CurrentElem[2]] # Last point

        # Get the triangles of the mesh containing each point
        Tr1 = TheLine.embed_triangles[CurrentElem[1]]
        Tr2 = TheLine.embed_triangles[CurrentElem[2]]

        if (HasIntersection[jndex]) # If the segment intersects two triangles
            Y = InterSectPoints[:, jndex]
            lg = Seglengths[:, jndex]

            ############
            # Triangle 1
            ############

            TheTriangle = TheMesh.elements[:,Tr1]

            _, stress_tensorP  = basis_func_stress_tensor_pt_val_2D_P1(p1, TheMesh, mu, lamda)
            # _, stress_tensorY  = basis_func_stress_tensor_pt_val_2D_P1(Y, TheMesh, nu, lambda)

            # get the normal unit vectors of the segment
            NormalVec = TheLine.element_normals[:,jndex]

            # Compute the stress tensors applied to the normal vector
            for iloop in 1:6
                NormStresses[:,iloop] = stress_tensorP[:,:,iloop]*NormalVec 
            end
           
            # Compute integral contributions (basis functions on the segment are
            IntegralZzz = lg[1]*(NormStresses'*Matrix(I,2,2)) 

            # Store the results in Bmat
            for i in 1:3
                Bmat[2 * TheTriangle[i] - 1, 2 * jndex - 1] += IntegralZzz[2i-1,1]
                Bmat[2 * TheTriangle[i] - 1, 2 * jndex] += IntegralZzz[2i-1,2]
                Bmat[2 * TheTriangle[i], 2 * jndex - 1] += IntegralZzz[2i,1] 
                Bmat[2 * TheTriangle[i], 2 * jndex] += IntegralZzz[2i,2]
            end

            ############
            # Triangle 2
            ############

            TheTriangle = TheMesh.elements[:,Tr2]

            _, stress_tensorP = basis_func_stress_tensor_pt_val_2D_P1(p2, TheMesh, mu, lamda)
            # _, stress_tensorY = basis_func_stress_tensor_pt_val_2D_P1(Y, TheMesh, nu, lambda)

           # get the normal unit vectors of the segment
           NormalVec = TheLine.element_normals[:,jndex]

           # Compute the stress tensors applied to the normal vector
           for iloop in 1:6
               NormStresses[:,iloop] = stress_tensorP[:,:,iloop]*NormalVec 
           end
          
           # Compute integral contributions (basis functions on the segment are
           IntegralZzz = lg[2]*(NormStresses'*Matrix(I,2,2)) 

           # Store the results in Bmat
           for i in 1:3
               Bmat[2 * TheTriangle[i] - 1, 2 * jndex - 1] += IntegralZzz[2i-1,1]
               Bmat[2 * TheTriangle[i] - 1, 2 * jndex] += IntegralZzz[2i-1,2]
               Bmat[2 * TheTriangle[i], 2 * jndex - 1] += IntegralZzz[2i,1] 
               Bmat[2 * TheTriangle[i], 2 * jndex] += IntegralZzz[2i,2]
           end

        else # Segment lies entirely in one triangle

            TheTriangle = TheMesh.elements[:, Tr1]


            # Since the degree of the elements is P1 we don't need to compute
            # the stress tensors at both end of the segment
            _, stress_tensorP1 = basis_func_stress_tensor_pt_val_2D_P1(p1, TheMesh, mu, lamda)

            # get the normal unit vectors of the segment
            NormalVec = TheLine.element_normals[:,jndex]

            # Compute the stress tensors applied to the normal vector
            for iloop in 1:6
                NormStresses[:,iloop] = stress_tensorP1[:,:,iloop]*NormalVec
            end

            # Compute integral contributions (basis functions on the segment are
            # (1,0) and (0,1) since we use P0 approx) via rectangle formula
            IntegralZzz = TheLine.element_measures[jndex]*(NormStresses'*Matrix(I,2,2)) # taille (6,2)

            # print(IntegralZzz)
            # println(TheTriangle)

           # Store the results in Bmat
           for i in 1:3
                Bmat[2 * TheTriangle[i] - 1, 2 * jndex - 1] += IntegralZzz[2i-1,1]
                Bmat[2 * TheTriangle[i] - 1, 2 * jndex] += IntegralZzz[2i-1,2]
                Bmat[2 * TheTriangle[i], 2 * jndex - 1] += IntegralZzz[2i,1] 
                Bmat[2 * TheTriangle[i], 2 * jndex] += IntegralZzz[2i,2]
           end

         end
    end

    return Bmat
end

# Function to interpolate the solution at point x inside a triangle
function interpolate(U, TheMesh::volcmesh, x::RealVec)
    
    TheTri = intriangle(x, TheMesh)

    # Get the coordinates of the nodes of the current triangle
    x1 =  TheMesh.elements[:,TheTri][1]
    x2 =  TheMesh.elements[:,TheTri][2]
    x3 =  TheMesh.elements[:,TheTri][3]

    _,phi,_ = basis_func_pt_val_2D_P1(x, TheMesh)
 
    # Linear interpolation using the barycentric coordinates
    u_x = phi[1] * U[x1] + phi[2]* U[x2] + phi[3] * U[x3]

    return u_x
end


function Interpolate_in_another_mesh(TheMesh::volcmesh, TheImmersedMesh::volcmesh, U_immersed)

    N = nbnodes(TheMesh)
    U_interpolate = zeros(N)

    for index = 1:N
        x = TheMesh.nodes[:, index]
        U_interpolate[index] = interpolate(U_immersed, TheImmersedMesh,  x)
    end

    return U_interpolate

end

function interpolate_2DOF(U, TheMesh::volcmesh, x::RealVec, dof::Index)
    
    TheTri = intriangle(x, TheMesh)

    # Get the coordinates of the nodes of the current triangle
    x1 =  TheMesh.elements[:,TheTri][1]
    x2 =  TheMesh.elements[:,TheTri][2]
    x3 =  TheMesh.elements[:,TheTri][3]

    _,phi,_ = basis_func_pt_val_2D_P1(x, TheMesh)
 
    # Linear interpolation using the barycentric coordinates
    u_x = phi[1] * U[x1 + dof - 1] + phi[2]* U[x2 + dof - 1] + phi[3] * U[x3 + dof - 1]

    return u_x
end

function Interpolate_in_another_mesh_2DOF(TheMesh::volcmesh, TheImmersedMesh::volcmesh, U_immersed)

    N = nbnodes(TheMesh)
    U_interpolate = zeros(2*N)

    for index = 1:N
        x = TheMesh.nodes[:, index]
        U_interpolate[2*index-1] = interpolate_2DOF(U_immersed, TheImmersedMesh,  x, 1)
        U_interpolate[2*index] = interpolate_2DOF(U_immersed, TheImmersedMesh,  x, 2)
    end

    return U_interpolate

end

