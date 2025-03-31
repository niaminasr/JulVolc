# This is the file meshtools.jl
# In this file the classes and types for mesh handling are defined
#
# In julvolc meshes are handled as a structure containing the number of nodes and 3 ls :matrixts, edges and elements
# (triangles in 2D and tetrahedra in 3D)
# Fractures are also defined in this file
#
# a fracture is basicly a mesh embedded in another with higher dimension
#
# (c) Olivier Bodart - 2021/09/22 --> 
#
# New functions by Niami Nasr Nov. 2024
# meshdivergencevector
# Crée un tableau contenant la divergence d'une fonction tensorielle évaluée aux nœuds d'un maillage.
#
# linevector
# Crée un tableau contenant la valeur d'une fonction évaluée aux nœuds des segments de ligne.
#
# midlinevector
# Crée un tableau contenant la valeur d'une fonction évaluée aux points milieux des segments de ligne.
#
# plot2Dlinewithnormals
# Trace tous les segments d'un objet linfrac2D sur le maillage d'enrobage.
#
# New function by O. Bodart Dec. 2024
#
# move2Dmesh : moves the points of a mesh according to a displacement field and stores the result
# into a new volcmesh object
#


include("commons.jl")

###########################################################################
###########################################################################
###########################################################################
# The mesh class in 2D or 3D
###########################################################################
###########################################################################
###########################################################################


mutable struct volcmesh
    #
    # Main data (the mesh is stored in (points, edges,elements) format)
    #
    nodes::RealMat
    edges::IndexMat
    elements::IndexMat
    #
    # Secondary Data (computed by called methods)
    #
    element_measures::RealVec
    edge_measures::RealVec
    # Flags
    #
    isempty::Bool
    elm_meas_flag::Bool
    edg_meas_flag::Bool
    #
    ## Default constructor(creates an empty mesh to be filled later)
    volcmesh() = new(@EmptyRealMat,@EmptyIntMat,@EmptyIntMat,@EmptyRealVec,@EmptyRealVec,true,false,false)
    #
    ## Constructor using p,e,t arrays as an input
    function volcmesh(points::RealMat,edges::IndexMat,triangles::IndexMat)
        # Constructor for the volcmesh class
        #     Optional input arguments are
        #     - nodes : points array i.e. the vertices of the tesselation (each column is a point)
        #     - edges : edges array (each column is a boundary edge in 2D, a triangle in 3D)
        #     - t : segments triangles or tetrahaedra array (nodes in trigonometric order in 2D
        #     each column has Dimension+1 node numbers and one region number)
        #     If the dimension of the arrays do not match an empty object is defined
        #     An edge consists in 2 or 3 node numbers plus one boundary number. Therefore
        #     each column in the edges array has Dimension+1 elements, where Dimension is
        #     the space dimension.
        #
        # WARNING : the data in edges and Triangles are automatically converted into integer values
        # (Index aka Int64)
        # Dimensions are checked and the points array is verified to contain Real values

        Dimension = size(points,1) # each column in points is a set of coordinates
        DimEdges = size(edges,1) - 1 # Edges has one more row than the space dimension (boundary number)
        print(Dimension,"  ",DimEdges)
        
        if (Dimension == DimEdges) && (eltype(points) == RealNum)
            cedges = convert(IndexMat,edges) 
            ctriangles = convert(IndexMat,triangles)
            return new(points,cedges,ctriangles,@EmptyRealVec,@EmptyRealVec,false,false,false)
        else
            println("ERROR -- volcmesh constructor\\ Incompatible dimensions or bad data types\\ Creating an empty mesh")
            return new(@EmptyRealMat,@EmptyIntMat,@EmptyIntMat,@EmptyRealVec,@EmptyRealVec,true,false,false)
        end
    end

    #
    # Constructor importing the mesh from a matlab or gmsh file
    #
    function volcmesh(FileName::String, FileType::Index)
        # Constructor for the volcmesh class from an external file
        #
        # Currently supported : Matlab
        # Soon supported : GMSH
        #
        # Input parameters :
        #
        # FileName : a string containing a valid (existing) file name
        # FileType : 1 for Matlab, 2 for gmsh,...
        #
        # No case or switch instruction existing in Julia at the moment (2021-09)
        # we proceed with brute force

        if (FileType == 1)
            (points,edges,triangles) = MatlabMeshImport(FileName)
        elseif (FileType == 2)
            (points,edges,triangles) = GmshMeshImport(FileName)
        else
            points = @EmptyRealMat
            edges = @EmptyIntMat
            triangles = @EmptyIntMat
        end

        if (length(points) == 0) || (length(edges) == 0) || (length(triangles) == 0) # if one array is empty, create an empty mesh
            return new(@EmptyRealMat,@EmptyIntMat,@EmptyIntMat,@EmptyRealVec,@EmptyRealVec,true,false,false) 
        else
            return new(points,edges,triangles,@EmptyRealVec,@EmptyRealVec,false,false,false)
        end
    end

end # definition of the volcmesh struct

###########################################################################
###########################################################################
#########                   Inner Data Methods              ###############
###########################################################################
###########################################################################


function nbnodes(TheMesh::volcmesh)::Index
    # Return the number of nodes in the provided mesh
    return(size(TheMesh.nodes,2))
end

function nbedges(TheMesh::volcmesh)::Index
    # Return the number of boundary segments in the provided mesh
    return(size(TheMesh.edges,2))
end

function nbelements(TheMesh::volcmesh)::Index
    # Return the number of elements in the provided mesh
    return(size(TheMesh.elements,2))
end

function dimension(TheMesh::volcmesh)::Index
    # return the space dimension of the mesh
    return(size(TheMesh.nodes,1))
end

###########################################################################
###########################################################################
######### Surfaces, volumes and lengths computations        ###############
###########################################################################
###########################################################################



function mk_element_measures(TheMesh::volcmesh)::Bool
    # Computation of the element measures in a mesh
    # for 2-D or 3-D meshes
    #
    # Input : TheMesh : a non empty mesh of type volcmesh

    if (dimension(TheMesh) == 2)
        return(element_surfaces(TheMesh))
    else
        return(element_volumes(TheMesh))
    end
end # of mk_element_measures


function mk_edge_measures(TheMesh::volcmesh)::Bool
    # Computation of the boundary esges measures in a mesh
    # for 2-D or 3-D meshes
    #
    # Input : TheMesh : a non empty mesh of type volcmesh

    if (dimension(TheMesh) == 2)
        return(edge_lengths(TheMesh))
    else
        return(edge_surfaces(TheMesh))
    end
end # of mk_edge_measures

###############################
# Elements volumes and surfaces
###############################

#
# Compute the areas of the triangles of a 2-D mesh
#
function element_surfaces(TheMesh::volcmesh)::Bool
    # This routine computes the surfaces of the triangles of a 2-D mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - Surfaces : array containing the surfaces of the triangles in the mesh
    #
    # First read the number of triangles and declare the output array
    #
    nbtriangles = nbelements(TheMesh)
    TheMesh.element_measures = zeros(RealNum,nbtriangles)
    #
    # Init. various variables to avoid memory reallocation in the loop
    # 
    TheTriangle = zeros(Index,3)
    x1 = zeros(RealNum,2)
    x2 = zeros(RealNum,2)
    x3 = zeros(RealNum,2) # Nodes of a triangle
    #edge_
    # Loop over triangles and build the elementary matrix
    #
    for iloop in 1:nbtriangles
        #
        # Get the coordinates of the nodes
        #
        TheTriangle = TheMesh.elements[:,iloop]
        x1 = TheMesh.nodes[:, TheTriangle[1]]
        x2 = TheMesh.nodes[:, TheTriangle[2]]
        x3 = TheMesh.nodes[:, TheTriangle[3]]
        #
        # Compute the surface of the triangle
        #
        TheMesh.element_measures[iloop] = 0.5*((x2[1] - x1[1])*(x3[2]- x1[2]) - (x3[1] - x1[1])*(x2[2] - x1[2]))   # This is |T|
        #
    end #For iloop (the loop over the triangles is finished)
    #
    # Set the associated flag to true
    #

    TheMesh.elm_meas_flag = true
    
end # of element_surfaces

#
# Compute the volumes of the tetrahaedra of a 3-D mesh
#
function element_volumes(TheMesh::volcmesh)::Bool
    print("Not available at the moment")
    return(false)
end

###############################
# Boundary surfaces and lengths
###############################


#
# Compute the lengths of the edges  of a 2-D mesh
#
function edge_lengths(TheMesh::volcmesh)::Bool
    # This routine computes the lengths of the edge segments in a 2D mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output: no
    # - lengths_vector : array containing the lengths of the boundary segments in the mesh
    # stored in the mesh data array edge_measures
    # Each column of the edges array contains the n° of begin and end nodes and the n° of the boundary part
    # it belongs to
    #
    # Read the number of boundary segments and initialise the length array
    #
    # 

    # Compute the differences between the start and end points of each boundary segment

    Xdiff = TheMesh.nodes[:,TheMesh.edges[2,:] ] - TheMesh.nodes[:,TheMesh.edges[1,:] ]
    
    # Compute the euclidean norm of each difference vector

    Ne = nbedges(TheMesh)   # also equal to `size(Xdiff,2)`
    lengths_vector = zeros(RealNum,Ne)  # build a column vector with Ne elements

    for iloop = 1:Ne # Euclidean norm (package LinearAlgebra) of each column ie along the first index ("à l'ancienne")
        lengths_vector[iloop] = LinearAlgebra.norm(Xdiff[:,iloop])
    end
    
    TheMesh.edge_measures = lengths_vector
    TheMesh.edg_meas_flag = (Ne > 0)    # if the edge array is empty, there are no boundary measures
    return(TheMesh.edg_meas_flag)
end # End of edge_lengths

#
# Compute the surface of the edges of a 3-D  tetrahedral mesh
#
function edge_surfaces(TheMesh)::Bool
    # This routine computes the surfaces of the edge triangles in a 3D mesh.
    # Input variables:
    # - TheMesh : mesh of class volcmesh
    # Output:
    # - surfaces : array containing the surfaces of the boundary triangles in the mesh
    # stored in the mesh data array edge_measures
    # Each column of the edges array contains XXXXXXXXXXXXXXX and the n° of the boundary part
    # it belongs to
    #
    # Read the number of boundary segments and initialise the length array
    #
    # Loop over edge segments and compute their length

    surfaces = np.array([])

    TheMesh.edge_measures = surfaces
    TheMesh.edg_meas_flag = false

    return(TheMesh.edg_meas_flag)

end # of edge_surfaces
#


############## BOUNDARY NORMALS
#############
######## WARNING : does not modify a volcmesh structure for the time being
######## This should change soon


function mk_edge_normals(TheMesh::volcmesh)::RealMat
    #
    # Check if the edge measures are present
    #
    if !TheMesh.edg_meas_flag
        mk_edge_measures(TheMesh)
    end

    TheEdge = zeros(Index, (3,1))   # this has nothing to do with U2 of course.

    # Normals = zeros(RealNum,(2,nbedges(TheMesh)) ###useless created and filled up later on

   
    # compute "tangential" vectors

    Xdiff = TheMesh.nodes[:,TheMesh.edges[2,:] ] - TheMesh.nodes[:,TheMesh.edges[1,:] ]

    # Normalize using the edge lengths
    # TheMesh.edge_measures is a column vector, so we transpose it
    # in order to divide each column of Xdiff by its length (2-norm)

    Xdiff = Xdiff ./ TheMesh.edge_measures'

    # Get normal vectors [x y] --> [y -x]

    Normals = [Xdiff[2,:]' ; -Xdiff[1,:]']

    return(Normals)

end



###########################################################################
###########################################################################
#########                         IO METHODS                ###############
###########################################################################
###########################################################################


function MatlabMeshImport(FileName::String)
    # Imports points, edges and triangles from a given .mat file
    # these arrays typically are built with the Matlab PDE Toolbox (c)
    # Therefore is is assumed that the arrays in the file are named
    # p,e,t (standing for points, edges,triangles)
    # First, check the existence of the file

    if isfile(FileName)
        values = MAT.matread(FileName) # The Values should be p, e and t
        if haskey(values,"e") && haskey(values,"p") && haskey(values,"t")
            points = get(values,"p",@EmptyRealMat)
            e = get(values,"e",@EmptyIntMat)
            t = get(values,"t",@EmptyIntMat)
            #
            # In 2-D, keep the lines number 1, 2 and 4 (Start node, End node, boundary number)
            # t should be unchanged since we keep the region number for later uses if it is present in the array
            #
            Dimension = size(points,1)
            if (Dimension == 2)
                edges = convert(IndexMat, [ e[1:2,:] ; e[5,:]' ]) # 5 ou 5 ???
                triangles = convert(IndexMat,t)
            elseif (Dimension ==3)
                # do nothing for the moment
                print("3-D not supported yet")
            end
        else
            points = @EmptyRealMat
            edges = @EmptyIntMat
            triangles = @EmptyIntMat
        end
    else
        points = @EmptyRealMat
        edges = @EmptyIntMat
        triangles = @EmptyIntMat
    end 

    return(points,edges,triangles)

end

###########################################################################
###########################################################################
#########                       2D Utilities                ###############
###########################################################################
###########################################################################

function meshzeros(TheMesh::volcmesh, ndof::Index = 1)::RealVec
    # Creates an null array with as many elements as nodes in the mesh.
    # Input Parameters : 
    # - TheMesh : volcmesh type mesh object
    # - [ndof] : number of degrees of freedom at each node [default = 1]
    # Output parameters
    # - an array with as many elements as there are nodes in the mesh
    #
    NbNodes = nbnodes(TheMesh)
    #
    Vec = zeros(RealNum, NbNodes*ndof)
    #
    return(Vec)
end # of meshzeros

function meshvector(TheMesh::volcmesh, funk, ndof::Index = 1)::RealVec
    # Creates an array which contains the value of a function evaluated at the nodes of a mesh.
    # Input Parameters : 
    # - TheMesh : volcmesh type mesh object
    # - [funk] : the function to be evaluated at the nodes [default fzero]
    #   if ndof is more than 1, funk has to return a RealVec of length ndof
    # - [ndof] : number of degrees of freedom at each node [default = 1]
    # Warning : if ndof > 1, funk should return a array with ndof*nbnodes elements.
    # If no function is provided, an array of zeros is created. In this case the user vould better call
    # the meshzeros function.
    # Output parameters
    # - an array with as many elements as there are nodes in the mesh
    #
    NbNodes = nbnodes(TheMesh)
    #
    Vec = zeros(RealNum, NbNodes*ndof)
    #
    index = 1
    for numnode = 1:NbNodes
        Vec[index:(index+ndof-1)] .= funk(TheMesh.nodes[:,numnode]) # use broadcasting in case ndof = 1
        index += ndof
    end
    #

    return(Vec)
end # of meshvector

function meshdivergencevector(TheMesh::volcmesh, divergence_func, ndof::Index = 1)::RealVec
    # Creates an array which contains the divergence of a tensor function evaluated at the nodes of a mesh.
    NbNodes = nbnodes(TheMesh)

    # Initialize the vector to store the divergence values
    Vec = zeros(RealNum, NbNodes * ndof)

    index = 1
    for numnode = 1:NbNodes
        # Evaluate the tensor at the current node
        tensor = divergence_func(TheMesh.nodes[:, numnode])

        # println(divergence_func(TheMesh.nodes[:, numnode]), "  ", numnode)

        # Store the divergence in the output vector 
        Vec[index:(index + ndof - 1)] .= tensor
        # println("Vecc", Vec[index], "  ",index, "  ", TheMesh.nodes[:, numnode])
        # println("Vecc",  Vec[index + ndof - 1], " ", index + ndof - 1, "  ", TheMesh.nodes[:, numnode])
        index += ndof
    end

    # return Vec
    return Vec
end # of meshdivergence



function intriangle(X::RealVec,TheMesh::volcmesh)
    # Checks if a given point is inside or on the boundary of one of the triangles of a mesh
    # If the triangle areas in the mesh are not computed, this will be done.
    #
    # INPUT VARIABLES
    # - X : a 2-D coordinate vector
    # - TheMesh : volcmesh type object (non empty)
    #
    # OUTPUT 
    # - TriangleZzz : number(s) of the triangle where the point X is (array of Index), 0 if the point is outside the mesh
    # If  one triangle is found, a single value of type Index is returned (the number of the triangle in the mesh),
    # if several have been found, an array is returned
    # else 0 is returned.
    # 
    # 
    # The method used consists in computing the normalized barycentric coordinates of the points X. They are given by
    #
    # alpha = S1 / T, beta = S2 / T and gamma = 1 - alpha - beta
    #
    # where S1 is the (signed) surface of the triangle (X T_2 T_3), S2 the (signed) surface of (T_1 X T_3)
    # and T_i is the vertex n°i of the triangle
    # WARNING : the vertices are supposed (as in all the routines) to be stored in trigonnometric order. 
    #
    # WARNING : the 0 value is returned in two cases : the point is is 0 triangle or more than 1 (meaning it is on a boundary)
    #
    # Initlalisations
    #
    # First check if the triangle areas have been computed

    if !(TheMesh.elm_meas_flag)
        mk_element_measures(TheMesh)
    end

    # Initialise the output array to empty

    TriangleZzz::IndexVec = []

    nbtriangles = nbelements(TheMesh)

    # Loop over the elements

    for iloop in 1:nbtriangles
    
        # Get the coordinates of the nodes

        TheTriangle = TheMesh.elements[:,iloop]
        x1 = TheMesh.nodes[:, TheTriangle[1]]
        x2 = TheMesh.nodes[:, TheTriangle[2]]
        x3 = TheMesh.nodes[:, TheTriangle[3]]

        # Compute the surfaces of the sub triangles

        S1 = 0.5*( (X[1] - x3[1])*(x2[2] - x3[2]) - (X[2] - x3[2])*(x2[1] - x3[1]) ) # surface of (P x2 x3)
        S2 = 0.5*((X[2] - x3[2])*(x1[1] - x3[1]) - (X[1] - x3[1])*(x1[2] - x3[2])) # surface of (x1 P x3)

        # Compute the barycentric coordinates

        alpha = S1 / TheMesh.element_measures[iloop]
        beta  = S2 / TheMesh.element_measures[iloop]
        gamma = 1 - alpha - beta
        
        #
        # Now check the signs of the three coordinates
        #   
        epsilon = 5.0e-10
        if (alpha >= -epsilon) && (beta >= -epsilon) && (gamma >= -epsilon) # Allow boundaries
            append!(TriangleZzz, iloop) # Point is inside the triangle
        end
    end

    # Build a convenient output (therefore no return type for the function)
    # The boolean flag is set to true if multiple triangles were found
    if (length(TriangleZzz) == 0)
        return(0)
    else
        return(TriangleZzz)
    end

    


end # of intriangle



function insidetriangle(X::RealVec,TheMesh::volcmesh)
    # Checks if a given point is strictly inside one of the triangles of a mesh
    # If the triangle areas in the mesh are not computed, this will be done.
    #
    # INPUT VARIABLES
    # - X : a 2-D coordinate vector
    # - TheMesh : volcmesh type object (non empty)
    #
    # OUTPUT 
    # - TriangleZzz : number(s) of the triangle where the point X is (array of Index), 0 if the point is outside the mesh
    # If  one triangle is found, a single value of type Index is returned (the number of the triangle in the mesh),
    # else 0 is returned.
    # 
    # 
    # The method used consists in computing the normalized barycentric coordinates of the points X. They are given by
    #
    # alpha = S1 / T, beta = S2 / T and gamma = 1 - alpha - beta
    #
    # where S1 is the (signed) surface of the triangle (X T_2 T_3), S2 the (signed) surface of (T_1 X T_3)
    # and T_i is the vertex n°i of the triangle
    # WARNING : the vertices are supposed (as in all the routines) to be stored in trigonnometric order. 
    #
    # WARNING : the 0 value is returned in two cases : the point is is 0 triangle or more than 1 (meaning it is on a boundary)
    #
    # Initlalisations
    #
    # First check if the triangle areas have been computed

    if !(TheMesh.elm_meas_flag)
        mk_element_measures(TheMesh)
    end

    # Initialise the output array to empty

    TriangleZzz::IndexVec = []

    nbtriangles = nbelements(TheMesh)

    # Loop over the elements

    for iloop in 1:nbtriangles
    
        # Get the coordinates of the nodes

        TheTriangle = TheMesh.elements[:,iloop]
        x1 = TheMesh.nodes[:, TheTriangle[1]]
        x2 = TheMesh.nodes[:, TheTriangle[2]]
        x3 = TheMesh.nodes[:, TheTriangle[3]]

        # Compute the surfaces of the sub triangles

        S1 = 0.5*( (X[1] - x3[1])*(x2[2] - x3[2]) - (X[2] - x3[2])*(x2[1] - x3[1]) ) # surface of (P x2 x3)
        S2 = 0.5*((X[2] - x3[2])*(x1[1] - x3[1]) - (X[1] - x3[1])*(x1[2] - x3[2])) # surface of (x1 P x3)

        # Compute the barycentric coordinates

        alpha = S1 / TheMesh.element_measures[iloop]
        beta  = S2 / TheMesh.element_measures[iloop]
        gamma = 1 - alpha - beta
        
        #
        # Now check the signs of the three coordinates
        #

        if (alpha > 0.0) && (beta > 0.0) && (gamma > 0.0)  # If boundaries of triangles were valid >= would replace >
            append!(TriangleZzz, iloop)
        end
    end

    # Build a convenient output (therefore no return type for the function)
    # The boolean flag is set to true if multiple triangles were found

    if (length(TriangleZzz) == 0)
        return(0)
    else
        return(TriangleZzz[1])
    end


end # of insidetriangle

function common_nodes(Elem1::Index,Elem2::Index,TheMesh::volcmesh)::IndexVec
    #
    # common_nodes(Elem1::Index,Elem2::Index,TheMesh::volcmesh)::IndexVec
    #
    # for two given elements, this function return the list of their common nodes
    #
    # Input
    # - Elem1 : number of the first element
    # - Elem2 : number of the second element
    # - TheMesh : mesh to be used
    #
    # Output
    # - an array of integers containing the numbers of the common nodes in the mesh
    # (empty if non common nodes)
    #
    # In 2D, if there are two common nodes then the triangular have one common side (thank you captain obvious).

    # First check that the elements really exist

    Nelem = nbelements(TheMesh)
    if Nelem == 0
        throw(ArgumentError("The provided mesh is empty"))
    end

    if (Elem1 < 1) || (Elem1 > Nelem)
        throw(ArgumentError("Elem1 is out of range"))
    end

    if (Elem2 < 1) || (Elem2 > Nelem)
        throw(ArgumentError("Elem2 is out of range"))
    end

    # Extract the node numbers for each element

    t1 = TheMesh.elements[1:3,Elem1]
    t2 = TheMesh.elements[1:3,Elem2]

    # Use function findall to get common numbers

    return(t1[findall(in(t2),t1)])

end

###########################################################################
###########################################################################
###########################################################################
# The fracture & reservoir class in 2D or 3D
###########################################################################
###########################################################################
###########################################################################

# 1D objects in a 2D environment.

###################################
#### The line class
###################################

### Linefrac2D : should be embedded in a 2D mesh
#
# a line consists in a set of (sorted) points making a set of segments
#
# the structure features a connectivity table (elements) to override problems
# in the unsorted case.
# a variable edges which will point to the nodes array since the edges
# of a line are its defining points
# As in JuliaLang the = operation makes two arrays be the same this is OK
# Edge measures will also be in the data but always void ans the flag edg_meas_flag
# will be set to True by default


mutable struct linefrac2D
    #
    # Main data (the mesh is stored in (points, edges,elements) format)
    #
    nodes::RealMat
    edges::IndexVec
    elements::IndexMat
    #
    # Reference to the containing mesh
    embed_mesh::volcmesh # The mesh itself
    embed_triangles::IndexMat # for each node indicates the triangles it belongs to
    #
    # Secondary Data (computed by called methods)
    #
    element_measures::RealVec
    element_normals::RealMat
    # Flags
    #
    isempty::Bool
    elm_meas_flag::Bool
    elm_norm_flag::Bool
    #
    # Extra flag for the reservoir case
    #
    isreservoir::Bool
    #
    # Basic constructor (returns void structure)

    linefrac2D() = new(@EmptyRealMat,@EmptyIntVec,@EmptyIntMat,volcmesh(), @EmptyIntMat,@EmptyRealVec,@EmptyRealMat,true,false,false,false)

    # Advanced Constructor #1 : sorted points

    function linefrac2D(points::RealMat,EmbeddingMesh::volcmesh)
        # Building a line fracture from its mesh points
        # The fracture will be given by all the line segments connecting
        # each point to its successor in the array
        # The points should be sorted and all be inside the mesh EmbededMesh
        #
        # WARNING : a line fracture is not a closed curve. Closed curves will be defined
        # ine the same framewrok via the reservoir2D creation function, setting the "isreservoir"
        # flag to true if this is the case
        # Inputs
        #
        # points : array of the points defining the edges of the segments
        # EmbeddingMesh : mesh in which all the points are to be. (of volcmesh type)
        #
        # Output : a linefrac structure with segments measures non computed
        #

        # INIT : testing the confirmity of the input data

        
        # check if 2D points
        if (size(points,1) != 2)
            throw(ArgumentError("Points should be 2D"))
        end

        # Now test if the first & last points coincide

        if isapprox(points[:,begin],points[:,end])
            throw(ArgumentError("Points make a closed Loop"))
        end

        # Check if the points are all in the provided mesh and store the
        # numbers of the triangles

        NumPoints = size(points,2)
   
   
        Triangles = zeros(Index,(NumPoints,1))
        for iloop in 1:NumPoints
            localtoto = intriangle(points[:,iloop],EmbeddingMesh)
            if (localtoto == 0) # the current point is not in the mesh
                throw(ArgumentError("Points must lie in the mesh"))
            else # the current point lies in at least one element
                Triangles[iloop] = localtoto[1]
            end
        end
    
        # Here we arrive with no error (no exception thrown)
        # Proceed with effective construction of the object

        # Build a useless connectivity table
        # (actually will be used for the computation of the element lengths
        # without switching cases)

        ElemZzz = zeros(Index,(2,NumPoints-1))
        for iloop in 1:(NumPoints-1)
            ElemZzz[:,iloop] = [iloop ; iloop+1]
        end

        Edges = collect(1:NumPoints)
    

        # Create a temporary structure with all the available data

        return(new(points,Edges,ElemZzz,EmbeddingMesh,Triangles,@EmptyRealVec,@EmptyRealMat,false,false,false,false))

    end

    # Advanced Constructor #2 : points and connectivity table

    function linefrac2D(points::RealMat,connect_table::IndexMat, EmbeddingMesh::volcmesh)
        # Building a line fracture from its mesh points
        # The fracture will be given by all the line segments connecting
        # each point to its successor in the array
        # The points should be not sorted and all be inside the mesh EmbededMesh
        #
        # WARNING : a line fracture is not a closed curve. Closed curves will be defined
        # as another object celles Reservoir2D
        # Inputs
        #
        # points : array of the points defining the edges of the segments
        #
        # EmbeddingMesh : mesh in which all the points are to be. (of volcmesh type)
        #
        # Output : a linefrac structure with segments measures non computed
        #

        # INIT : testing the confirmity of the input data

        
        # check if 2D points
        if (size(points,1) != 2)
            throw(ArgumentError("Points should be 2D"))
        end

        # Now test if we have a closed loop (brute force : is the same point twice in the list)

        NumPoints = size(points,2) # total number of points

        for iloop in 1:(NumPoints-1)
            for jloop in (iloop+1):NumPoints
                if isapprox(points[:,iloop], points[:,jloop])
                    throw(ArgumentError("Points make a closed Loop"))
                end
            end
        end

        # Check if the points are all in the provided mesh and store the
        # numbers of the triangles

   
        Triangles = zeros(Index,(NumPoints,1))
        for iloop in 1:NumPoints
            localtoto = intriangle(points[:,iloop],EmbeddingMesh)
            if (localtoto == 0) # the current point is not in the mesh
                throw(ArgumentError("Points must lie in the mesh"))
            else # the current point lies in at least one element
                Triangles[iloop] = localtoto[1]
            end
        end
    
        # Here we arrive with no error (no exception thrown)
        # Proceed with effective construction of the object

        # Build a useless connectivity table
        # (actually will be used for the computation of the element lengths
        # without switching cases)

        ElemZzz = connect_table
        Edges = collect(1:NumPoints)

        # Create a temporary structure with all the available data

        return(new(points,Edges,ElemZzz,EmbeddingMesh,Triangles,@EmptyRealVec,@EmptyRealMat,false,false,false,false))

    end
end

######## Associated methods

# number of points

function nbnodes(TheLine::linefrac2D)::Index
    # Return the number of nodes in the provided linefrac2D object
    return(size(TheLine.nodes,2))
end

function nbedges(TheLine::linefrac2D)::Index
    # Return the number of boundary nodes in the provided linefrac2D object
    # i.e. the number of nodes
    return(length(TheLine.edges))
end

function nbelements(TheLine::linefrac2D)::Index
    # Return the number of elements in the provided linefrac2D object
    return(size(TheLine.elements,2))
end

function dimension(TheLine::linefrac2D)::Index
    # return the space dimension of the embedding mesh for a linefrac2D object
    return(size(TheLine.nodes,1))
end

# Computing lengths of the segments

function mk_segment_lengths(TheLine::linefrac2D)::Bool
    # This routine computes the lengths of the segments of a linefrac2D object.
    # Input variables:
    # - TheLine : linefrac2D object
    # Output: Bool set to false if the object is empty
    

    # check if there are actually elements in the line

    Ne = nbelements(TheLine)

    if (Ne > 0)

        # Compute the differences between the start and end points of each segment
        # Each segment is givenby the connectivity table elements

        Xdiff = TheLine.nodes[:,TheLine.elements[2,:] ] - TheLine.nodes[:,TheLine.elements[1,:] ]
        
        # Compute the euclidean norm of each difference vector

        lengths_vector = zeros(RealNum,Ne)  # build a column vector with Ne elements

        for iloop = 1:Ne # Euclidean norm (package LinearAlgebra) of each column ie along the first index ("à l'ancienne")
            lengths_vector[iloop] = LinearAlgebra.norm(Xdiff[:,iloop])
        end
        
        TheLine.element_measures = lengths_vector
    end

    TheLine.elm_meas_flag = (Ne > 0)    # if the edge array is empty, there are no boundary measures
    return(TheLine.elm_meas_flag)
end # End of edge_lengths

function mk_segment_normals(TheLine::linefrac2D)::Bool
    # This routine builds the unitary normals of each segment of a linefrac2D object.
    # Input variables:
    # - Line : linefrac2D object
    # Output: bool indicating the success of the routine
    #
    # if TheLine is empty no computation will be done
    # if the segment lengths are not computed, thine method computes them and stores them in the data
    #
    # HENCE : this routine actually computes segment lenghts AND normals

    # Read the number of segments

    Ne = nbelements(TheLine)

    if (Ne > 0)

        # Compute the differences between the start and end points of each segment
        # Each segment is givenby the connectivity table elements

        Xdiff = TheLine.nodes[:,TheLine.elements[2,:] ] - TheLine.nodes[:,TheLine.elements[1,:] ]

        # First check if the lengths need to be computed

        if !TheLine.elm_meas_flag
                
            
            # Compute the euclidean norm of each difference vector

            lengths_vector = zeros(RealNum,Ne)  # build a line vector with Ne elements

            for iloop = 1:Ne # Euclidean norm (package LinearAlgebra) of each column ie along the first index ("à l'ancienne")
                lengths_vector[iloop] = LinearAlgebra.norm(Xdiff[:,iloop])
            end
            
            TheLine.element_measures = lengths_vector
            TheLine.elm_meas_flag = true
        end

        # Now compute the normals : Note that Xdiff contains actually the tangent vector to each segment
        #
        # First swap lines in Xdiff and normalize
        # Second multiply the firstsecond row by -1

        TheLine.element_normals = Xdiff[[2,1],:] ./ TheLine.element_measures'
        TheLine.element_normals[2,:] = - TheLine.element_normals[2,:]
        
    end

    TheLine.elm_norm_flag = (Ne > 0)    # if the edge array is empty, there are no boundary measures
    return(TheLine.elm_norm_flag)
end # End of edge_lengths


# Remeshing a line to avoid too distant vertices

function local_remesh_2Dline(TheLine::linefrac2D)
    # Given the index PointNum of a point in the object add (by dichotomy) one point
    # in the segment number ElemNum as a start point
    #
    # The connectivity table "elements" will be used to cover the case when the points are unsorted.
    #
    # The routine works as follows : for each segment the triangles are checked :
    # - if the two triangles are the same, nothing is done
    # - if the two triangles are neighbours, nothing is done
    # - if the two triangles are not neighbours a point is added at the center of the segment. 
    #
    # In the last case if the point does not solve the problem i.e. two of the nodes are still
    # too far from each other, repeat the process

    Ne = nbelements(TheLine)

    if (Ne == 0)
        throw(ArgumentError("LineFrac2D object is empty"))
    end

    index = 1
    while (index <= Ne)  # scan over all the segments
                        # since the number of elements may change during the process
                        # a while loop is necessary

        CurrentElem = TheLine.elements[:,index]

        # Start and end of the current segment
        p1 = TheLine.nodes[:,CurrentElem[1]] # First point
        p2 = TheLine.nodes[:,CurrentElem[2]] # Last point

        # Examine the common nodes of the two points


        # Get the triangles of the mesh containing each point

        Tr1 = TheLine.embed_triangles[CurrentElem[1]]
        Tr2 = TheLine.embed_triangles[CurrentElem[2]]

        # Get the common nodes of the two triangles

        CN = common_nodes(Tr1,Tr2,TheLine.embed_mesh)

        if (length(CN) < 2)  # Refine until the two ends of the segment are close enough

            # Seek a point in a neighbouring or the same element via dichotomy

            pnew = copy(p2)
            Trnew = Tr2     # to avoid the locality of Trnew in the next block
            while (length(CN) < 2)
                pend = copy(pnew)
                pnew = (p1 + pend) * 0.5 # Insert the center of the segment in the linefrac2D object
                Triank = intriangle(pnew, TheLine.embed_mesh)
                Trnew = Triank[1]
                CN = common_nodes(Tr1,Trnew,TheLine.embed_mesh)
            end

            # Add the new point and modify the connectivity table

            TheLine.nodes = hcat(TheLine.nodes , pnew)
            Ne = +(Ne , 1) # we have Ne segments thus Ne+1 nodes

            TheLine.elements[2,index] = Ne+1    # The new end of the segment is the last 
                                                # point added to the structure
            TheLine.elements = hcat(TheLine.elements[:,1:index], [Ne+1;CurrentElem[2]],
                                    TheLine.elements[:,index+1:end])

            TheLine.embed_triangles = vcat(TheLine.embed_triangles, Trnew)


        end

        index = +(index,1) # Move to the next segment
    end

    # (Re)compute the segment lengths and if necessary the new normals

    mk_segment_lengths(TheLine)
    if TheLine.elm_norm_flag
        mk_segment_normals(TheLine)
    end
    
end


# Computing the intersections of the segments of a linefrac2D object with the edges of
# the triangles of its embedding mesh
#
# WARNING : we suppose here that the fracture has been refined so that each segment
# intersects an edge of the mesh
#
# The routine will also compute the effective lengths the segments have inside each triangle


function frac_mesh_intersect(TheLine::linefrac2D)::Tuple{Vector{Bool},RealMat,RealMat}
    # This functions computes for each segment of the TheLine, the points where they
    # intersect the common edge of their triangles in the embedding mesh
    # The line is suposed to be properly meshed in the way (see local_remesh_2Dline)
    #
    # Input(s) :
    # - TheLine : a linefrac2D object
    #
    # Output(s) :
    # - HasIntersection : array of booleans. If an element is set to true then the corresponding segment lies across
    #                                        two triangles. If not the segment lies in one element only
    # - InterSectPoints : float array of 2 x (Number of segments). Coordinates of the intersection points
    # - Seglengths :  float array of 2 x (Number of segments). Respective lengths in each triangle.
    #

    # For each segment, three cases will be considered :
    # - there is a full side as a neighbour --> compute the intersection
    # - there is one point in common --> the intersection is a vertex of the mesh
    # - the segment lies inside one triangle

    # Init.

    Ne = nbelements(TheLine)

    if (Ne == 0)
        throw(ArgumentError("LineFrac2D object is empty"))
    end

    HasIntersection = Vector{Bool}(undef, Ne)
    InterSectPoints = zeros(RealNum, (2,Ne))
    Seglengths = zeros(RealNum, (2,Ne))

    # Check if the degment lengths have been computed

    if !TheLine.elm_meas_flag
        mk_segment_lengths(TheLine)
    end

    # Loop over segments

    for index = 1:Ne

        CurrentElem = TheLine.elements[:,index]

        # Start and end of the current segment
        p1 = TheLine.nodes[:,CurrentElem[1]] # First point
        p2 = TheLine.nodes[:,CurrentElem[2]] # Last point

        # Get the triangles of the mesh containing each point

        Tr1 = TheLine.embed_triangles[CurrentElem[1]]
        Tr2 = TheLine.embed_triangles[CurrentElem[2]]

        # Get the common nodes of the two triangles

        CN = common_nodes(Tr1,Tr2,TheLine.embed_mesh)

        # Check the 3 different cases

        if (length(CN) == 3) # The segment lies in one triangle
            
            HasIntersection[index] = false
            Seglengths[1,index] = TheLine.element_measures[index]
            Seglengths[2,index] = TheLine.element_measures[index]

        elseif (length(CN) ==1) # The intersection is a vertex of the mesh
           
            HasIntersection[index] = true
            InterSectPoints[:,index] = TheLine.embed_mesh.nodes[:,CN[1]]
            Seglengths[:,index] = 
                    [LinearAlgebra.norm(InterSectPoints[:,index]-p1) ; LinearAlgebra.norm(p2 - InterSectPoints[:,index])]

        else # a side of one triangle is the intersection
           
            HasIntersection[index] = true
            
            # Compute the intersection point
            S1 = TheLine.embed_mesh.nodes[:,CN[1]]
            S2 = TheLine.embed_mesh.nodes[:,CN[2]]


            d = p2 - p1
            dprime = S2 - S1
            Y = [(S1[1] - p1[1]) ; (S1[2] - p1[2])] # RHS of the system
            Delta = dprime[1]*d[2] - d[1]*dprime[2] # Determinant of the (Cramer) system 

            t = (dprime[1]*Y[2] - dprime[2]*Y[1]) / Delta # Step to "walk from p1" to the intersection

            InterSectPoints[:,index] = p1 + t*d
            Seglengths[:,index] = 
            [LinearAlgebra.norm(InterSectPoints[:,index]-p1) ; LinearAlgebra.norm(p2 - InterSectPoints[:,index])]

        end
    end

    return((HasIntersection,InterSectPoints,Seglengths))
end


#######################################################################
#
#  The reservoir case : last point == first point
#
#######################################################################

# A reservoir will be defined by a set of points like a fracture. The first 
# point should not be repeated at the end of the list and will be added by the creation
# routine.

function reservoir2D(points::RealMat,EmbeddingMesh::volcmesh)::linefrac2D

    # First create a fracture object
    # Remember the last point in the list must not be the first

    TheRes = linefrac2D(points,EmbeddingMesh)

    # read the connectivity table to get the first and last node numbers

    numfirstnode = TheRes.elements[1,1]
    numlastnode = TheRes.elements[2,end]
    nbpoints = nbnodes(TheRes)


    TheRes.nodes = hcat(TheRes.nodes,TheRes.nodes[:,numfirstnode])  # add first point at the end
    TheRes.elements = hcat(TheRes.elements, [numlastnode ; nbpoints+1])

    TheRes.embed_triangles = vcat(TheRes.embed_triangles,TheRes.embed_triangles[1])


    TheRes.isreservoir = true

    return(TheRes)

end

function reservoir2D(points::RealMat,connect_table::IndexMat,EmbeddingMesh::volcmesh)::linefrac2D

    # Create a reservoir from a list of points and a connectivity table


    # First create a fracture object
    # Remember the last point in the list must not be the first

    TheRes = linefrac2D(points,connect_table,EmbeddingMesh)

    # read the connectivity table to get the first and last node numbers

    numfirstnode = TheRes.elements[1,1]
    numlastnode = TheRes.elements[2,end]
    nbpoints = nbnodes(TheRes)

    TheRes.nodes = hcat(TheRes.nodes,TheRes.nodes[:,numfirstnode])  # add first point at the end
    TheRes.elements = hcat(TheRes.elements, [numlastnode ; nbpoints+1])

    TheRes.embed_triangles = vcat(TheRes.embed_triangles,TheRes.embed_triangles[1])


    TheRes.isreservoir = true

    return(TheRes)

end


function reservoir2D(TheLine::linefrac2D)

    # Transform a fracture into a reservoir

    # read the connectivity table to get the first and last node numbers

    numfirstnode = TheLine.elements[1,1]
    numlastnode = TheLine.elements[2,end]
    nbpoints = nbnodes(TheLine)

    TheLine.nodes = hcat(TheLine.nodes,TheLine.nodes[:,numfirstnode])  # add first point at the end
    TheLine.elements = hcat(TheLine.elements, [numlastnode ; nbpoints+1])

    TheLine.embed_triangles = vcat(TheLine.embed_triangles,TheLine.embed_triangles[1])


    TheLine.isreservoir = true

end

function linevector(TheLine::linefrac2D, funk, ndof::Index = 1)::RealVec
    # Creates an array which contains the value of a function evaluated at the nodes of the line segments.
    # Input Parameters: 
    # - TheLine : linefrac2D type line object
    # - funk : the function to be evaluated at the nodes
    # - ndof : number of degrees of freedom at each node [default = 1]
    # Output parameters:
    # - an array with as many elements as there are nodes in the line

    NbNodes = nbelements(TheLine)# Get the number of nodes
    
    Vec = zeros(RealNum, NbNodes * ndof)  # Initialize the result vector

    index = 1
    for node in 1:NbNodes
        # Evaluate the function at the current node
        Vec[index:(index + ndof - 1)] .= funk(TheLine.nodes[:, node])
        index += ndof
    end
  
    return Vec

end

function midlinevector(TheLine::linefrac2D, funk, ndof::Index = 1)::RealVec
    # Creates an array which contains the value of a function evaluated at the midpoints of the line segments.
    # Input Parameters: 
    # - TheLine : linefrac2D type line object
    # - funk : the function to be evaluated at the midpoints
    # - ndof : number of degrees of freedom at each node [default = 1]
    # Output parameters:
    # - an array with as many elements as there are midpoints in the line

    NbElements = nbelements(TheLine)  # Get the number of elements
    NbMidpoints = NbElements  # One midpoint per element
    
    Vec = zeros(RealNum, NbMidpoints * ndof)  # Initialize the result vector

    midpoints = zeros(2, NbElements)
    Vec = zeros(2*NbElements)
    for iloop = 1:NbElements

        # Calculate the midpoint for this segment
        midpoints[:, iloop] = (TheLine.nodes[:, TheLine.elements[1, iloop]] + TheLine.nodes[:, TheLine.elements[2, iloop]]) / 2
        
        # Evaluate the function at the midpoint
        funk_value = funk(midpoints[:, iloop]) 

        # Multiply the function value by the normal
        Vec[(iloop - 1) * ndof + 1:(iloop - 1) * ndof + ndof] .=  TheLine.element_normals[:, iloop]
        # println(Vec[(iloop - 1) * ndof + 1:(iloop - 1) * ndof + ndof])
    end

    return Vec
end


##################################################################################################
##################################################################################################
##########################
##########################  Graphic utilities
##########################
########################## (using the Plots package)
##########################
##################################################################################################
##################################################################################################
using Plots

function plot2Dmeshsave(TheMesh::volcmesh)
    #
    # Plot all the triangles of a 2D mesh
    #
    Nelem = nbelements(TheMesh)
    if Nelem == 0
        throw(ArgumentError("The provided mesh is empty"))
    end

    # if we really have a mesh proceed to plot
    NbElem = nbelements(TheMesh)

    # First element, we add the first point at the end to plot closed triangles
    X = zeros(RealNum, 2, 4)
    TheElem = TheMesh.elements[1:3, 1]
    X = [TheMesh.nodes[:, TheElem[1]] TheMesh.nodes[:, TheElem[2]] TheMesh.nodes[:, TheElem[3]] TheMesh.nodes[:, TheElem[1]]]
    plot(X[1, :], X[2, :], lc=:blue, legend=false)

    # Other elements
    for iloop in 2:NbElem
        TheElem = TheMesh.elements[1:3, iloop]
        X = [TheMesh.nodes[:, TheElem[1]] TheMesh.nodes[:, TheElem[2]] TheMesh.nodes[:, TheElem[3]] TheMesh.nodes[:, TheElem[1]]]
        plot!(X[1, :], X[2, :], lc=:blue, legend=false)
    end

    # Display the nodes
    scatter!(TheMesh.nodes[1, begin:end], TheMesh.nodes[2, begin:end], mc=:red, markersize=1.7, legend=false)

    # Save the plot as an image (e.g., PNG file)
    savefig("mesh_plot.png")  # You can change the file name and format (e.g., ".pdf", ".jpg", etc.)

    println("Plot saved as 'mesh_plot.png'")
end



function plot2Dmesh(TheMesh::volcmesh)
    #
    # Plot all the triangles of a 2D mesh
    #
    Nelem = nbelements(TheMesh)
    if Nelem == 0
        throw(ArgumentError("The provided mesh is empty"))
    end

    # if we really have a mesh procees to plot

    NbElem = nbelements(TheMesh)

    # First element, we add the firstpoint at the end in order to plot closed triangles

    X = zeros(RealNum,2,4)
    TheElem = TheMesh.elements[1:3,1]
    X = [ TheMesh.nodes[:,TheElem[1]]  TheMesh.nodes[:,TheElem[2]]  TheMesh.nodes[:,TheElem[3]] TheMesh.nodes[:,TheElem[1]]]
    plot(X[1,:],X[2,:],lc=:blue,legend=false)


    # Other elements

    for iloop in 2:NbElem # always make it simple in the first place
        TheElem = TheMesh.elements[1:3,iloop]
        X = [ TheMesh.nodes[:,TheElem[1]]  TheMesh.nodes[:,TheElem[2]]  TheMesh.nodes[:,TheElem[3]] TheMesh.nodes[:,TheElem[1]]]
        plot!(X[1,:],X[2,:],lc=:blue,legend=false)

    end
    
    # First display the nodes

   display(scatter!(TheMesh.nodes[1,begin:end],TheMesh.nodes[2,begin:end],mc=:red,markersize=1.7,legend=false))

    



# That's all folks ! (O. Bodart, Dec. 2023)
end

function plot2Dmesh!(TheMesh::volcmesh) # plot a mesh over an existing plot
    #
    # Plot all the triangles of a 2D mesh
    #
    Nelem = nbelements(TheMesh)
    if Nelem == 0
        throw(ArgumentError("The provided mesh is empty"))
    end

    # Process to plotting if the mesh if nonempty

    NbElem = nbelements(TheMesh)

    # First element, we add the firstpoint at the end in order to plot closed triangles

    X = zeros(RealNum,2,4)
    TheElem = TheMesh.elements[1:3,1]
    X = [ TheMesh.nodes[:,TheElem[1]]  TheMesh.nodes[:,TheElem[2]]  TheMesh.nodes[:,TheElem[3]] TheMesh.nodes[:,TheElem[1]]]
    plot!(X[1,:],X[2,:],lc=:blue,legend=false)


    # Other elements

    for iloop in 2:NbElem # always make it simple in the first place
        TheElem = TheMesh.elements[1:3,iloop]
        X = [ TheMesh.nodes[:,TheElem[1]]  TheMesh.nodes[:,TheElem[2]]  TheMesh.nodes[:,TheElem[3]] TheMesh.nodes[:,TheElem[1]]]
        plot!(X[1,:],X[2,:],lc=:blue,legend=false)

    end
    
    # First display the nodes

    display(scatter!(TheMesh.nodes[1,begin:end],TheMesh.nodes[2,begin:end],mc=:red,markersize=1.7,legend=false))



# That's all folks ! (O. Bodart, Dec. 2023)
end

function plot2Dmesh(TheMesh::volcmesh, FileName::String)
    #
    # Plot all the triangles of a 2D mesh
    #
    Nelem = nbelements(TheMesh)
    if Nelem == 0
        throw(ArgumentError("The provided mesh is empty"))
    end

    # if we really have a mesh procees to plot

    NbElem = nbelements(TheMesh)

    
    # First element, we add the firstpoint at the end in order to plot closed triangles

    X = zeros(RealNum,2,4)
    TheElem = TheMesh.elements[1:3,1]
    X = [ TheMesh.nodes[:,TheElem[1]]  TheMesh.nodes[:,TheElem[2]]  TheMesh.nodes[:,TheElem[3]] TheMesh.nodes[:,TheElem[1]]]
    plot(X[1,:],X[2,:],lc=:blue,legend=false)


    # Other elements

    for iloop in 2:NbElem # always make it simple in the first place
        TheElem = TheMesh.elements[1:3,iloop]
        X = [ TheMesh.nodes[:,TheElem[1]]  TheMesh.nodes[:,TheElem[2]]  TheMesh.nodes[:,TheElem[3]] TheMesh.nodes[:,TheElem[1]]]
        plot!(X[1,:],X[2,:],lc=:blue,legend=false)

    end
    

    # Last, display the nodes in red colour

    scatter!(TheMesh.nodes[1,begin:end],TheMesh.nodes[2,begin:end],color=:red,markersize=1.7,legend=false)



    savefig(FileName) # To be improved 

# That's all folks ! (O. Bodart, Dec. 2023)
end


function plot2Dline(TheLine::linefrac2D)
    #
    # Plot all the segments of a linfrac2D object over the embedding mesh

    Ne = nbelements(TheLine)
    mid_point = zeros(2, Ne)
    mid_point_x = zeros(Ne)
    mid_point_y = zeros(Ne)
    normal_vec = zeros(2,Ne)
    normal_vec_x = zeros(Ne)
    normal_vec_y = zeros(Ne)

    TheSegment = TheLine.elements[:,1]
    X = [TheLine.nodes[:,TheSegment[1]] TheLine.nodes[:,TheSegment[2]]]
    display(plot(X[1,:],X[2,:],lc=:red,legend=false))

    for iloop = 2:Ne
        #
        TheSegment = TheLine.elements[:,iloop]
        X = [TheLine.nodes[:,TheSegment[1]] TheLine.nodes[:,TheSegment[2]]] 
        plot!(X[1,:],X[2,:],lc=:red,legend=false)
    end

    # Scatter plot for nodes
    scatter!(TheLine.nodes[1, begin:end], TheLine.nodes[2, begin:end], mc=:red, markersize=2, legend=false)

    # Plot the embedding mesh over the figure
    plot2Dmesh!(TheLine.embed_mesh)
end

function plot2Dlinewithnormals(TheLine::linefrac2D)
    #
    # Plot all the segments of a linfrac2D object over the embedding mesh

    Ne = nbelements(TheLine)
    mid_point = zeros(2, Ne)
    mid_point_x = zeros(Ne)
    mid_point_y = zeros(Ne)
    normal_vec = zeros(2,Ne)
    normal_vec_x = zeros(Ne)
    normal_vec_y = zeros(Ne)

    TheSegment = TheLine.elements[:,1]
    X = [TheLine.nodes[:,TheSegment[1]] TheLine.nodes[:,TheSegment[2]]]
    display(plot(X[1,:],X[2,:],lc=:red,legend=false))

    for iloop = 2:Ne
        #
        TheSegment = TheLine.elements[:,iloop]
        X = [TheLine.nodes[:,TheSegment[1]] TheLine.nodes[:,TheSegment[2]]] 
        plot!(X[1,:],X[2,:],lc=:red,legend=false)

        # Calculate the midpoint
        mid_point[:, iloop] = (TheLine.nodes[:, TheSegment[1]] + TheLine.nodes[:, TheSegment[2]]) / 2 

        mid_point_x[iloop] = mid_point[1, iloop]
        mid_point_y[iloop] = mid_point[2, iloop]

         # Get the normal for the current segment
        normal_vec_x[iloop] =  1.0*TheLine.element_normals[1, iloop]
        normal_vec_y[iloop] =  1.0*TheLine.element_normals[2, iloop]
    end

    display(scatter!(mid_point_x, mid_point_y, mc=:blue, markersize=2.5,legend=false))

    display(quiver(mid_point_x, mid_point_y, 
            quiver=(normal_vec_x, normal_vec_y), 
            color=:black, 
            legend=false, 
            arrowsize=0.1,  # Adjust this value for smaller arrow sizes
            linewidth=0.5,  # Adjust line width if needed
            xlabel="Index", 
            ylabel="Displacement", 
            title="Displacement Vector Plot"))


    # Scatter plot for nodes
    scatter!(TheLine.nodes[1, begin:end], TheLine.nodes[2, begin:end], mc=:red, markersize=2, legend=false)


    # Plot the embedding mesh over the figure
    # plot2Dmesh!(TheLine.embed_mesh)
end


#### Moving the points of a mesh according to a 2D dispacement field

function move2Dmesh(TheMesh::volcmesh,U::RealVec,factor::RealNum = 1.)::volcmesh
    # Moves the points of a Mesh deformed by a vector field
    #
    #

    # Create a duplicate of the original mesh

    TempMesh = deepcopy(TheMesh) # we don't want to lose the original mesh

    # Move the pointz of the mesh

    indexvec = 1 # we will scan the meshvector U via this index
    for iloop in 1:nbnodes(TempMesh)
        #
        TempMesh.nodes[:,iloop] += factor*U[indexvec:(indexvec+1)]
        indexvec +=2
    end

    # exit with the result

    return(TempMesh)
end