# Julvolc type definitions
#

# O. Bodart 2021/09/22

################
# Basic common types aliases
################

RealNum = Float64
RealMat = Matrix{RealNum}
RealVec = Vector{RealNum}

macro EmptyRealMat()
    return(Array{RealNum}(undef,0, 0))
end

macro EmptyRealVec()
    return(Array{RealNum}(undef,0))
end


Index = Int64
IndexMat = Matrix{Index}
IndexVec = Vector{Index}


macro EmptyIntMat()
    return(Array{Index}(undef,0, 0))
end

macro EmptyIntVec()
    return(Array{Index}(undef,0))
end
