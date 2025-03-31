include("commons.jl")


function write_vector_to_tex(vec::Vector, filename::String)
    open(filename, "w") do file
        # Write the vector in LaTeX format
        write(file, "\\begin{bmatrix}\n")
        for v in vec
            write(file, "$v \n")
        end
        write(file, "\\end{bmatrix}\n")
    end
end

function read_vector_from_tex(filename::String)
    vec = Float64[]  # Start with an empty vector to store numbers
    open(filename, "r") do file
        # Read the file line by line
        for line in eachline(file)
            # Try to parse the number and store it in the vector
            try
                push!(vec, parse(Float64, strip(line)))  # Parse each line as a Float64
            catch e
                continue  # Skip lines that don't contain valid numbers
            end
        end
    end
    return vec
end