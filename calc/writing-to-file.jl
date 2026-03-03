# The general shape for which this is tested
shape = (20, 10, 3)


a = collect(range(1, 600; step=1))

b = rand(range(0, 1; step=0.01), shape)

a = reshape(a, shape)

c = a .+ b

"""
Write an array to a file
"""
function write_to_file(filename::AbstractString, tensor::Array{T,3}) where T<:Number
    file = open(filename, "w")

    for islow in axes(tensor, 3)
        if islow == 1
        else
            for ifast in axes(tensor, 1)
                for imedium in axes(tensor, 2)
                    write(file, string(tensor[ifast, imedium, islow]) * " ")
                end
                write(file, "\n")
            end

            for ifast in axes(tensor, 1)
                for imedium in axes(tensor, 2)
                    write(file, string(tensor[ifast, imedium, islow]) * " ")
                end
                write(file, "\n")
            end
        end

    end
    write(file, "\n")

    close(file)
end

write_to_file("test.data", c)

"""
Read from a file into an array
"""
function read_from_file!(filename::AbstractString, tensor::Array{T,3}) where T<:Number
    file = open(filename, "r")

    for islow in axes(tensor, 3)
        if islow == 1
            readline(file)
        else
            readline(file)
            readline(file)
        end
        for ifast in axes(tensor, 1)
            tensor[ifast, :, islow] = parse.(T, split(readline(file)))
        end
    end
    close(file)
end

# g = Array{Float64}(undef, shape)
# read_from_file!("test.data", g)
#
# # Is written and read data the same?
# println(isequal(c, g))

