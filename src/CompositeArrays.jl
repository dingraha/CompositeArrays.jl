module CompositeArrays

export CompositeVector, NamedCompositeVector, CompositeMatrix, NamedCompositeMatrix, update!

#
# Abstract types. I guess these need to be defined before they're used.
#
abstract type AbstractCompositeVector{T, N} <: AbstractArray{T, 1} end
abstract type AbstractCompositeMatrix{T} <: AbstractArray{T, 2} end

#
# Utilities.
#

function composite_index(offsets::Array{Int, 1}, i::Int)
    j = 0
    for offset in offsets
        if i <= offset
            break
        end
        j += 1
    end
    k = i - offsets[j]
    return j, k
end

function composite_index(cv::AbstractCompositeVector, i::Int)
    return composite_index(cv.offsets, i)
end

function calc_offsets(sizes::Array{NTuple{N, Int}, 1}) where {N}
    offsets = [0]
    n = length(sizes)
    for i in 1:n
        o = offsets[i] + prod(sizes[i])
        push!(offsets, o)
    end
    return offsets
end

function calc_offsets(data::Array{Array{T, N}, 1}) where {T,N}
    return calc_offsets(size.(data))
end

function calc_offsets(sizes::Array{Int, 1})
    return calc_offsets(Tuple.(sizes))
end

#
# Composite Vectors
#

mutable struct CompositeVector{T, N, A} <: AbstractCompositeVector{T, N}
    sizes::Vector{NTuple{N, Int}}
    offsets::Vector{Int}
    data::Vector{A}
    function CompositeVector{T, N, A}(sizes::Vector{NTuple{N, Int}}) where {T, N, A<:AbstractArray{T, N}}
        data = Vector{A}(undef, length(sizes))
        return new(copy(sizes), calc_offsets(sizes), data)
    end
end

function CompositeVector(::Type{A}, sizes::Vector{NTuple{N, Int}}) where {T, N, A<:AbstractArray{T, N}}
    return CompositeVector{T, N, A}(sizes)
end

# Default `T` to `Float64` and `A` to `Array{Float64, N}`.
function CompositeVector(sizes::Vector{NTuple{N, Int}}) where {N}
    return CompositeVector{Float64, N, Array{Float64, N}}(sizes)
end

function CompositeVector(data::Vector{A}) where {T, N, A<:AbstractArray{T, N}}
    cv = CompositeVector{T, N, A}(size.(data))
    update!(cv, data)
    return cv
end

Base.IndexStyle(::Type{<:AbstractCompositeVector}) = IndexLinear()

Base.size(a::AbstractCompositeVector) = (a.offsets[end],)

function Base.getindex(a::AbstractCompositeVector, i::Int)
    j, k = composite_index(a, i)
    return a.data[j][k]
end

function Base.setindex!(a::AbstractCompositeVector, v, i::Int)
    j, k = composite_index(a, i)
    a.data[j][k] = v
end

function Base.similar(a::CompositeVector{T, N, A}, ::Type{S}) where {T, N, A, S}
    b = CompositeVector{S, N, Array{S, N}}(a.sizes)

    data = Vector{Array{S, N}}()
    for sz in b.sizes
        push!(data, Array{S}(undef, sz...))
    end
    update!(b, data)

    return b
end

function update!(a::CompositeVector, data)
    for (i, expected_size) in enumerate(a.sizes)
        if size(data[i]) == expected_size
            a.data[i] = data[i]
        else
            throw(ArgumentError("sizes of data argument not equal to CompositeVector sizes"))
        end
    end
    return nothing
end

mutable struct NamedCompositeVector{T, N, A} <: AbstractCompositeVector{T, N}
    sizes::Vector{NTuple{N, Int}}
    names::Vector{String}
    offsets::Array{Int}
    name2idx::Dict{String, Int}
    data::Vector{A}
    ddata::Dict{String, A}
    function NamedCompositeVector{T, N, A}(sizes::Vector{NTuple{N, Int}}, names) where {T, N, A<:AbstractArray{T, N}}
        if length(names) != length(sizes)
            raise(ArgumentError("length of names and sizes arguments differ"))
        end
        name2idx = Dict(name=>i for (i, name) in enumerate(names))
        data = Vector{A}(undef, length(sizes))
        ddata = Dict{String, A}()
        return new(copy(sizes), copy(names), calc_offsets(sizes), name2idx, data, ddata)
    end
end

function NamedCompositeVector(::Type{A}, sizes::Vector{NTuple{N, Int}}, names) where {T, N, A<:AbstractArray{T, N}}
    return NamedCompositeVector{T, N, A}(sizes, names)
end

function NamedCompositeVector(data::AbstractVector{A}, names) where {T, N, A<:AbstractArray{T, N}}
    a = NamedCompositeVector{T, N, A}(size.(data), names)
    update!(a, data)
    return a
end


function NamedCompositeVector(ddata::AbstractDict{String, A}) where {T, N, A<:AbstractArray{T, N}}
    names = sort([k for k in keys(ddata)])
    data = [ddata[name] for name in names]
    a = NamedCompositeVector{T, N, A}(size.(data), names)
    update!(a, data)
    return a
end


function update!(a::NamedCompositeVector{T, N, A}, data::AbstractVector{A}) where {T, N, A<:AbstractArray{T, N}}
    for (i, expected_size) in enumerate(a.sizes)
        if size(data[i]) == expected_size
            a.data[i] = data[i]
            a.ddata[a.names[i]] = data[i]
        else
            throw(ArgumentError("sizes of data argument not equal to NamedCompositeVector sizes"))
        end
    end
    return nothing
end

function update!(a::NamedCompositeVector{T, N, A}, ddata::AbstractDict{String, A}) where {T, N, A<:AbstractArray{T, N}}
    # Get a Vector of the new data in the correct order.
    data = [ddata[name] for name in a.names]

    # Do it.
    update!(a, data)

    return nothing
end

function Base.similar(a::NamedCompositeVector{T, N, A}, ::Type{S}) where {T, N, A, S}

    # Need to get a list of names in the right order.
    idx2name = Dict(v=>k for (k, v) in a.name2idx)
    names = [idx2name[i] for i in 1:length(idx2name)]

    b = NamedCompositeVector{S, N, Array{S, N}}(a.sizes, names)

    data = Vector{Array{S, N}}()
    for sz in b.sizes
        push!(data, Array{S}(undef, sz...))
    end
    update!(b, data)

    return b
end

#
# Composite Matrix
#

mutable struct CompositeMatrix{T, N, A} <: AbstractCompositeMatrix{T}
    row_sizes::Vector{NTuple{N, Int}}
    col_sizes::Vector{NTuple{N, Int}}
    row_offsets::Vector{Int}
    col_offsets::Vector{Int}
    data::Matrix{A}
    function CompositeMatrix{T, N, A}(row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}) where {T, N, A<:AbstractMatrix{T}}
        data = Matrix{A}(undef, length(row_sizes), length(col_sizes))
        return new(copy(row_sizes), copy(col_sizes), calc_offsets(row_sizes), calc_offsets(col_sizes), data)
    end
end

function CompositeMatrix(::Type{A}, row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}) where {T, N, A<:AbstractMatrix{T}}
    return CompositeMatrix{T, N, A}(row_sizes, col_sizes)
end

function CompositeMatrix(row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}, data::AbstractMatrix{A}) where {T, N, A<:AbstractMatrix{T}}
    a = CompositeMatrix{T, N, A}(row_sizes, col_sizes)
    update!(a, data)
    return a
end

Base.IndexStyle(::Type{<:AbstractCompositeMatrix}) = IndexCartesian()

Base.size(a::AbstractCompositeMatrix) = (a.row_offsets[end], a.col_offsets[end])

function Base.getindex(a::AbstractCompositeMatrix, I::Vararg{Int, 2})
    rj, rk = composite_index(a.row_offsets, I[1])
    cj, ck = composite_index(a.col_offsets, I[2])
    return a.data[rj, cj][rk, ck]
end

function Base.setindex!(a::AbstractCompositeMatrix, v, I::Vararg{Int, 2})
    rj, rk = composite_index(a.row_offsets, I[1])
    cj, ck = composite_index(a.col_offsets, I[2])
    a.data[rj, cj][rk, ck] = v
end

function Base.similar(a::CompositeMatrix{T, N, A}, ::Type{S}) where {T, N, A, S}
    b = CompositeMatrix{S, N, Matrix{S}}(a.row_sizes, a.col_sizes)
    data = Matrix{Matrix{S}}(undef, length(b.row_sizes), length(b.col_sizes))
    for j in 1:length(b.col_sizes)
        col_size = prod(b.col_sizes[j])
        for i in 1:length(b.row_sizes)
            row_size = prod(b.row_sizes[i])
            data[i, j] = Matrix{S}(undef, row_size, col_size)
        end
    end
    return b
end

function update!(a::CompositeMatrix{T, N, A}, data::Matrix{A}) where {T, N, A<:AbstractMatrix{T}}
    # Need to iterate over data, checking that each entry is consistent with
    # the coresponding entry in row_sizes and col_sizes.
    for j in 1:length(a.col_sizes)
        col_size = prod(a.col_sizes[j])
        for i in 1:length(a.row_sizes)
            row_size = prod(a.row_sizes[i])
            expected_size = (row_size, col_size)
            if size(data[i, j]) == expected_size
                a.data[i, j] = data[i, j]
            else
                throw(ArgumentError("entry at row $(i), col $(i) in data argument does not match expected size $(expected_size)"))
            end
        end
    end

    return nothing
end

mutable struct NamedCompositeMatrix{T, N, A} <: AbstractCompositeMatrix{T}
    row_sizes::Vector{NTuple{N, Int}}
    col_sizes::Vector{NTuple{N, Int}}
    row_names::Vector{String}
    col_names::Vector{String}
    row_offsets::Vector{Int}
    col_offsets::Vector{Int}
    row_name2idx::Dict{String, Int}
    col_name2idx::Dict{String, Int}
    data::Matrix{A}
    ddata::Dict{Tuple{String, String}, A}
    function NamedCompositeMatrix{T, N, A}(row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}, row_names, col_names) where {T, N, A<:AbstractMatrix{T}}
        if length(row_names) != length(row_sizes)
            raise(ArgumentError("length of row_names and row_sizes arguments differ"))
        end
        if length(col_names) != length(col_sizes)
            raise(ArgumentError("length of col_names and col_sizes arguments differ"))
        end
        row_name2idx = Dict(name=>i for (i, name) in enumerate(row_names))
        col_name2idx = Dict(name=>i for (i, name) in enumerate(col_names))

        data = Matrix{A}(undef, length(row_sizes), length(col_sizes))
        ddata = Dict{Tuple{String, String}, A}()
                
        return new(copy(row_sizes), copy(col_sizes), copy(row_names), copy(col_names), calc_offsets(row_sizes), calc_offsets(col_sizes), row_name2idx, col_name2idx, data, ddata)
    end
end

function NamedCompositeMatrix(::Type{A}, row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}, row_names, col_names) where {T, N, A<:AbstractMatrix{T}}
    return NamedCompositeMatrix{T, N, A}(row_sizes, col_sizes, row_names, col_names)
end

function NamedCompositeMatrix(row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}, row_names, col_names, data::AbstractMatrix{A}) where {T, N, A<:AbstractMatrix{T}}
    a = NamedCompositeMatrix{T, N, A}(row_sizes, col_sizes, row_names, col_names)
    update!(a, data)
    return a
end

function NamedCompositeMatrix(row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}, row_names, col_names, ddata::AbstractDict{Tuple{String, String}, A}) where {T, N, A<:AbstractMatrix{T}}
    a = NamedCompositeMatrix{T, N, A}(row_sizes, col_sizes, row_names, col_names)
    update!(a, ddata)
    return a
end

# Most of this is the same as CompositeMatrix. I think I could make
# CompositeMatrix's version work with AbstractCompositeMatrix, then just call
# that from here, and update the ddata attribute after that.
function update!(a::NamedCompositeMatrix{T, N, A}, data::Matrix{A}) where {T, N, A<:AbstractMatrix{T}}
    # Need to iterate over data, checking that each entry is consistent with
    # the coresponding entry in row_sizes and col_sizes.
    for j in 1:length(a.col_sizes)
        col_size = prod(a.col_sizes[j])
        for i in 1:length(a.row_sizes)
            row_size = prod(a.row_sizes[i])
            expected_size = (row_size, col_size)
            if size(data[i, j]) == expected_size
                a.data[i, j] = data[i, j]
                a.ddata[a.row_names[i], a.col_names[j]] = data[i, j]
            else
                throw(ArgumentError("entry at row $(i), col $(i) in data argument does not match expected size $(expected_size)"))
            end
        end
    end

    return nothing
end

function update!(a::NamedCompositeMatrix{T, N, A}, ddata::Dict{Tuple{String, String}, A}) where {T, N, A<:AbstractMatrix{T}}
    # Need to iterate over ddata, checking that each entry is consistent with
    # the coresponding entry in row_sizes and col_sizes.
    for ((row_name, col_name), val) in ddata
        i = a.row_name2idx[row_name]
        row_size = prod(a.row_sizes[i])

        j = a.col_name2idx[col_name]
        col_size = prod(a.col_sizes[j])

        expected_size = (row_size, col_size)

        if size(val) == expected_size
            a.data[i, j] = val
            a.ddata[row_name, col_name] = val
        else
            throw(ArgumentError("entry at row $(row_name), col $(col_name) in data argument does not match expected size $(expected_size)"))
        end
    end

    return nothing
end

function Base.similar(a::NamedCompositeMatrix{T, N, A}, ::Type{S}) where {T, N, A, S}
    b = NamedCompositeMatrix{S, N, Matrix{S}}(a.row_sizes, a.col_sizes, a.row_names, a.col_names)
    data = Matrix{Matrix{S}}(undef, length(b.row_sizes), length(b.col_sizes))
    for j in 1:length(b.col_sizes)
        col_size = prod(b.col_sizes[j])
        for i in 1:length(b.row_sizes)
            row_size = prod(b.row_sizes[i])
            data[i, j] = Matrix{S}(undef, row_size, col_size)
        end
    end

    update!(b, data)

    return b
end

end # module
