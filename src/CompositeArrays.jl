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
        return new(copy(sizes), calc_offsets(sizes))
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
    if ! all(size.(data) == a.sizes)
        throw(ArgumentError("sizes of data argument not equal to CompositeVector sizes"))
    end
    a.data = data
    return nothing
end

mutable struct NamedCompositeVector{T, N, A} <: AbstractCompositeVector{T, N}
    sizes::Vector{NTuple{N, Int}}
    offsets::Array{Int}
    name2idx::Dict{String, Int}
    data::Vector{A}
    ddata::AbstractDict{String, A}
    function NamedCompositeVector{T, N, A}(sizes::Vector{NTuple{N, Int}}, names) where {T, N, A<:AbstractArray{T, N}}
        name2idx = Dict(name=>i for (i, name) in enumerate(names))
        return new(copy(sizes), calc_offsets(sizes), name2idx)
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
    if ! all(size.(data) == a.sizes)
        throw(ArgumentError("sizes of data argument not equal to NamedCompositeVector sizes"))
    end
    a.data = data
    a.ddata = Dict(name=>a.data[idx] for (name, idx) in a.name2idx)
    return nothing
end

function update!(a::NamedCompositeVector{T, N, A}, ddata::AbstractDict{String, A}) where {T, N, A<:AbstractArray{T, N}}
    # Need to get a list of names in the right order.
    idx2name = Dict(v=>k for (k, v) in a.name2idx)
    names = [idx2name[i] for i in 1:length(idx2name)]

    # Get a Vector of the new data in the correct order.
    data = [ddata[name] for name in names]

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

#function NamedCompositeVector(data::Array{Array{T, N}, 1}, names) where {T, N}
#    n = length(data)
#    if n != length(names)
#        throw(ArgumentError("data and names arguments should be the same length."))
#    end

#    ddata = Dict{String, Array{T, N}}()
#    name2idx = Dict{String, Int}()
#    for i in 1:n
#        ddata[names[i]] = data[i]
#        name2idx[names[i]] = i
#    end

#    offsets = calc_offsets(data)

#    return NamedCompositeVector(data, ddata, offsets, name2idx)
#end

#function NamedCompositeVector(ddata::Dict{String, Array{T, N}}) where {T, N}

#    i = 0
#    name2idx = Dict{String, Int}()
#    data = Array{T, N}[]
#    for (k, v) in ddata
#        i += 1
#        name2idx[k] = i
#        push!(data, v)
#    end

#    offsets = calc_offsets(data)

#    return NamedCompositeVector(data, ddata, offsets, name2idx)
#end

## function update!(ncv::NamedCompositeVector{T, N}, ddata::AbstractDict{String, <:AbstractArray{T, N}}) where {T, N}
#function update!(ncv::NamedCompositeVector{T, N}, ddata::AbstractDict{String, <:AbstractArray}) where {T, N}
#    for (name, a) in ncv.ddata
#        a_new = ddata[name]

#        if ! (a_new === a)
#            if size(a_new) != size(a)
#                ArgumentError("new data for variable $(name) has size $(size(a_new)), but original has size $(size(a))")
#            end

#            ncv.ddata[name] = a_new
#            ncv.data[ncv.name2idx[name]] = a_new
#        end

#    end

#    return nothing
#end

#function Base.similar(ncv::NamedCompositeVector{T, N}, ::Type{S}) where {T,N,S}
#    ddata = Dict{String, Array{S, N}}()
#    for (k, v) in ncv.ddata
#        ddata[k] = similar(v, S)
#    end
#    return NamedCompositeVector(ddata)
#end

#
# Composite Matrix
#

mutable struct CompositeMatrix{T, N, A} <: AbstractCompositeMatrix{T}
    row_sizes::Vector{NTuple{N, Int}}
    col_sizes::Vector{NTuple{N, Int}}
    row_offsets::Vector{Int}
    col_offsets::Vector{Int}
    data::AbstractMatrix{A}
    function CompositeMatrix{T, N, A}(row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}) where {T, N, A<:AbstractMatrix{T}}
        return new(copy(row_sizes), copy(col_sizes), calc_offsets(row_sizes), calc_offsets(col_sizes))
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
            if ! (size(data[i, j]) == expected_size)
                throw(ArgumentError("entry at row $(i), col $(i) in data argument does not match expected size $(expected_size)"))
            end
        end
    end
    a.data = data
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
    data::AbstractMatrix{A}
    ddata::AbstractDict{Tuple{String, String}, A}
    function NamedCompositeMatrix{T, N, A}(row_sizes::Vector{NTuple{N, Int}}, col_sizes::Vector{NTuple{N, Int}}, row_names, col_names) where {T, N, A<:AbstractMatrix{T}}
        if length(row_names) != length(row_sizes)
            raise(ArgumentError("length of row_names and row_sizes arguments differ"))
        end
        if length(col_names) != length(col_sizes)
            raise(ArgumentError("length of col_names and col_sizes arguments differ"))
        end
        row_name2idx = Dict(name=>i for (i, name) in enumerate(row_names))
        col_name2idx = Dict(name=>i for (i, name) in enumerate(col_names))
        return new(copy(row_sizes), copy(col_sizes), copy(row_names), copy(col_names), calc_offsets(row_sizes), calc_offsets(col_sizes), row_name2idx, col_name2idx)
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
            if ! (size(data[i, j]) == expected_size)
                throw(ArgumentError("entry at row $(i), col $(i) in data argument does not match expected size $(expected_size)"))
            end
        end
    end
    a.data = data

    ddata = Dict{Tuple{String, String}, A}()
    for (col_name, col_idx) in a.col_name2idx
        for (row_name, row_idx) in a.row_name2idx
            ddata[(row_name, col_name)] = data[row_idx, col_idx]
        end
    end
    a.ddata = ddata

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

    # ddata = Dict{Tuple{String, String}, Matrix{S}}()
    # for (col_name, col_idx) in b.col_name2idx
    #     for (row_name, row_idx) in b.row_name2idx
    #         ddata[(row_name, col_name)] = b.data[row_idx, col_idx]
    #     end
    # end
    # b.ddata = ddata

    return b
end

#function CompositeMatrix{T}(row_sizes::Array{NTuple{N, Int}, 1}, col_sizes::Array{NTuple{N, Int}, 1}) where {T, N}
#    row_offsets = calc_offsets(row_sizes)
#    col_offsets = calc_offsets(col_sizes)

#    m = length(row_sizes)
#    n = length(col_sizes)
#    data = Array{Array{T, 2}, 2}(undef, m, n)
#    for j in 1:n
#        for i in 1:m
#            data[i, j] = Array{T}(undef, prod(row_sizes[i]), prod(col_sizes[j]))
#        end
#    end

#    return CompositeMatrix{T}(data, row_offsets, col_offsets)
#end

#Base.size(cm::AbstractCompositeMatrix) = (cm.row_offsets[end], cm.col_offsets[end])

#Base.IndexStyle(::Type{<:AbstractCompositeMatrix}) = IndexCartesian()

#function Base.getindex(cm::AbstractCompositeMatrix, I::Vararg{Int, 2})
#    rj, rk = composite_index(cm.row_offsets, I[1])
#    cj, ck = composite_index(cm.col_offsets, I[2])
#    return cm.data[rj, cj][rk, ck]
#end

#function Base.setindex!(cm::AbstractCompositeMatrix, v, I::Vararg{Int, 2})
#    rj, rk = composite_index(cm.row_offsets, I[1])
#    cj, ck = composite_index(cm.col_offsets, I[2])
#    cm.data[rj, cj][rk, ck] = v
#end

#function Base.similar(cm::CompositeMatrix, ::Type{S}) where {S}
#    return CompositeMatrix(similar(cm.data, S), copy(cm.row_offsets), copy(cm.col_offsets))
#end

#mutable struct NamedCompositeMatrix{T} <: AbstractCompositeMatrix{T}
#    data::Array{<:AbstractArray{T, 2}, 2}
#    ddata::Dict{Tuple{String, String}, <:AbstractArray{T, 2}}
#    row_offsets::Array{Int, 1}
#    col_offsets::Array{Int, 1}
#    name2idx::Dict{Tuple{String, String}, Tuple{Int, Int}}
#end

#function NamedCompositeMatrix(ddata::Dict{Tuple{String, String}, Array{T, 2}}) where {T}

#    # I need a sorted list of the row and col keys. And I also need to get the
#    # size of each row and col.
#    rows = Set{String}()
#    cols = Set{String}()
#    row_sizes_all = Dict{String, Array{Int, 1}}()
#    col_sizes_all = Dict{String, Array{Int, 1}}()
#    for ((row_name, col_name), val) in ddata
#        push!(rows, row_name)
#        push!(cols, col_name)

#        m, n = size(val)

#        if ! haskey(row_sizes_all, row_name)
#            row_sizes_all[row_name] = Array{Int, 1}()
#        end
#        push!(row_sizes_all[row_name], m)

#        if ! haskey(col_sizes_all, col_name)
#            col_sizes_all[col_name] = Array{Int, 1}()
#        end
#        push!(col_sizes_all[col_name], n)
#    end
#    rows = sort(collect(rows))
#    cols = sort(collect(cols))

#    # Check that the row and col sizes are all the same
#    row_sizes = Array{Int, 1}()
#    for k in rows
#        v = row_sizes_all[k]
#        if all(v .== v[1])
#            push!(row_sizes, v[1])
#        else
#            println("bad row sizes for $(k): $(v)")
#        end
#    end

#    col_sizes = Array{Int, 1}()
#    for k in cols
#        v = col_sizes_all[k]
#        if all(v .== v[1])
#            push!(col_sizes, v[1])
#        else
#            println("bad col sizes for $(k): $(v)")
#        end
#    end

#    # Now initialize the data.
#    m = length(rows)
#    n = length(cols)
#    data = Array{Array{T, 2}, 2}(undef, m, n)
#    name2idx = Dict{Tuple{String, String}, Tuple{Int, Int}}()
#    for j in 1:n
#        for i in 1:m
#            row, col = rows[i], cols[j]
#            # Assuming that every possible row/col combination exists in ddata.
#            data[i, j] = ddata[(row, col)]
#            name2idx[(row, col)] = (i, j)
#        end
#    end

#    row_offsets = calc_offsets(row_sizes)
#    col_offsets = calc_offsets(col_sizes)

#    return NamedCompositeMatrix(data, ddata, row_offsets, col_offsets, name2idx)
#end

#function Base.similar(ncm::NamedCompositeMatrix{T}, ::Type{S}) where {T, S}
#    ddata = Dict{Tuple{String, String}, Array{S, 2}}()
#    for (k, v) in ncm.ddata
#        ddata[k] = similar(v, S)
#    end
#    return NamedCompositeMatrix(ddata)
#end

## function update!(ncm::NamedCompositeMatrix{T}, ddata::AbstractDict{Tuple{String, String}, <:AbstractArray{T, 2}}) where {T}
#function update!(ncm::NamedCompositeMatrix{T}, ddata::AbstractDict{Tuple{String, String}, <:AbstractArray}) where {T}
#    for (name, a) in ncm.ddata
#        a_new = ddata[name]

#        if ! (a_new === a)
#            println("updating $(name)")
#            if size(a_new) != size(a)
#                ArgumentError("new data for variable $(name) has size $(size(a_new)), but original has size $(size(a))")
#            end

#            ncm.ddata[name] = a_new
#            ncm.data[ncm.name2idx[name]...] = a_new
#        else
#            println("not updating $(name)")
#        end

#    end

#    for (name, a) in ncm.ddata
#        a_new = ddata[name]
#        println("a_new === a: $(a_new === a)")
#        println("typeof(a_new): $(typeof(a_new))")
#        println("typeof(a): $(typeof(a))")
#    end

#    return nothing
#end

end # module
