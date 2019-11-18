ncm = NamedCompositeMatrix(Matrix{Float64}, [(1, 2), (3, 4)], [(2, 3), (4, 5)], ["y1", "y2"], ["x1", "x2"])
@Test.test size(ncm) == (1*2 + 3*4, 2*3 + 4*5)

data = Matrix{Matrix{Float64}}(undef, 2, 2)
data[1, 1] = reshape(collect(1:(1*2*2*3)), 1*2, 2*3)
data[2, 1] = reshape(collect(1:(3*4*2*3)), 3*4, 2*3)
data[1, 2] = reshape(collect(1:(1*2*4*5)), 1*2, 4*5)
data[2, 2] = reshape(collect(1:(3*4*4*5)), 3*4, 4*5)
update!(ncm, data)

@Test.test ncm[:, 1] == cat(data[1, 1][:, 1], data[2, 1][:, 1], dims=1)
@Test.test ncm[:, 2] == cat(data[1, 1][:, 2], data[2, 1][:, 2], dims=1)
@Test.test ncm[:, 9] == cat(data[1, 2][:, 3], data[2, 2][:, 3], dims=1)
@Test.test ncm[1, :] == cat(data[1, 1][1, :], data[1, 2][1, :], dims=1)
@Test.test ncm[2, :] == cat(data[1, 1][2, :], data[1, 2][2, :], dims=1)
@Test.test ncm[8, :] == cat(data[2, 1][6, :], data[2, 2][6, :], dims=1)

ncm_int = similar(ncm, Int64)
@Test.test size(ncm_int) == size(ncm)

ncm_from_data = NamedCompositeMatrix([(1, 2), (3, 4)], [(2, 3), (4, 5)], ["y1", "y2"], ["x1", "x2"], data)
@Test.test ncm_from_data == ncm

@Test.test ncm.ddata[("y1", "x1")] == data[1, 1]
@Test.test ncm.ddata[("y1", "x2")] == data[1, 2]
@Test.test ncm.ddata[("y2", "x1")] == data[2, 1]
@Test.test ncm.ddata[("y2", "x2")] == data[2, 2]

ddata = Dict(("y1", "x1")=>data[1, 1],
             ("y1", "x2")=>data[1, 2],
             ("y2", "x1")=>data[2, 1],
             ("y2", "x2")=>data[2, 2])

ncm_from_ddata = NamedCompositeMatrix([(1, 2), (3, 4)], [(2, 3), (4, 5)], ["y1", "y2"], ["x1", "x2"], ddata)
@test ncm_from_ddata == ncm
