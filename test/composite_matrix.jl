cm = CompositeMatrix(Matrix{Float64}, [(1, 2), (3, 4)], [(2, 3), (4, 5)])
@Test.test size(cm) == (1*2 + 3*4, 2*3 + 4*5)

data = Matrix{Matrix{Float64}}(undef, 2, 2)
data[1, 1] = reshape(collect(1:(1*2*2*3)), 1*2, 2*3)
data[2, 1] = reshape(collect(1:(3*4*2*3)), 3*4, 2*3)
data[1, 2] = reshape(collect(1:(1*2*4*5)), 1*2, 4*5)
data[2, 2] = reshape(collect(1:(3*4*4*5)), 3*4, 4*5)

update!(cm, data)
@Test.test cm[:, 1] == cat(data[1, 1][:, 1], data[2, 1][:, 1], dims=1)
@Test.test cm[:, 2] == cat(data[1, 1][:, 2], data[2, 1][:, 2], dims=1)
@Test.test cm[:, 9] == cat(data[1, 2][:, 3], data[2, 2][:, 3], dims=1)
@Test.test cm[1, :] == cat(data[1, 1][1, :], data[1, 2][1, :], dims=1)
@Test.test cm[2, :] == cat(data[1, 1][2, :], data[1, 2][2, :], dims=1)
@Test.test cm[8, :] == cat(data[2, 1][6, :], data[2, 2][6, :], dims=1)

cm_int = similar(cm, Int64)
@Test.test size(cm_int) == size(cm)

cm_from_data = CompositeMatrix([(1, 2), (3, 4)], [(2, 3), (4, 5)], data)
@Test.test cm_from_data == cm
