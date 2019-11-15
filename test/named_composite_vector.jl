x = convert(Array{Float64}, reshape(collect(1:5*2), 5, 2))
y = convert(Array{Float64}, reshape(collect(1:5*1), 5, 1))

ncv = NamedCompositeVector(Array{Float64, 2}, size.([x, y]), ["x", "y"])
update!(ncv, [x, y])
@Test.test ncv.ddata["x"] == x
@Test.test ncv.ddata["y"] == y
@Test.test size(ncv) == (5*2 + 5*1,)
@Test.test ncv[10] == 10
@Test.test ncv[11] == 1
@Test.test ncv[1:10] == 1:10
@Test.test ncv[10+1:10+5] == 1:5
ncv[2] = 8.0
@Test.test ncv[2] == 8.0
ncv[12] = 55.0
@Test.test ncv[12] == 55.0

x_bad = convert(Array{Float64}, reshape(collect(1:5*3), 5, 3))
@Test.test_throws ArgumentError update!(ncv, [x_bad, y])

ncv_int = similar(ncv, Int64)
@Test.test size(ncv_int) == size(ncv)
@Test.test ncv_int.sizes == ncv.sizes
@Test.test ncv_int.offsets == ncv.offsets
@Test.test ncv_int.name2idx == ncv.name2idx
for i in eachindex(ncv_int)
    ncv_int[i] = -i
    @Test.test ncv_int[i] == -i
end

ncv_from_data = NamedCompositeVector([x, y], ["x", "y"])
@Test.test ncv_from_data == ncv

ncv_from_ddata = NamedCompositeVector(Dict("x"=>x, "y"=>y))
@Test.test ncv_from_ddata == ncv

update!(ncv_int, Dict("x"=>convert(Array{Int64}, x.+1), "y"=>convert(Array{Int64}, y.+1)))
@Test.test ncv_int == ncv .+ 1
