x = convert(Array{Float64}, reshape(collect(1:5*2), 5, 2))
y = convert(Array{Float64}, reshape(collect(1:5*1), 5, 1))

cv = CompositeVector(size.([x, y]))
update!(cv, [x, y])
@Test.test size(cv) == (5*2 + 5*1,)
@Test.test cv[10] == 10
@Test.test cv[11] == 1
@Test.test cv[1:10] == 1:10
@Test.test cv[10+1:10+5] == 1:5
cv[2] = 8.0
@Test.test cv[2] == 8.0
cv[12] = 55.0
@Test.test cv[12] == 55.0

x_bad = convert(Array{Float64}, reshape(collect(1:5*3), 5, 3))
@Test.test_throws ArgumentError update!(cv, [x_bad, y])

cv_int = similar(cv, Int64)
@Test.test size(cv_int) == size(cv)
@Test.test cv_int.sizes == cv.sizes
@Test.test cv_int.offsets == cv.offsets
for i in eachindex(cv_int)
    cv_int[i] = -i
    @Test.test cv_int[i] == -i
end

cv_from_data = CompositeVector([x, y])
@Test.test cv_from_data == cv
