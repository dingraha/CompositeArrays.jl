using ForwardDiff

function f(x1, x2)
    return @. 2*x1^2 + x2^3
end

function f_d(x1, x2)

    function f_wrap(x)
        return f(x.data[1], x.data[2])
    end

    x_cv = CompositeVector(size.([x1, x2]))
    update!(x_cv, [x1, x2])

    return ForwardDiff.jacobian(f_wrap, x_cv)
end

@Test.test f_d([2.0], [3.0]) == [8.0 27.0]

function f!(y1, y2, x1, x2)
    @. y1 = 2*x1^2 + x2^3
    @. y2 = 0.5*x1^2 - x2^3
    return nothing
end

function f_d!(y1, y2, x1, x2)

    function f_wrap!(y, x)
        return f!(y.data[1], y.data[2], x.data[1], x.data[2])
    end

    x_cv = CompositeVector([x1, x2])
    y_cv = CompositeVector([y1, y2])

    return ForwardDiff.jacobian(f_wrap!, y_cv, x_cv)
end

@Test.test f_d!([0.], [0.], [2.0], [3.0]) == [8.0  27.0;
                                              2.0 -27.0]

function g!(inputs, outputs)
    x1 = inputs["x1"]
    x2 = inputs["x2"]
    y1 = outputs["y1"]
    y2 = outputs["y2"]
    @. y1 = 2*x1^2 + x2^3
    @. y2 = 0.5*x1^2 - x2^3
    return nothing
end

function g_d!(inputs, outputs)

    function g_wrap!(y, x)
        return g!(x.ddata, y.ddata)
    end

    x_ncv = NamedCompositeVector(inputs)
    y_ncv = NamedCompositeVector(outputs)

    return ForwardDiff.jacobian(g_wrap!, y_ncv, x_ncv)
end

@Test.test g_d!(
           Dict("x1"=>[2.0], "x2"=>[3.0]),
           Dict("y1"=>[0.0], "y2"=>[0.0])) == [8.0  27.0; 2.0 -27.0]

