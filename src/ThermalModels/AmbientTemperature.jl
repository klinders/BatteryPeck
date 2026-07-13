using DataInterpolations
using DataFrames

function AmbientTemperature(;name, T::DataFrame)
    @parameters t
    @named src = Interpolation(LinearInterpolation, T.temperature, T.time)
    @named clk = ContinuousClock()
    @named output = RealOutput()

    eqs = [
        connect(clk.output, src.input)
        connect(output, src.output)
    ]
    return System(eqs, t; systems=[output, clk, src], name)
end