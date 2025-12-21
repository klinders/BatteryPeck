using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using Unitful

function WangEtAl2014(; name)


    @named temperature = RealInput()
    @named current = RealInput()

    @independent_variables t

    D = Differential(t)

    @constants begin
        Q = 3.6 # Nominal capacity in Ah

        # Wang et al. (2014) model parameters
        a = 8.61E-6
        b = -5.13E-3
        c = 7.63E-1
        d = -6.7E-3
        e = 2.35
        f = 14876
        Ea = 24500
        R = 8.314
    end

    @variables begin
        ah_throughput(t) = 0.0 # Ah throughput
        calendar_loss(t) = 0.0 # Calendar loss in Ah
        cyclic_loss(t) = 0.0 # Cyclic loss in Ah
        capacity_loss(t) # Total capacity loss in Ah
        c_rate(t)
        T(t)
        I(t)
    end
    
    eqs = [
        T ~ temperature.u
        I ~ current.u

        c_rate ~ abs(I)/Q # C-rate calculation
        D(ah_throughput) ~ abs(I)

        # Cyclic loss calculation
        D(cyclic_loss) ~ (a*(T^2) + b*T + c)*exp((d*T + e)*c_rate) * ah_throughput

        # Calendar loss calculation
        D(calendar_loss) ~  f*exp(-Ea/(R*T))*(1/(2))*((t*24*3600+1e-6)^-0.5)

        # Total capacity loss
        capacity_loss ~ cyclic_loss + calendar_loss
    ]

    model = System(eqs, t, name=name)

    compose(model, temperature, current)

end