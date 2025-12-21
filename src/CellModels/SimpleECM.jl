using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

@mtkmodel CellVoltage begin
    @extend OnePort()
    @parameters begin
        Q = 3.6 # Capacity in Ampere-hours
        Voc0 = 3.7 # Open-circuit voltage in Volts
        a = 0.5 # Coefficient for voltage equation
    end
    @variables begin
        z(t) = 1.0 # State of charge (SOC) as a variable
    end
    # Equations
    @equations begin
        D(z) ~ -i / (Q*3600)
        v ~ Voc0 + a * z # Example equation for voltage source
    end
end

@mtkmodel SimpleECM begin
    @extend OnePort()
    @parameters begin
        Q = 3.6 # capacity in Ampere-hours
        R0 = 0.01 # Resistance in Ohms
        R1 = 0.01 # Resistance in Ohms
        C1 = 500.0 # Capacitance in Farads
    end
    @components begin
        # Resistors and Capacitors
        r0 = Resistor(R=R0)
        r1 = Resistor(R=R1)
        c1 = Capacitor(C=C1,v=0.0)
        # Voltage source
        voc = CellVoltage(Q=Q)
    end

    # Equations
    @equations begin
        connect(r0.n, p)
        connect(r0.p, r1.n)
        connect(r0.p, c1.n)
        connect(r1.p, voc.p)
        connect(c1.p, voc.p)
        connect(voc.n, n)
    end
end


# Example usage:
# using DifferentialEquations
# sys = SingleECM(R=0.01, C=1000.0, V_oc=3.7)
# u0 = [V => 3.7, I => 0.0, V_R => 0.0, V_C => 0.0, Q => 3600.0]
# prob = ODEProblem(sys, u0, (0.0, 100.0))
# sol = solve(prob)
