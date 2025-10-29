using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical


function SingleCellPack(;name="SingleCellPack", params=Chen2020(), Qcell=5, configuration=(96,3))

    @named p = Pin()
    @named n = Pin()
    @named current = RealInput()
    @named source = Current() # Current source
    @named cell = SPMe(params) # Battery cell model

    @equations begin
        connect(current.u, source.I)
        connect(source.p, cell.p)
        # [connect(battery[i].n, battery[i+1].p) for i in 1:Ncell-1]...
        connect(source.n, cell.n)
    end
end