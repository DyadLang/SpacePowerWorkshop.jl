@mtkmodel MPPT begin
    @parameters begin
        
    end
    @variables begin
        I_in(t)
    end
    @components begin
        p = __JSML__Pin()
        n = __JSML__Pin()
    end
    @equations begin
        v ~ p.v - n.v
        I_in ~ p.i
        p.i + n.i ~ 0
    end
end