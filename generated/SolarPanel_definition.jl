### DO NOT EDIT THIS FILE
### This file is auto-generated by the Dyad command-line compiler.
### If you edit this code it is likely to get overwritten.
### Instead, update the Dyad source code and regenerate this file


@doc Markdown.doc"""
   SolarPanel(; name)
"""
@component function SolarPanel(; name)

  ### Symbolic Parameters
  __params = Any[]

  ### Variables
  __vars = Any[]

  ### Constants
  __constants = Any[]

  ### Components
  __systems = ODESystem[]
  push!(__systems, @named cell = SpacePowerWorkshop.PVCell())
  push!(__systems, @named converter = SpacePowerWorkshop.DCDC_MPPT())
  push!(__systems, @named i = ElectricalComponents.CurrentSensor())
  push!(__systems, @named ground = ElectricalComponents.Ground())
  push!(__systems, @named temp = SpacePowerWorkshop.TempSensor())

  ### Defaults
  __defaults = Dict()

  ### Initialization Equations
  __initialization_eqs = []

  ### Equations
  __eqs = Equation[]
  push!(__eqs, cell.G ~ temp.G_eff)
  push!(__eqs, cell.T_reading ~ temp.T)
  push!(__eqs, connect(ground.g, cell.n, converter.n))
  push!(__eqs, connect(cell.p, i.n))
  push!(__eqs, connect(i.p, converter.p))
  push!(__eqs, cell.Vt ~ converter.Vt)

  # Return completely constructed ODESystem
  return ODESystem(__eqs, t, __vars, __params; systems=__systems, defaults=__defaults, name, initialization_eqs=__initialization_eqs)
end
export SolarPanel

Base.show(io::IO, a::MIME"image/svg+xml", t::typeof(SolarPanel)) = print(io,
  """<div style="height: 100%; width: 100%; background-color: white"><div style="margin: auto; height: 500px; width: 500px; padding: 200px"><svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 1000 1000"
    overflow="visible" shape-rendering="geometricPrecision" text-rendering="geometricPrecision">
      <defs>
        <filter id='red-shadow' color-interpolation-filters="sRGB"><feDropShadow dx="0" dy="0" stdDeviation="100" flood-color="#ff0000" flood-opacity="0.5"/></filter>
        <filter id='green-shadow' color-interpolation-filters="sRGB"><feDropShadow dx="0" dy="0" stdDeviation="100" flood-color="#00ff00" flood-opacity="0.5"/></filter>
        <filter id='blue-shadow' color-interpolation-filters="sRGB"><feDropShadow dx="0" dy="0" stdDeviation="100" flood-color="#0000ff" flood-opacity="0.5"/></filter>
        <filter id='drop-shadow' color-interpolation-filters="sRGB"><feDropShadow dx="0" dy="0" stdDeviation="40" flood-opacity="0.5"/></filter>
      </defs>
    
      </svg></div></div>""")
