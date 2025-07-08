using Roots            # Julia example

q = 1.602e-19
Ns = 1
k = 1.38e-23
qok = 11605.0
T0 = 301.15
Isc = 27*18e-3
Im = 27*17.5e-3
Vm = 2.406
Voc = 2.72
β = -5.6e-3
α = 10e-6

λ   = q/(Ns*k*T0)
y   = Im/Isc
ν   = Vm/Voc           # exponent
f(s) = (Im / Isc) - exp(λ*Vm/s) + (Isc - Im)/Isc*exp(λ*Voc/s)

f1(s) = exp(λ*Vm/s)
f2(s) = (Isc - Im)/Isc*exp(λ*Voc/s)

s0  = 3.37             # sensible initial guess
A   = find_zero(f, s0) # positive root of Eq. (2)

f(A)

# eq 8
Irs = Isc/(exp(λ*Voc/A) - 1)
Irs = (Isc - Im)/(exp(λ * Vm / A) - 1)

Is(G, ΔT) = (exp(abs(β)*ΔT*λ/A)*G*(Isc + α*ΔT)) / ((G*Isc/Irs + 1)^(T0/(ΔT + T0)) - exp((abs(β)*ΔT*λ)/A))

sat(v) = exp(λ*v/A) - 1