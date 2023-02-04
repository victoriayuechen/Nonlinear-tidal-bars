using Plots

L = 100
B = 20
T = 100

H = 0.5
u0 = 1.0

At = (2*π/T)
Ax = π/(1.5*L)
Ay = π/(B)

Bt = (2*π/T)
Bx = π/(1.5*L)
By = π/B

Ct = (2*π/T)
Cx = π/(1.5*L)
Cy = π/B

Sh = 0.05
Sx = 0.1
Sy = 0.1

#plotting

step = 0.1
partx::Int64 = L / step
party::Int64 = B/ step

x = range(0, L, partx)
y = range(0, B, party)

ti = 90

ζ(x, y, t) = Sh*sin(Ay*y)*cos(Ax*(x-L/2))*sin(At*t) + H
ux(x, y, t) = u0+Sx*sin(Bt*t)*cos(Bx*(x-L/2))*cos(By*(y-B/2))
uy(x, y, t) = Sy*sin(Ct*t)*cos(Cx*(x-L/2))*sin(Cy*y)

Fh = Float64[ζ(ix,iy, ti) for ix in x, iy in y]' #convert f(x,y) to an array
Fux = Float64[ux(ix,iy, ti) for ix in x, iy in y]'
Fuy = Float64[uy(ix,iy, ti) for ix in x, iy in y]'

plot(x, y, Fuy,st=:heatmap)
plot!(xlabel="x",ylabel="y")