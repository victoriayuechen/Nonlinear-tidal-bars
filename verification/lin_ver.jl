using Pkg
Pkg.activate(".")

using Gridap
using Gridap.CellData
using WriteVTK
using LineSearches: BackTracking
using GridapGmsh
using Plots

#dir = "swe-solver/output_linear_swe/verification"
#dir_plots = "swe-solver/output_linear_swe/verification/plots"
dir = "output"

function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end

function linear_SWE(order,degree, dt, T, ds, L, B)

    #Parameters
    #B = 10 #y
    #L = 10 #x
    latitude = 52
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    #f=0
    g = 9.81
    H = 0.5 #Constant layer depth at rest
    #T = 50
    t0 = 0

    #dir_pvd = "swe-solver/output_linear_swe/verification/energy"
    dir_pvd = "verification/output"
    #Make model (option 1)

    partitionx::Int64 = round(Int64, L/ds)
    partitiony::Int64 = round(Int64, B/ds)
    domain = (0,L,0,B)
    partition = (partitionx,partitiony)

    model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))  

    ##''''''''''''''Make labels''''''''''''''##
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["bottom","top"]

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)

    # manufactured solutions

    H = 0.5
    u0 = 0.0#1.0

    A3 = (2*π/T)
    A1 = π/(1*L)
    A2 = π/(B)

    B3 = (2*π/T)
    B1 = π/(1*L)
    B2 = π/B

    C3 = (2*π/T)
    C1 = π/(1*L)
    C2 = π/B

    Sh = 0.05
    Su1 = 0.1
    Su2 = 0.1

    # ζ(x, y, t) = Sh*sin(Ay*y)*cos(Ax*(x-L/2))*sin(At*t) + H
    # ux(x, y, t) = u0+Sx*sin(Bt*t)*cos(Bx*(x-L/2))*cos(By*(y-B/2))
    # uy(x, y, t) = Sy*sin(Ct*t)*cos(Cx*(x-L/2))*sin(Cy*y)

    um(x,t::Real) = VectorValue(u0+Su1*sin(B3*t)*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2)), Su2*sin(C3*t)*cos(C1*(x[1]-L/2))*sin(C2*x[2]))
    um(t::Real) = x -> um(x,t)

    ζm(x,t::Real) = Sh*sin(A2*x[2])*cos(A1*(x[1]-L/2))*sin(A3*t) + H
    ζm(t::Real) = x -> ζm(x,t)

    fh(x,t::Real) = H*(C2*Su2*cos(C1*(x[1]-L/2))*cos(C2*x[2])*sin(C3*t)-B1*Su1*sin(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t))+A3*Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*cos(A3*t)
    fh(t::Real) = x -> fh(x,t)

    fu(x,t::Real) = VectorValue(-Su2*f*cos(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)+B3*Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*cos(B3*t)-A1*Sh*g*sin(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t),
    C3*Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*cos(C3*t)+f*(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0)+A2*Sh*g*cos(A1*(x[1]-L/2))*cos(A2*x[2])*sin(A3*t))
    fu(t::Real) = x -> fu(x,t)
    
    #Make reference spaces

    reffe_rt = ReferenceFE(raviart_thomas, Float64, order)
    V = TestFESpace(model,reffe_rt,conformity=:HDiv,dirichlet_tags=DC)#, dirichlet_masks=[(false, true), (false, true)])
    udc(x,t::Real) = VectorValue(0.0, 0.0)
    udc(t::Real) = x -> udc(x,t)
    U = TransientTrialFESpace(V, [um, um])

    reffe_lg = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lg,conformity=:L2)
    hdc(x,t::Real) = 0.0
    hdc(t::Real) = x -> hdc(x,t)
    P = TransientTrialFESpace(Q)

    Y = MultiFieldFESpace([V,Q])
    X = TransientMultiFieldFESpace([U,P])


    #res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) + (∇⋅(w))*g*h + ∂t(h)*ϕ + ϕ*H*(∇⋅(u))-fu(t)⋅w-fh(t)*ϕ)dΩ # - f_u(t)⋅w - f_h(t)*ϕ
    res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) + (∇⋅(w))*g*h + ∂t(h)*ϕ + ϕ*H*(∇⋅(u))-fu(t)⋅w-fh(t)*ϕ)dΩ
    jac(t,(u,h),(du,dh),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dh + ϕ*(H*(∇⋅(du))))dΩ
    jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ

    op = TransientFEOperator(res,jac,jac_t,X,Y)
    #nls = NLSolver(show_trace=true,linesearch=BackTracking())
    nls = LUSolver()

    h_init(x, t::Real) = H + 0.5 * x[1]#0.05*exp(-0.01*(x[1]-0.5)^2 -0.01*(x[2]-2)^2)
    h_init(t::Real) = x -> h_init(x,t)

    u_init(x, t::Real) = VectorValue(0.0, 0.0)
    u_init(t::Real) = x -> u_init(x,t)

    x0 = interpolate_everywhere([um(t0), ζm(t0)], X(t0))
    #x0 = interpolate_everywhere([u_init(t0), h_init(t0)], X(t0))
    
    ode_solver = ThetaMethod(nls,dt,0.5)
    s_h = solve(ode_solver,op,x0,0.0,T)
    

    eul2 = Float64[] # L2 norm error array
    ehl2 = Float64[]
    t_base = Float64[] # time base


    if isdir(dir_pvd)
        output_file = paraview_collection(joinpath(dir_pvd,"new"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
            for (x,t) in s_h
                u_h,h_h = x
                
                eh = (h_h - ζm(t))
                eu = (u_h - um(t))

                push!(ehl2, (sqrt(sum( ∫( eh*eh )*dΩ ))) )
                push!(eul2, (sqrt(sum( ∫( eu⋅eu )*dΩ ))) )
                push!(t_base, t)
                #pvd[t] = createvtk(Ω,joinpath(dir_pvd,"new$t.vtu"),cellfields=["u"=>u_h,"h"=>h_h, "hm"=>ζm(t), "eh"=>eh, "eu"=>eu, "um"=>um(t), "fu"=>fu(t), "fh"=>fh(t)])
                    #"hm"=>h(t), "eh"=>eh, "eu"=>eu, "eh_rel"=>h_grid, "eu_rel"=>u_grid, "fu"=>f_u(t), "fh"=>f_h(t), "um"=>u(t)])
                println("done $t/$T")
            end
        end
    else
        mkdir(dir_pvd)
        output_file = paraview_collection(joinpath(dir_pvd,"new"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
            for (x,t) in s_h
                u_h,h_h = x
                
                eh = (h_h - h(t))
                eu = (u_h - u(t))
                
                #pvd[t] = createvtk(Ω,joinpath(dir_pvd,"new$t.vtu"),cellfields=["u"=>u_h,"h"=>h_h, "hm"=>ζm(t), "eh"=>eh, "eu"=>eu, "um"=>um(t)])
                    #"hm"=>h(t), "eh"=>eh, "eu"=>eu, "eh_rel"=>h_grid, "eu_rel"=>u_grid, "fu"=>f_u(t), "fh"=>f_h(t), "um"=>u(t)])
                println("done $t/$T")
            end
        end
    end

(eul2, ehl2, t_base)
end

#linear_SWE(1, 4, 0.05, 10, 1)

order = 1
degree = 4
Tend = 50
dt = 0.5
ds = 0.01
L = 2
B = 2


#linear_SWE(order, degree, dt, Tend, ds)


eul2 = Float64[]
ehl2 = Float64[]
eul2s = Float64[]
ehl2s = Float64[]
t_base = Float64[]
ts = Float64[]

# Time convergence
time_steps = [1.6, 1.4, 1.2, 1.0 ,0.8, 0.6, 0.4, 0.2, 0.1]
for t in time_steps

    stamp = round(Int64, Tend/t)
    eul2, ehl2, t_base = linear_SWE(order, degree, t, Tend, ds, L, B)
    

    push!(eul2s, sum(eul2)/stamp)
    push!(ehl2s, sum(ehl2)/stamp)
    push!(ts, t)

end

print(ts)
print(eul2s, ehl2s)

# plot(ts,[eul2s, ehl2s],
# label=["L2_u" "L2_h"],
# shape=:auto,
# xlabel="time step",ylabel="averaged error norm",
# yscale = :log, xscale = :log,
# legend = :topleft)
# savefig(joinpath(dir_plots,"time_conv_log.pdf"))



# Spatial convergence
# m_size = [0.64, 0.32, 0.16, 0.08, 0.04, 0.02, 0.01]
# x_base = Float64[]
# dt = 0.1
# for ds in m_size

#     n = Tend/dt
#     # model_name = "swe-solver/meshes/1x4periodic$m.msh"
#     eul2, ehl2, t_base = linear_SWE(order, degree, dt, Tend, ds, L, B)
    
#     push!(eul2s, sum(eul2)/n)
#     push!(ehl2s, sum(ehl2)/n)
#     push!(x_base, L/ds)

# end

# p4 = plot(x_base,[eul2s, ehl2s],
# label=["L2_u" "L2_h"],
# shape=:auto,
# xlabel="mesh size",ylabel="averaged error norm",
# yscale = :log, xscale = :log,
# legend = :topleft)
# savefig(p4, joinpath(dir_plots,"s_conv_log.pdf"))
# print(x_base)
# print(eul2s, ehl2s)