using Pkg
Pkg.activate(".")

using Gridap
using Gridap.CellData
using WriteVTK
using LineSearches: BackTracking
using GridapGmsh
using Plots
# include("mesh_generator.jl")
# using .MyMeshGenerator


dir = "swe-solver/output_linear_swe/verification"
dir_plots = "swe-solver/output_linear_swe/verification/plots"

function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end


function linear_SWE(order,degree, dt, tstamp, model_nr)

    #Parameters
    W = 4.0 #Channel width
    L = 1.0 #Channel Length
    latitude = 52
    η = 7.29e-5
    #f = 2*η*sin(latitude*(π/180))
    f=0
    g = 9.81
    H = 0.5 #Constant layer depth at rest
    Tend = 50
    t0 = 0

    dir_pvd = "swe-solver/output_linear_swe/verification/energy"

    #Make model (option 1)
    # domain = (0,W,0,L)
    # partition = (20,20)

    #Make model (option 2)
    # modelfile = "swe-solver/meshes/1x2periodic_fine.msh"
    # generate_rectangle_mesh(Float32(W), Float32(L), modelfile, "rectangle", Float32(0.01), true)
    # model = GmshDiscreteModel(modelfile)

    #Load up the model
    model = GmshDiscreteModel("swe-solver/meshes/1x4periodic$model_nr.msh")

    #Make labels
    # labels = get_face_labeling(model)
    # add_tag_from_tags!(labels,"bottom",[1,2,5])
    # add_tag_from_tags!(labels,"left",[7])
    # add_tag_from_tags!(labels,"right",[8])
    # add_tag_from_tags!(labels,"top",[3,4,6])
    # add_tag_from_tags!(labels,"inside",[9])
    # writevtk(model, "verification/output/model")
    DC = ["left","right"]

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)

    Γ = BoundaryTriangulation(model,tags=["top", "bottom"])
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)

    # manufactured solutions

    # u(x,t::Real) = VectorValue(sin((π * t)/(2 * Tend)) * (x[1] + x[2]), sin((π * t)/(2 * Tend)) * (x[1] + x[2]))
    # u(t::Real) = x -> u(x,t)

    # h(x,t::Real) = sin((π * t)/(2 * Tend)) * (x[1] + x[2])
    # h(t::Real) = x -> h(x,t)

    # ug(x, t::Real) = VectorValue(sin((π * t)/(2 * Tend)), sin((π * t)/(2 * Tend)))
    # ug(t::Real) = x -> ug(x,t)

    # u(x,t::Real) = VectorValue(sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * sin((π * x[2])/(L)),
    #                              sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * sin((π * x[2])/(L)))
    # u(t::Real) = x -> u(x,t)

    # h(x,t::Real) = sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * sin((π * x[2])/(L))
    # h(t::Real) = x -> h(x,t)

    A = (π/W)
    B = (π/L)
    C = ((4 * π)/(1 * Tend))

    D = (π/W)
    E = (π/(7 * L))
    F = A
    G = B

    yref = L/2

    u(x,t::Real) = VectorValue(0.1 * sin(A*x[1])*sin(B*x[2])*sin(C*t) + 1.0,
                                0.1 * sin(A*x[1])*cos(E*(x[2]-yref))*sin(C*t) + 1.0)
    u(t::Real) = x -> u(x,t)

    h(x,t::Real) = 0.1 * sin(F*x[1])*sin(G*x[2])*sin(C*t) + H
    h(t::Real) = x -> h(x,t)

    u_norm(x,t::Real) = sqrt(0.1 * (sin(A*x[1])*sin(B*x[2])*sin(C*t) + 1.0)^2 + 
                                0.1 * (sin(A*x[1])*cos(E*(x[2]-yref))*sin(C*t) + 1.0)^2)
    u_norm(t::Real) = x -> u_norm(x,t)

    
    #corresponding force functions
    # f_u(x,t::Real) = VectorValue((π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * (x[1] + x[2]) - f * sin((π * t)/(2 * Tend)) * (x[1] + x[2]) + g * sin((π * t)/(2 * Tend)),
    # (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * (x[1] + x[2]) + f * sin((π * t)/(2 * Tend)) * (x[1] + x[2]) + g * sin((π * t)/(2 * Tend)))
    # f_u(t::Real) = x -> f_u(x,t)

    # f_h(x,t::Real) = (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * (x[1] + x[2]) + H * 2 * sin((π * t)/(2 * Tend))
    # f_h(t::Real) = x -> f_h(x,t)

    # f_u(x,t::Real) = VectorValue(((π * g)/B) * sin((π * t)/(2 * Tend)) * cos((π * x[1])/(B)) * sin((π * x[2])/(L)) + 
    #                             (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * sin((π * x[2])/(L)),
    #                             ((π * g)/L) * sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * cos((π * x[2])/(L)) + 
    #                             (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * sin((π * x[2])/(L)))
    # f_u(t::Real) = x -> f_u(x,t)

    # f_h(x,t::Real) = H * ((π/B1) * sin((π * t)/(2 * Tend)) * cos((π * x[1])/(B1)) * sin((π * x[2])/(L)) + 
    #                         (π/L) * sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * cos((π * x[2])/(L))) +
    #                         (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * sin((π * x[2])/(L))
    # f_h(t::Real) = x -> f_h(x,t)

    f_u(x,t::Real) = 0.1 * VectorValue(F*g*cos(F*x[1])*sin(G*x[2])*sin(C*t)+C*sin(A*x[1])*sin(B*x[2])*cos(C*t),
                                G*g*sin(F*x[1])*cos(G*x[2])*sin(C*t)+C*sin(A*x[1])*cos(E*(x[2]-yref))*cos(C*t))
    f_u(t::Real) = x -> f_u(x,t)

    f_h(x,t::Real) =  H*(A*cos(A*x[1])*sin(B*x[2])*sin(C*t)-E*sin(A*x[1])*sin(E*(x[2]-yref))*sin(C*t))+C*sin(F*x[1])*sin(G*x[2])*cos(C*t)
    f_h(t::Real) = x -> f_h(x,t)

    #Make reference spaces

    reffe_rt = ReferenceFE(raviart_thomas, Float64, order)
    V = TestFESpace(model,reffe_rt,conformity=:HDiv,dirichlet_tags=DC)#, dirichlet_masks=[(false, true), (false, true)])
    udc(x,t::Real) = VectorValue(0.0, 0.0)
    udc(t::Real) = x -> udc(x,t)
    U = TransientTrialFESpace(V, [udc, udc])

    reffe_lg = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lg,conformity=:L2)
    hdc(x,t::Real) = 0.0
    hdc(t::Real) = x -> hdc(x,t)
    P = TransientTrialFESpace(Q)

    Y = MultiFieldFESpace([V,Q])
    X = TransientMultiFieldFESpace([U,P])


    res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) + (∇⋅(w))*g*h + ∂t(h)*ϕ + ϕ*H*(∇⋅(u)) )dΩ # - f_u(t)⋅w - f_h(t)*ϕ
    jac(t,(u,h),(du,dh),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dh + ϕ*(H*(∇⋅(du))))dΩ
    jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ

    op = TransientFEOperator(res,jac,jac_t,X,Y)
    #nls = NLSolver(show_trace=true,linesearch=BackTracking())
    nls = LUSolver()

    h_init(x, t::Real) = H + 0.5 * x[1]#0.05*exp(-0.01*(x[1]-0.5)^2 -0.01*(x[2]-2)^2)
    h_init(t::Real) = x -> h_init(x,t)

    u_init(x, t::Real) = VectorValue(0.0, 0.0)
    u_init(t::Real) = x -> u_init(x,t)

    #x0 = interpolate_everywhere([u(t0), h(t0)], X(t0))
    x0 = interpolate_everywhere([u_init(t0), h_init(t0)], X(t0))
    
    ode_solver = ThetaMethod(nls,dt,0.5)
    s_h = solve(ode_solver,op,x0,0.0,Tend)
    

eul2 = Float64[] # L2 norm error array
ehl2 = Float64[]
h_c = Float64[]
u_c = Float64[]
energy = Float64[]

t_base = Float64[] # time base

# (x, t) = s_h(0.0)
# u_h, h_h = x

# initial_energy = 0.5 * (g * sum( ∫( h_h*h_h )*dΩ ) + H * sum( ∫( u_h⋅u_h )*dΩ ))

if isdir(dir_pvd)
    output_file = paraview_collection(joinpath(dir_pvd,"energy"))do pvd
        #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
        for (x,t) in s_h
            u_h,h_h = x
            
            eh = (h_h - h(t))
            eu = (u_h - u(t))
            
            h_grid = (h_h - h(t)) / h(t)
            u_grid = (u_h - u(t)) / u_norm(t)

            h_cell = sum(∫( h_grid*h_grid )*dΩ)
            u_cell = sum(∫( u_grid⋅u_grid )*dΩ)

            push!(ehl2, (sqrt(sum( ∫( eh*eh )*dΩ ))) )
            push!(eul2, (sqrt(sum( ∫( eu⋅eu )*dΩ ))) )
            push!(h_c, h_cell)
            push!(u_c, u_cell)
            push!(energy, (0.5 * (g * sum( ∫( h_h*h_h )*dΩ ) + H * sum( ∫( u_h⋅u_h )*dΩ ))))
            # push!(eh1, sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ )))

            push!(t_base, t)
            pvd[t] = createvtk(Ω,joinpath(dir_pvd,"linear$t.vtu"),cellfields=["u"=>u_h,"h"=>h_h])
                #"hm"=>h(t), "eh"=>eh, "eu"=>eu, "eh_rel"=>h_grid, "eu_rel"=>u_grid, "fu"=>f_u(t), "fh"=>f_h(t), "um"=>u(t)])
            println("done $t/$Tend")
        end
    end
else
    mkdir(dir_pvd)
    output_file = paraview_collection(joinpath(dir_pvd,"energy"))do pvd
        #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
        for (x,t) in s_h
            u_h,h_h = x
            
            eh = (h_h - h(t))
            eu = (u_h - u(t))
            
            h_grid = (h_h - h(t)) / h(t)
            u_grid = (u_h - u(t)) / u_norm(t)

            h_cell = sum(∫( h_grid*h_grid )*dΩ)
            u_cell = sum(∫( u_grid⋅u_grid )*dΩ)

            push!(ehl2, (sqrt(sum( ∫( eh*eh )*dΩ ))) )
            push!(eul2, (sqrt(sum( ∫( eu⋅eu )*dΩ ))) )
            push!(h_c, h_cell)
            push!(u_c, u_cell)
            push!(energy, (0.5 * (g * sum( ∫( h_h*h_h )*dΩ ) + H * sum( ∫( u_h⋅u_h )*dΩ ))))
            # push!(eh1, sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ )))

            push!(t_base, t)
            pvd[t] = createvtk(Ω,joinpath(dir_pvd,"linear$t.vtu"),cellfields=["u"=>u_h,"h"=>h_h])
                #"hm"=>h(t), "eh"=>eh, "eu"=>eu, "eh_rel"=>h_grid, "eu_rel"=>u_grid, "fu"=>f_u(t), "fh"=>f_h(t), "um"=>u(t)])
            println("done $t/$Tend")
        end
    end
end

# p1 = plot(t_base,[eul2, ehl2, h_c],
# label=["L2_u" "L2_h"],
# shape=:auto,
# xlabel="time",ylabel="error L2 norm")
# #savefig(p1, joinpath(dir_plots,"norm.pdf"))

# p2 = plot(t_base,[h_c, u_c],
# label=["h_rel" "u_rel"],
# shape=:auto,
# xlabel="time",ylabel="relative error")
# #savefig(p2, joinpath(dir_plots,"norm_rel.pdf"))

p1 = plot(t_base,[energy],
label=["total_energy"],
shape=:auto,
xlabel="time",ylabel="total energy")
savefig(p1, joinpath(dir_plots,"energy.pdf"))

(eul2, ehl2, t_base)

end


order = 2
degree = 4
Tend = 50
dt = 0.5
stamp = 2

modelnr_def = 100
linear_SWE(order, degree, dt, stamp, modelnr_def)


eul2 = Float64[]
ehl2 = Float64[]
eul2s = Float64[]
ehl2s = Float64[]
t_base = Float64[]
ts = Float64[]

# Time convergence
# time_steps = [1.6, 0.8, 0.4, 0.2, 0.1]
# for t in time_steps

#     stamp = round(Int64, Tend/t)
#     eul2, ehl2, t_base = linear_SWE(order, degree, t, stamp, 50)
    

#     push!(eul2s, sum(eul2)/stamp)
#     push!(ehl2s, sum(ehl2)/stamp)
#     push!(ts, t)

# end

# plot(ts,[eul2s, ehl2s],
# label=["L2_u" "L2_h"],
# shape=:auto,
# xlabel="time step",ylabel="averaged error norm",
# yscale = :log, xscale = :log,
# legend = :topleft)
# savefig(joinpath(dir_plots,"time_conv_log.pdf"))

# print(ts)
# print(eul2s, ehl2s)

# Spatial convergence
# meshes = [1500, 1000, 500, 100, 50, 10, 5]
# x_base = Float64[]

# for m in meshes

#     n = Tend/dt
#     # model_name = "swe-solver/meshes/1x4periodic$m.msh"
#     eul2, ehl2, t_base = linear_SWE(1, 4, dt, m, m)
    
#     push!(eul2s, sum(eul2)/n)
#     push!(ehl2s, sum(ehl2)/n)
#     push!(x_base, m/1000.0)

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


# Bdr degree comparison

# degree = 2
# eul2, ehl2, t_base = linear_SWE(1, degree, dt, degree, model_str_def)
# plot(t_base,[eul2],
# label=["U_$degree"],
# shape=:auto,
# xlabel="time",ylabel="error")
# savefig(joinpath(dir_plots,"bdr_degree_new.pdf"))




