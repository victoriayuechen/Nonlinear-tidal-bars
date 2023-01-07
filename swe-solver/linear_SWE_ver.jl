using Pkg
Pkg.activate(".")

using Gridap
using Gridap.CellData
using WriteVTK
using LineSearches: BackTracking
# using GridapGmsh
# include("mesh_generator.jl")
# using .MyMeshGenerator



function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end

# function uh(u₀,h₀,X,Y,dΩ)
#     a((u,h),(w,ϕ)) = ∫(w⋅u +ϕ*h)dΩ
#     b((w,ϕ)) = ∫(w⋅u₀ + ϕ*h₀)dΩ
#     solve(AffineFEOperator(a,b,X,Y))
# end


function linear_SWE(order,degree)

    #Parameters
    B = 1 #Channel width
    L = 1 #Channel Length
    latitude = 52
    η = 7.29e-5
    #f = 2*η*sin(latitude*(π/180))
    f=0
    g = 9.81
    H = 0.5 #Constant layer depth at rest
    Tend = 50
    t0 = 0
    B1 = 1.5 * B


    #Make model
    domain = (0,B,0,L)
    partition = (20,20)
    # Generate the model
    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false))

    #Make labels
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    writevtk(model, "verification/output/model")
    DC = ["left","right","top","bottom"]
    dir = "verification/output"
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)

    # manufactured solutions

    # u(x,t::Real) = VectorValue(sin((π * t)/(2 * Tend)) * (x[1] + x[2]), sin((π * t)/(2 * Tend)) * (x[1] + x[2]))
    # u(t::Real) = x -> u(x,t)

    # h(x,t::Real) = sin((π * t)/(2 * Tend)) * (x[1] + x[2])
    # h(t::Real) = x -> h(x,t)

    # ug(x, t::Real) = VectorValue(sin((π * t)/(2 * Tend)), sin((π * t)/(2 * Tend)))
    # ug(t::Real) = x -> ug(x,t)

    u(x,t::Real) = VectorValue(sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * sin((π * x[2])/(L)),
                                 sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * sin((π * x[2])/(L)))
    u(t::Real) = x -> u(x,t)

    h(x,t::Real) = sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * sin((π * x[2])/(L))
    h(t::Real) = x -> h(x,t)

    # ug(x, t::Real) = VectorValue(sin((π * t)/(2 * Tend)), sin((π * t)/(2 * Tend)))
    # ug(t::Real) = x -> ug(x,t)

    #corresponding force functions
    # f_u(x,t::Real) = VectorValue((π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * (x[1] + x[2]) - f * sin((π * t)/(2 * Tend)) * (x[1] + x[2]) + g * sin((π * t)/(2 * Tend)),
    # (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * (x[1] + x[2]) + f * sin((π * t)/(2 * Tend)) * (x[1] + x[2]) + g * sin((π * t)/(2 * Tend)))
    # f_u(t::Real) = x -> f_u(x,t)

    # f_h(x,t::Real) = (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * (x[1] + x[2]) + H * 2 * sin((π * t)/(2 * Tend))
    # f_h(t::Real) = x -> f_h(x,t)

    f_u(x,t::Real) = VectorValue(((π * g)/B) * sin((π * t)/(2 * Tend)) * cos((π * x[1])/(B)) * sin((π * x[2])/(L)) + 
                                (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * sin((π * x[2])/(L)),
                                ((π * g)/L) * sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * cos((π * x[2])/(L)) + 
                                (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * sin((π * x[2])/(L)))
    f_u(t::Real) = x -> f_u(x,t)

    f_h(x,t::Real) = H * ((π/B1) * sin((π * t)/(2 * Tend)) * cos((π * x[1])/(B1)) * sin((π * x[2])/(L)) + 
                            (π/L) * sin((π * t)/(2 * Tend)) * sin((π * x[1])/(B1)) * cos((π * x[2])/(L))) +
                            (π/(2 * Tend)) * cos((π * t)/(2 * Tend)) * sin((π * x[1])/(B)) * sin((π * x[2])/(L))
    f_h(t::Real) = x -> f_h(x,t)
    #Make reference spaces

    reffe_rt = ReferenceFE(raviart_thomas,Float64, order)
    V = TestFESpace(model,reffe_rt,conformity=:HDiv,dirichlet_tags=DC)
    udc(x,t::Real) = VectorValue(1.0, 0.0)
    udc(t::Real) = x -> udc(x,t)
    U = TransientTrialFESpace(V, [u, u, u, u])

    reffe_lg = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lg,conformity=:H1, dirichlet_tags=DC)
    hdc(x,t::Real) = 0.0
    hdc(t::Real) = x -> hdc(x,t)
    P = TransientTrialFESpace(Q, [hdc, hdc, hdc, hdc])

    Y = MultiFieldFESpace([V,Q])
    X = TransientMultiFieldFESpace([U,P])

    # a1(u,v) = ∫(v⋅u)dΩ
    # l1(v) = ∫(v⋅u₀)dΩ
    # un = solve(TransientFEOperator(a1,l1,U,V))

    # a2(u,v) = ∫(v*u)dΩ
    # l2(v) = ∫(v*h₀)dΩ
    # hn = solve(TransientFEOperator(a2,l2,P,Q))

    # #Lu_m(x,t::Real) = fu1_res(x, t)
    # Lu_m(x,t::Real) = VectorValue(π*cos(π*t)*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])^2*sin(π*x[2])*π*cos(π*t/2)/2 + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*cos(π*x[2])/2 + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*sin(π*t)/2 - f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*cos(π*x[1]) - f*sin(π*x[1])*sin(π*x[2])^2*sin(π*t/2) - f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t) + g*π*cos(π*x[1]), π*cos(π*t)*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*cos(π*x[1])/2 + sin(π*x[1])*sin(π*x[2])^2*π*cos(π*t/2)/2 + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*sin(π*t)/2 + f*sin(π*x[1])^2*sin(π*x[2])*sin(π*t/2) + f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*cos(π*x[2]) + f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t) + g*π*cos(π*x[2]))
    # Lu_m(t::Real) = x -> Lu_m(x, t)

    # #Lh_m(x, t::Real) = fh_res(x[1], x[2], t)
    # Lh_m(x, t::Real) = π*cos(π*t) + 2*H*π*cos(π*x[1])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + H*π*cos(π*x[1])*sin(π*x[2])*sin(π*t/2)*cos(π*x[2]) + H*π*cos(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t) + 2*H*π*cos(π*x[2])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + H*sin(π*x[1])*π*cos(π*x[2])*sin(π*t/2)*cos(π*x[1]) + H*sin(π*x[1])*π*cos(π*x[2])*sin(π*t/2)*sin(π*t)
    # Lh_m(t::Real) = x -> Lh_m(x, t)

    # u_m(x, t::Real) = VectorValue(sin(π*x[1])^2*sin(π*x[2])*sin(π*t/2) + cos(π*x[2])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t), cos(π*x[1])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])^2*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t))
    # u_m(t::Real) = x -> u_m(x, t)

    # h_m(x, t::Real) = 0.5 + sin(π*x[1]) #+ sin(π*x[2]) + sin(π*t)
    # h_m(t::Real) = x -> h_m(x, t)

    # uhn = uh(un,hn,X,Y,dΩ)
    # un,hn = uhn
    #forcefunc(t) = VectorValue(0.0,0.0*0.5*cos(π*(1/5)*t))

    #res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) - (∇⋅(w))*g*h + ∂t(h)*ϕ + ϕ*H*(∇⋅(u)) - f_u(t)⋅w - f_h(t)*ϕ)dΩ #- Lu_m(t)⋅w - Lh_m(t)*ϕ
    res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) + (∇(h) ⋅ w)*g + ∂t(h)*ϕ - (∇(ϕ)⋅(u*H)) - f_u(t)⋅w - f_h(t)*ϕ)dΩ # + f_u(t)⋅w + f_h(t)*ϕ
    #jac(t,(u,h),(du,dh),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dh + ϕ*(H*(∇⋅(du))))dΩ
    #jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ

    op = TransientFEOperator(res,X,Y)
    #nls = NLSolver(show_trace=true,linesearch=BackTracking())
    nls = LUSolver()

    #h_init(x) = 0.5 + 0.05*exp(-0.01*(x[1]-50)^2 -0.01*(x[2]-25)^2)
    #h_init(t::Real) = x -> h_init(x,t)

    x0 = interpolate_everywhere([u(t0), h(t0)], X(t0))
    
    ode_solver = ThetaMethod(nls,2,0.5)
    s_h = solve(ode_solver,op,x0,0.0,Tend)


    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
            for (x,t) in s_h
                u_h,h_h = x
                eh = h_h - h(t)
                eu = u_h - u(t)
                pvd[t] = createvtk(Ω,joinpath(dir,"linear$t.vtu"),cellfields=["u"=>u_h,"h"=>h_h, "fu"=>f_u(t), "fh"=>f_h(t), "um"=>u(t),
                 "hm"=>h(t), "eh"=>eh, "eu"=>eu])#, "um"=>u(t), "hm"=>h(t) , "eh"=>eh, "eu"=>eu,
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear")) do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
            for (x,t) in s_h
                u_h,h_h = x
                eh = h_h - h(t)
                eu = u_h - u(t)
                pvd[t] = createvtk(Ω,joinpath(dir,"linear$t.vtu"),cellfields=["u"=>u_h,"h"=>h_h, "eh"=>eh, "eu"=>eu, "um"=>u(t), "hm"=>h(t)])#, "um"=>u(t), "hm"=>h(t)
                println("done $t/$Tend")
            end
        end
    end
end

# function h₀((x,y))
#     h = 0.5 + 0.05*exp(-0.01*(x-50)^2 -0.01*(y-25)^2)
#     h
# end


# function u₀((x,y))
#     u = VectorValue(0.0,0.0)
#     u
# end

# function fu1_res((x, y), t)
#     H_const = H
#     u_res = VectorValue(π*cos(π*t) - 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y) - H_const*π*cos(π*x)*sin(π*y)*exp(t)*cos(π*y) + H_const*π*cos(π*x)*sin(π*y)*cos(π*y) - H_const*π*cos(π*x)*sin(π*y)*sin(π*t)*exp(t) + H_const*π*cos(π*x)*sin(π*y)*sin(π*t) - 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y) - H_const*sin(π*x)*π*cos(π*y)*exp(t)*cos(π*x) + H_const*sin(π*x)*π*cos(π*y)*cos(π*x) - H_const*sin(π*x)*π*cos(π*y)*sin(π*t)*exp(t) + H_const*sin(π*x)*π*cos(π*y)*sin(π*t), S2 = -π*cos(π*t)*sin(π*x)*sin(π*y)*exp(t) + π*cos(π*t)*sin(π*x)*sin(π*y) - sin(π*x)*sin(π*y)*exp(t)*cos(π*x) - sin(π*x)*sin(π*y)^2*exp(t) - sin(π*x)*sin(π*y)*exp(t)*sin(π*t) - f*sin(π*x)^2*sin(π*y)*exp(t) + f*sin(π*x)^2*sin(π*y) - f*sin(π*x)*sin(π*y)*exp(t)*cos(π*y) + f*sin(π*x)*sin(π*y)*cos(π*y) - f*sin(π*x)*sin(π*y)*sin(π*t)*exp(t) + f*sin(π*x)*sin(π*y)*sin(π*t) + g*π*cos(π*y))
#     u_res
# end

# function fh_res(x, y, t)
#     H_const = H
#     h_res = π*cos(π*t) - 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y) - H_const*π*cos(π*x)*sin(π*y)*exp(t)*cos(π*y) + H_const*π*cos(π*x)*sin(π*y)*cos(π*y) - H_const*π*cos(π*x)*sin(π*y)*sin(π*t)*exp(t) + H_const*π*cos(π*x)*sin(π*y)*sin(π*t) - 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y) - H_const*sin(π*x)*π*cos(π*y)*exp(t)*cos(π*x) + H_const*sin(π*x)*π*cos(π*y)*cos(π*x) - H_const*sin(π*x)*π*cos(π*y)*sin(π*t)*exp(t) + H_const*sin(π*x)*π*cos(π*y)*sin(π*t)
#     h_res
# end

# function u_man(t, (x, y))
#     u_man = VectorValue((sin(π*x) + cos(π*y) + sin(π*t))*sin(π*x)*sin(π*y)*(-exp(t) + 1))

# end
linear_SWE(1,2)