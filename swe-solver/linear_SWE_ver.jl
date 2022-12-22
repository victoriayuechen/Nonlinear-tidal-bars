using Pkg
Pkg.activate(".")

using Gridap
using WriteVTK
using LineSearches: BackTracking
# using GridapGmsh
# include("mesh_generator.jl")
# using .MyMeshGenerator



function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end

function uh(u₀,h₀,X,Y,dΩ)
    a((u,h),(w,ϕ)) = ∫(w⋅u +ϕ*h)dΩ
    b((w,ϕ)) = ∫(w⋅u₀ + ϕ*h₀)dΩ
    solve(AffineFEOperator(a,b,X,Y))
end


function linear_SWE(order,degree,h₀,u₀)

    #Parameters
    B = 100 #Channel width
    L = 100 #Channel Length
    latitude = 52
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    g = 9.81
    H = 0.5 #Constant layer depth at rest


    #Make model
    domain = (0,B,0,L)
    partition = (50,50)
    # Generate the model
    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,true))

    #Make labels
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["left","right"]
    dir = "swe-solver/output_linear_swe"
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)



    #Make reference spaces

    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
    V = TestFESpace(model,reffe_rt,conformity=:HDiv,dirichlet_tags=DC)
    udc(x,t::Real) = VectorValue(0.0,0.0)
    udc = x -> udc(x,t)
    U = TransientTrialFESpace(V)

    reffe_lg = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lg,conformity=:L2)
    P = TransientTrialFESpace(Q)

    Y = MultiFieldFESpace([V,Q])
    X = MultiFieldFESpace([U,P])

    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))

    #Lu_m(x,t::Real) = fu1_res(x, t)
    Lu_m(x,t::Real) = VectorValue(π*cos(π*t)*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])^2*sin(π*x[2])*π*cos(π*t/2)/2 + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*cos(π*x[2])/2 + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*sin(π*t)/2 - f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*cos(π*x[1]) - f*sin(π*x[1])*sin(π*x[2])^2*sin(π*t/2) - f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t) + g*π*cos(π*x[1]), π*cos(π*t)*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*cos(π*x[1])/2 + sin(π*x[1])*sin(π*x[2])^2*π*cos(π*t/2)/2 + sin(π*x[1])*sin(π*x[2])*π*cos(π*t/2)*sin(π*t)/2 + f*sin(π*x[1])^2*sin(π*x[2])*sin(π*t/2) + f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*cos(π*x[2]) + f*sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t) + g*π*cos(π*x[2]))
    Lu_m(t::Real) = x -> Lu_m(x, t)

    #Lh_m(x, t::Real) = fh_res(x[1], x[2], t)
    Lh_m(x, t::Real) = π*cos(π*t) + 2*H*π*cos(π*x[1])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + H*π*cos(π*x[1])*sin(π*x[2])*sin(π*t/2)*cos(π*x[2]) + H*π*cos(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t) + 2*H*π*cos(π*x[2])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + H*sin(π*x[1])*π*cos(π*x[2])*sin(π*t/2)*cos(π*x[1]) + H*sin(π*x[1])*π*cos(π*x[2])*sin(π*t/2)*sin(π*t)
    Lh_m(t::Real) = x -> Lh_m(x, t)

    u_m(x, t::Real) = VectorValue(sin(π*x[1])^2*sin(π*x[2])*sin(π*t/2) + cos(π*x[2])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t), cos(π*x[1])*sin(π*x[1])*sin(π*x[2])*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])^2*sin(π*t/2) + sin(π*x[1])*sin(π*x[2])*sin(π*t/2)*sin(π*t))
    u_m(t::Real) = x -> u_m(x, t)

    h_m(x, t::Real) = 0.5 + sin(π*x[1]) #+ sin(π*x[2]) + sin(π*t)
    h_m(t::Real) = x -> h_m(x, t)

    uhn = uh(un,hn,X,Y,dΩ)
    un,hn = uhn
    #forcefunc(t) = VectorValue(0.0,0.0*0.5*cos(π*(1/5)*t))

    res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) - (∇⋅(w))*g*h + ∂t(h)*ϕ + ϕ*H*(∇⋅(u)))dΩ #- Lu_m(t)⋅w - Lh_m(t)*ϕ
    jac(t,(u,h),(du,dh),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dh + ϕ*(H*(∇⋅(du))))dΩ
    jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ

    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,linesearch=BackTracking())
    Tend = 200
    ode_solver = ThetaMethod(nls,2,0.5)
    x = solve(ode_solver,op,uhn,0.0,Tend)

    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_topo"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
            for (x,t) in x
                u,h = x
                eh = h - h_m(t)
                eu = u - u_m(t)
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"h"=>(h), "eh"=>eh, "eu"=>eu, "um"=>u_m(t), "hm"=>h_m(t)])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_topo")) do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn)])
            for (x,t) in x
                u,h = x
                eh = abs(h - h_m(t))
                eu = abs(u - u_m(t))
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"h"=>(h), "eh"=>eh, "eu"=>eu, "um"=>u_m(t), "hm"=>h_m(t)])
                println("done $t/$Tend")
            end
        end
    end
end

function h₀((x,y))
    h = 0.5 + 0.05*exp(-0.01*(x-50)^2 -0.01*(y-25)^2)
    h
end


function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function fu1_res((x, y), t)
    H_const = H
    u_res = VectorValue(π*cos(π*t) - 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y) - H_const*π*cos(π*x)*sin(π*y)*exp(t)*cos(π*y) + H_const*π*cos(π*x)*sin(π*y)*cos(π*y) - H_const*π*cos(π*x)*sin(π*y)*sin(π*t)*exp(t) + H_const*π*cos(π*x)*sin(π*y)*sin(π*t) - 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y) - H_const*sin(π*x)*π*cos(π*y)*exp(t)*cos(π*x) + H_const*sin(π*x)*π*cos(π*y)*cos(π*x) - H_const*sin(π*x)*π*cos(π*y)*sin(π*t)*exp(t) + H_const*sin(π*x)*π*cos(π*y)*sin(π*t), S2 = -π*cos(π*t)*sin(π*x)*sin(π*y)*exp(t) + π*cos(π*t)*sin(π*x)*sin(π*y) - sin(π*x)*sin(π*y)*exp(t)*cos(π*x) - sin(π*x)*sin(π*y)^2*exp(t) - sin(π*x)*sin(π*y)*exp(t)*sin(π*t) - f*sin(π*x)^2*sin(π*y)*exp(t) + f*sin(π*x)^2*sin(π*y) - f*sin(π*x)*sin(π*y)*exp(t)*cos(π*y) + f*sin(π*x)*sin(π*y)*cos(π*y) - f*sin(π*x)*sin(π*y)*sin(π*t)*exp(t) + f*sin(π*x)*sin(π*y)*sin(π*t) + g*π*cos(π*y))
    u_res
end

function fh_res(x, y, t)
    H_const = H
    h_res = π*cos(π*t) - 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y) - H_const*π*cos(π*x)*sin(π*y)*exp(t)*cos(π*y) + H_const*π*cos(π*x)*sin(π*y)*cos(π*y) - H_const*π*cos(π*x)*sin(π*y)*sin(π*t)*exp(t) + H_const*π*cos(π*x)*sin(π*y)*sin(π*t) - 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y) - H_const*sin(π*x)*π*cos(π*y)*exp(t)*cos(π*x) + H_const*sin(π*x)*π*cos(π*y)*cos(π*x) - H_const*sin(π*x)*π*cos(π*y)*sin(π*t)*exp(t) + H_const*sin(π*x)*π*cos(π*y)*sin(π*t)
    h_res
end

function u_man(t, (x, y))
    u_man = VectorValue((sin(π*x) + cos(π*y) + sin(π*t))*sin(π*x)*sin(π*y)*(-exp(t) + 1))

end
linear_SWE(0,3,h₀,u₀)