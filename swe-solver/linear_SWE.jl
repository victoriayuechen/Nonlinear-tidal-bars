module Linear_SWE_solver
using Gridap
using WriteVTK
using LineSearches: BackTracking
using GridapGmsh
# include("mesh_generator.jl")
# using .MyMeshGenerator

export run_linear_SWE

function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end

function uζ(u₀,ζ₀,X,Y,dΩ)
    a((u,ζ),(w,ϕ)) = ∫(w⋅u +ϕ*ζ)dΩ
    b((w,ϕ)) = ∫(w⋅u₀ + ϕ*ζ₀)dΩ
    solve(AffineFEOperator(a,b,X,Y))
end


function run_linear_SWE(order,degree,ζ₀,u₀,model)

    #Parameters
    latitude = 52
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    g = 9.81
    H = 0.5 #Constant layer depth at rest

    # Generate the model
    DC = ["left","right"]
    dir = "swe-solver/output_linear_swe"
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    Γ = BoundaryTriangulation(model,tags=["top","bottom"])
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)



    #Make reference spaces

    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)#
    V = TestFESpace(model,reffe_rt,conformity=:HDiv,dirichlet_tags=DC)#
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
    l2(v) = ∫(v*ζ₀)dΩ
    ζn = solve(AffineFEOperator(a2,l2,P,Q))

    uζn = uζ(un,ζn,X,Y,dΩ)
    un,ζn = uζn
    #forcefunc(t) = VectorValue(0.0,0.0*0.5*cos(π*(1/5)*t))

    res(t,(u,ζ),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) - (∇⋅(w))*g*ζ + ∂t(ζ)*ϕ + ϕ*H*(∇⋅(u)))dΩ
    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dζ + ϕ*(H*(∇⋅(du))))dΩ
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ

    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,linesearch=BackTracking())
    Tend = 50
    ode_solver = ThetaMethod(nls,0.5,0.5)
    x = solve(ode_solver,op,uζn,0.0,Tend)

    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_topo"))do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ+H)])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_topo")) do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ+H)])
                println("done $t/$Tend")
            end
        end
    end
end
end