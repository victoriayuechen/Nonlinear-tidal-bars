module MyLinearSWE
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


function run_linear_SWE(order,degree,ζ₀,u₀,forcefunc,Tend,dt,model,H,DC,dir,latitude,filename,tcapture)

    #Parameters 
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    g = 9.81

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())


    #Make reference spaces
    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)#
    V = TestFESpace(model,reffe_rt,dirichlet_tags=DC,conformity=:HDiv)#
    udc(x,t::Real) = VectorValue(0.0,0.0)
    udc(t::Real) = x -> udc(x,t)
    U = TransientTrialFESpace(V,[udc,udc])

    reffe_lg = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lg,conformity=:L2)
    P = TransientTrialFESpace(Q)

    Y = MultiFieldFESpace([V,Q])
    X = TransientMultiFieldFESpace([U,P])


    x0 = interpolate_everywhere([u₀,ζ₀],X(0.0))
    un,ζn = x0

    forcefunc_x(t) = x -> forcefunc(x,t)

    res(t,(u,ζ),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) - (∇⋅(w))*g*ζ + ∂t(ζ)*ϕ + ϕ*H*(∇⋅(u)) - forcefunc_x(t)⋅w)dΩ

    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dζ + ϕ*(H*(∇⋅(du))))dΩ
    
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ


    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = LUSolver()
    ode_solver = ThetaMethod(nls,dt,0.5)#RungeKutta(nls,dt,Symbol("BE_1_0_1"))
    x = solve(ode_solver,op,x0,0.0,Tend)

    if isdir(dir) 
        output_file = paraview_collection(joinpath(dir,"linear_SWE_$(filename)"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_SWE_$(filename)_0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                if (round(t % tcapture,digits=1) == 0.0 || round(t % tcapture,digits=1) == 1.0)
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"linear_SWE_$(filename)_$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ)])
                    println("done $t/$Tend")
                end
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_SWE_$(filename)")) do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_SWE_$(filename)_0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                if (round(t % tcapture,digits=1) == 0.0 || round(t % tcapture,digits=1) == 1.0)
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"linear_SWE_$(filename)_$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ)])
                    println("done $t/$Tend")
                end
            end
        end
    end
end
end

