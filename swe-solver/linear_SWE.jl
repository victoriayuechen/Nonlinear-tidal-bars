using Gridap
using WriteVTK
using LineSearches: BackTracking
using GridapGmsh
# include("mesh_generator.jl")
# using .MyMeshGenerator

export run_linear_SWE
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


function run_linear_SWE(order,degree,ζ₀,u₀,forcefunc,Tend,dt,model,H,DC,dir,latitude)

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
    uda(x,t::Real) = VectorValue(1*sin(π*(1/3600)*t),0.0)
    udc(t::Real) = x -> udc(x,t)
    uda(t::Real) = x -> uda(x,t)
    U = TransientTrialFESpace(V,[udc,udc,udc,uda])

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
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,x0,0.0,Tend)

    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_SWE"))do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_SWE0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_SWE$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ)])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_SWE")) do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_SWE0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_SWE$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ)])
                println("done $t/$Tend")
            end
        end
    end
end

function ζ₀((x,y))
    h =0.0*exp(-0.1*(x-50)^2 -0.1*(y-30)^2)
    h
end


function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function forcefunc((x,y),t)
    f = VectorValue(0.0,0.0*cos(π*(1/5)*t))
    f
end

order = 1
degree = 4
model = GmshDiscreteModel("swe-solver/meshes/schelde_mesh.msh")
DC = ["shore_top","shore_bot","island","outlet_left"]
dir = "swe-solver/output_linear_swe"
H = 0.5
latitude = 52
Tend = 60*60*6
dt = 120


run_linear_SWE(order,degree,ζ₀,u₀,forcefunc,Tend,dt,model,H,DC,dir,latitude)