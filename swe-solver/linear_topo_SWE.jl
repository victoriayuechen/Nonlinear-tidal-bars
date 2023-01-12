using Pkg
Pkg.activate(".")

using Gridap
using WriteVTK
using LineSearches: BackTracking


function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end


function linear_SWE(order,degree,ζ₀,u₀,topography,forcefunc,Tend,dt,model,H,DC,dir,latitude)
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    g = 9.81

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
    X = TransientMultiFieldFESpace([U,P])

    b = interpolate_everywhere(topography,P(0.0))
    x0 = interpolate_everywhere([u₀,ζ₀],X(0.0))
    un, ζn = x0

    forcefunc_x(t) = x -> forcefunc(x,t)

    res(t,(u,ζ),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅(f*(perp∘(u))) - (∇⋅(w))*g*ζ + ∂t(ζ)*ϕ + ϕ*(H-b)*(∇⋅(u))- ∇(b)⋅u - w⋅forcefunc_x(t))dΩ
    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dζ + ϕ*((H-b)*(∇⋅(du))) - ∇(b)⋅u)dΩ
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ

    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = LUSolver()
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,x0,0.0,Tend)

    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_topo"))do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ)])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"linear_topo")) do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"linear_topo0.0.vtu"),cellfields=["u"=>un,"ζ"=>(ζn)])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"ζ"=>(ζ)])
                println("done $t/$Tend")
            end
        end
    end
end

function h₀((x,y))
    h = 0.05*exp(-0.01*(x-50)^2 -0.01*(y-25)^2)
    h
end

function forcefunc((x,y),t)
    func = VectorValue(0.0,0.0)
    func
end
function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function topography((x,y))
    b = 0.4*exp(-0.01*(y-50)^2)
    b
end

Tend = 50
dt = 0.5
order = 1
degree = 4
B = 100
L = 100
partition = (50,50)
model = CartesianDiscreteModel((0,B,0,L),partition;isperiodic=(false,true))
#Make labels
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])
DC = ["left","right"]
dir = "swe-solver/output_linear_topo_swe"
H = 0.5
latitude = 52

linear_SWE(order,degree,h₀,u₀,topography,forcefunc,Tend,dt,model,H,DC,dir,latitude)