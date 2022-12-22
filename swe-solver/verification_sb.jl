using Pkg
Pkg.activate(".")

using Gridap
using WriteVTK
using LineSearches: BackTracking


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

    h_m(x, t::Real) = 0.5 + sin(π*x[1]) + sin(π*x[2]) + sin(π*t)
    h_m(t::Real) = x -> h_m(x, t)

    uhn = uh(un,hn,X,Y,dΩ)
    un,hn = uhn
