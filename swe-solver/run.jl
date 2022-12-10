using Pkg
Pkg.activate(".")

using Gridap

#Parameters

B = 1000    #[m]
H = 3       #[m]
L = 10000   #[m]
dy = 50     #[m]
dx = 100    #[m]
dt = 5      #[s]
g = 9.81    #[m/s^2]
cd = 0.0025 #[-]




function Shallow_water_theta_newton(
    order,degree,h₀,u₀,f₀,topography)
    #Create model
    domain = (B,L)
    partition = (B/dy,L/dx)

    model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))

    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    nuemanntaggs = ["left","right"]


    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    Γ = BoundaryTriangulation(model,tags=nuemanntaggs)
    dΓ = Measure(Γ,degree)

    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
    V = FESpace(model,reffe_rt;conformity=:HDiv)
    U = TrailFESpace(V)

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = FESpace(model,reffe_lgn;conformity=:L2)
    P = TrailFESpace(Q)

    # reffe_lgn = ReferenceFE(lagrangian,Float64,order+1)
    # S = FESpace(model,reffe_lgn;conformity=:H1)
    # R = TrailFESpace(S)

    X = MultiFieldFESpace([V,Q,V])#∇u, ∇h, Du
    Y = MultiFieldFESpace([U,P,U])

    b = interpolate_everywhere(topography,P)
    E = [0 -1; 1 0]
    #Create initial solutions
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2 = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))




    for step=1:N
        function residual((Δu,Δh,F),(v,q,v2))
            one_m_θ = (1-θ)
            uiΔu = un + one_m_θ*Δu
            hiΔh = hn + one_m_θ*Δh
            hbiΔh = hn + b + one_m_θ*Δh
    
end
