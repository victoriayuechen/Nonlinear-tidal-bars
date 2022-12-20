function APVM_run(order,degree,h₀,u₀,topography,dir,periodic::Bool,Tend,dt)
    #Parameters
    B = 100
    L = 100
    partition = (50,50) #Split into X by X  cells
    latitude = 52 #Latitude of the model being analysed
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    cd = 0.025
    g = 9.81
    T0 = 0.0
    τ = dt*0.5

    #model = GmshDiscreteModel("swe-solver/meshes/10x10periodic2.msh")
    DC = ["left","right"]
    BC ="boundary"
    # #Make model
    domain = (0,B,0,L)
    if periodic
        model = CartesianDiscreteModel(domain,partition;isperiodic=(false,true))
    else
        model = CartesianDiscreteModel(domain,partition)
    end

    #Make labels
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["left","right"]

    #Make triangulations and boundaries
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    if periodic
        Γ = BoundaryTriangulation(model,tags=DC)
    else
        Γ = BoundaryTriangulation(model,tags=DC)
    end
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)
    
    #Make FE spaces
    if periodic
        reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)
        U = TransientTrialFESpace(V)
    else
        reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)
        U = TransientTrialFESpace(V)
    end

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn;conformity=:L2)
    P = TransientTrialFESpace(Q)

    reffe_lgn = ReferenceFE(lagrangian, Float64, order+1)
    S = TestFESpace(model, reffe_lgn;conformity=:H1)
    R = TransientTrialFESpace(S)

    H1MM, RTMM, L2MM, H1MMchol, RTMMchol, L2MMchol = setup_and_factorize_mass_matrices(dΩ,R,S,U,V,P,Q)

    Y = MultiFieldFESpace([V,Q,S,V])
    X = MultiFieldFESpace([U,P,R,U])

    #Create initial solutions
    #u, velocity
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))
    #h, fluid depth
    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))
    #b, topography
    a3(u,v) = ∫(v*u)dΩ
    l3(v) = ∫(v*topography)dΩ
    b = solve(AffineFEOperator(a3,l3,P,Q))
    #Constant FE space with the coriolis parameter, possibly can be removed.
    a4(u,v)=∫(v*u)dΩ
    l4(v)=∫(v*f)dΩ
    fn=solve(AffineFEOperator(a4,l4,R,S))

    #Build and compute initial potential vorticity
    q₀ = clone_fe_function(P,fn)
    compute_potential_vorticity!(q₀,H1MM,H1MMchol,dΩ,R,S,hn,un,fn)

    #Build and compute initial mass flux
    F₀ = clone_fe_function(V,un)
    compute_mass_flux!(F₀,dΩ,V,RTMMchol,un*hn)
    
    #Solve whole coupled system with the initial values
    uhn = uh(un,hn,F₀,q₀,X,Y,dΩ)
    un,hn,q,F= uhn

    #Forcing function on u(t)
    forcfunc(t) = VectorValue(0.0,0.25*cos((1/5)*π*t))#Has to be changed correct forcing function at some point

    #Norm operator
    norm(u) = meas∘(u) + 1e-14
    dnorm(u,du) = u ⋅ du / norm(u)

    #Define residual, jacobian in space and time
    res(t,(u,h,q,F),(w,ϕ,ϕ2,w2)) = ∫(∂t(u)⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(F)))⋅w - (∇⋅(w))*(g*(h+b) + 0.5*(u⋅u)) + ∂t(h)*ϕ + w2⋅(F-u*h) + ϕ2*(q*h) +  (perp∘(∇(ϕ2)))⋅u - (ϕ2*fn) + ϕ*(∇⋅(F))- forcfunc(t)⋅w + ((cd*(norm(u))*u)/(h+b+1e-14))⋅w)dΩ + ∫((g*(h+b) + 0.5*(u⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(u))*ϕ2)dΓ 
    jac(t,(u,h,q,F),(du,dh,dq,dF),(w,ϕ,ϕ2,w2)) = ∫(((dq - τ*(u⋅∇(dq) + du⋅∇(q)))*(perp∘(F)))⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(dF)))⋅w - (∇⋅(w))*(g*(dh) + (du⋅u)) + w2⋅(dF - du*h - dh*u) + ϕ2*(q*dh) + ϕ2*(dq*h) + du⋅(perp∘(∇(ϕ2))) + ϕ*(∇⋅(dF)) + (cd*(dnorm(u,du)*u)/(h+b+1e-14) + cd*(norm(u)/(h+b+1e-14)*du) + cd*(norm(u)*u*h)/((h+b+1e-14)*(h+b+1e-14)))⋅w )dΩ+ ∫((g*(dh) + (du⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(du))*ϕ2)dΓ
    jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ

    #Define operators and solvers
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,linesearch=BackTracking())
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,uhn,T0,Tend)

    #Output results
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonlinear_topo"))do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"nonlinear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn+b),"b"=>b])
            for (x,t) in x
                u,h,F = x
                pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"h"=>(h+b),"b"=>b])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonlinear_topo")) do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"nonlinear_topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn+b),"b"=>b])
            for (x,t) in x
                u,h,F = x
                pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"h"=>(h+b),"b"=>b])
                println("done $t/$Tend")
            end
        end
    end
end
end


# #Variable functions to be used to setup model, used for local tests
# function h₀((x,y))
#     hout = -topography((x,y)) +  0.5 + 0.05*exp(-0.01*(x-50)^2 -0.01*(y-25)^2)
#     hout
# end

# function topography((x,y))
#     p = 0.0
#     p
# end

# function u₀((x,y))
#     u = VectorValue(0.0,0.0)
#     u
# end
