using Gridap
using Parameters
includet("SWQ_Parameters_APVM.jl")


function Make_model(Param::Parameter)
    #''''''''''''''Make model''''''''''''''##
    @unpack B, L, x_points, y_points, order, degree = Param
    domain = (0,L,0,B)                      #Original x and y flipped
    partition = (x_points,y_points)         #Partition

    ##''''''''''''''Generate the model''''''''''''''##
    model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))  

    ##''''''''''''''Make labels''''''''''''''##
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["bottom","top"]
    
    ##''''''''''''''Define triangulation''''''''''''''##
    #Make triangulations and boundaries
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    periodic = true
    if periodic
        Γ = BoundaryTriangulation(model,tags=DC)
    else
        Γ = BoundaryTriangulation(model,tags=DC)
    end
    global nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)
  
    ##''''''''''''''Define function spaces''''''''''''''##
    udc(x,t::Real) = VectorValue(0.0,0.0)
    udc(t::Real) = x -> udc(x,t)


    #Make FE spaces
    if periodic
        reffe_rt = ReferenceFE(raviart_thomas,Float64,order)#
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)#
        U = TransientTrialFESpace(V,[udc,udc])
    else
        reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)#
        U = TransientTrialFESpace(V,[udc,udc])
    end

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn;conformity=:L2)#
    P = TransientTrialFESpace(Q)

    reffe_lgn = ReferenceFE(lagrangian, Float64, order+1)
    S = TestFESpace(model, reffe_lgn;conformity=:H1)#
    R = TransientTrialFESpace(S)


    Y = MultiFieldFESpace([V,Q,S,V])
    X = TransientMultiFieldFESpace([U,P,R,U])
    return Ω, dΩ, dΓ, Y, X, P, R, S, U, V
end