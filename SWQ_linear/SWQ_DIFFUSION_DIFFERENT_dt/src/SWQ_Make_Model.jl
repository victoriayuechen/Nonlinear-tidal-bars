using Gridap
using Parameters
includet("SWQ_Parameters.jl")


function Make_model(Param::Parameter)
    #''''''''''''''Make model''''''''''''''##
    @unpack B, L, x_points, y_points, order, degree, periodic_x, periodic_y, dirichlet_mask_x_1, dirichlet_mask_y_1, dirichlet_mask_x_2, dirichlet_mask_y_2 = Param
    domain = (0,L,0,B)                      #Original x and y flipped
    partition = (x_points,y_points)         #Partition

    ##''''''''''''''Generate the model''''''''''''''##
    model = CartesianDiscreteModel(domain,partition;isperiodic=(periodic_x,periodic_y))  

    ##''''''''''''''Make labels''''''''''''''##
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["bottom","top"]
    
    ##''''''''''''''Define triangulation''''''''''''''##
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)
  
    ##''''''''''''''Define function spaces''''''''''''''##
    v_boundary(x) = VectorValue(0.0,0.0)
    reffe_rt = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    V = TestFESpace(model,reffe_rt,dirichlet_tags=DC,dirichlet_masks=[(dirichlet_mask_x_1,dirichlet_mask_y_1),(dirichlet_mask_x_2,dirichlet_mask_y_2)])  #Zero at the x-boundaries
    U = TransientTrialFESpace(V)
    
    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn)
    P = TransientTrialFESpace(Q) 

    Y = MultiFieldFESpace([V,Q])  #u, ζ
    X = TransientMultiFieldFESpace([U,P])
    return Ω, dΩ, dΓ, Y, X, P, nΓ
end