using Parameters

@with_kw struct Parameter
    ##''''''''''''''Domain''''''''''''''##
    L::Float64 = 10000.0
    B::Float64 = 1000.0
    H::Float64 = 3.0

    ##''''''''''''''Grid''''''''''''''##
    x_points::Int = 100
    y_points::Int = 20
    order::Int = 1
    degree::Int = 4

    ##''''''''''''''Boundary Condition''''''''''''''##
    periodic_x::Bool = true
    periodic_y::Bool = false
    dirichlet_mask_x_1::Bool = false
    dirichlet_mask_y_1::Bool = true
    dirichlet_mask_x_2::Bool = false
    dirichlet_mask_y_2::Bool = true

    ##''''''''''''''Initial_solutions''''''''''''''##
    u_x::Float64 = 0.0
    u_y::Float64 = 0.0
    ζₙ::Float64 = 0.0

    ##''''''''''''''Timestepping''''''''''''''##
    dt::Real = 5 
    dt_spinup::Real = 10
    Tstart::Real = 0.0
    Tend_spinup::Real = 44700 
    Tend::Real = 44700
    theta::Float64 = 0.5
    T_save::Real = 50
    tolerance::Float64 = 2e-7

    ##''''''''''''''Physical parameter''''''''''''''##
    η::Float64 = 7.29e-5                             #angular speed of Earth rotation        (s^(-1))
    latitude::Float64 = 50.0
    f::Float64 = 2*η*sin(latitude*(π/180))           #coriolis parameter                     (s^(-1))
    g::Float64 = 9.81                                #Gravitational constant                 (ms^(-2))
    U_start::Float64 = 0.5                           #Background current amplitude           (ms^(-1))
    σ::Float64 = 2*pi/44700                          #Tidal frequency                        (s^(-1))
    cD::Float64 = 0.0025                             #Drag coefficient                       ()

    ##''''''''''''''Stabilization Parameters''''''''''''''##
    α::Float64 = 1e-6                                #Based on ν
    ν::Float64 = 1e-6                                #From Anna Louka

    ##''''''''''''''Saving directory and name''''''''''''''##
    dir::String = "./test"
    name::String = "test"
    CSVname::String = "CSV.csv"


    ##''''''''''''''Saving''''''''''''''##
    save_CSV::Bool = true                                 #Want to save to CSV? Set true
    nx_start::Real = 0                                    #Start point for savind CSV y-direction
    ny_start::Real = 25                                   #Start point for savind CSV y-direction
    nx::Real = 100                                        #Number of points for saving CSV x-direction
    ny::Real = 20                                         #Number of points for saving CSV y-direction


    ##''''''''''''''control parameters''''''''''''''##
    convection_on::Bool = true                      #Convection term 
    coriolis_on::Bool = true                        #Coriolis term
    friction_on::Bool = true                        #Friction term
    linear_friction_on::Bool = false                #Linear Friction term
    gravitional_on::Bool = true                     #Gravitational term
    forcing_on::Bool = true                         #Forcing function on
    momentum_on::Bool = true                        #Momentum equation on
    linear_on::Bool = false                         #Linear momentum equation
    boundary_on::Bool = true                        #Boundary term
    stabilization_ζ::Bool = true                    #Stabilization term ζ
    stabilization_u::Bool = true                    #Stabilization term u 
    show_iterations::Bool = false                   #Show the iterations
    
end


