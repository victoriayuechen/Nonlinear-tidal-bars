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

    ##''''''''''''''Initial_solutions''''''''''''''##
    u_x::Float64 = 0.0
    u_y::Float64 = 0.0
    ζₙ::Float64 = 0.0

    ##''''''''''''''Timestepping''''''''''''''##
    dt1::Real = 5 #10
    dt2::Real = 10
    dt3::Real = 10
    Tstart::Real = 0.0
    Tend1::Real = 44700
    Tend2::Real = 30000
    Tend3::Real = 40000
    theta1::Float64 = 0.5
    theta2::Float64 = 0.5
    theta3::Float64 = 0.5
    T_save::Real = 50

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

    ##''''''''''''''Saving directory''''''''''''''##
    dir = "./test"
end


