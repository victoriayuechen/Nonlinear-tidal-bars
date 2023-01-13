using Parameters

@with_kw struct Parameter
    ##''''''''''''''Domain''''''''''''''##
    L::Float64 = 10000.0
    B::Float64 = 1000.0
    H::Float64 = 3.0

    ##''''''''''''''Grid''''''''''''''##
    x_points = 20
    y_points = 100
    order = 1
    degree = 3

    ##''''''''''''''Timestepping''''''''''''''##
    dt::Real = 5
    Tstart::Real = 0.0
    Tend::Real = 44700
    theta::Float64= 0.5
    T_save::Real = 50

    ##''''''''''''''Physical parameter''''''''''''''##
    η::Float64 = 7.29e-5                             #angular speed of Earth rotation        (s^(-1))
    latitude::Float64 = 52.0
    f::Float64 = 2*η*sin(latitude*(π/180))           #coriolis parameter                     (s^(-1))
    g::Float64 = 9.81                                #Gravitational constant                 (ms^(-2))
    U_start::Float64 = 0.5                           #Background current amplitude           (ms^(-1))
    σ::Float64 = 2*pi/44700                          #Tidal frequency                        (s^(-1))
    cD::Float64 = 0.0025                             #Drag coefficient                       ()

    ##''''''''''''''Stabilization Parameters''''''''''''''##
    α::Float64 = 1e-6                                #Based on ν
    ν::Float64 = 1e-6                                #From Anna Louka
end


