using WriteVTK
using Gridap
using CSV
using DataFrames
using DelimitedFiles

##''''''''''''''Save the beauty''''''''''''''##
function writing_output(dir, name, x, Ω, Tend, P, namecsv)
    x_hep = []
    y_hep = []
    i_grid = 0
    j_grid = 25
    while i_grid <= 10000
        append!(x_hep, i_grid)
        i_grid += 100
    end
    while j_grid <= 1000
        append!(y_hep, j_grid)
        j_grid += 50
    end
    probe = [Point(i, j) for i in x_hep, j in y_hep]
    
    lDa = zeros(Float64, 1, length(probe))
    prbDa = DataFrame(lDa, :auto)

    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,name))do pvd
            for (x,t) in x
                u,D = x
                # Di = interpolate_everywhere((D-3+h), P(0.0))
                # push!(prbDa, Di.(probe)) 
                if abs(t - Tend) < 1e-5
                    println("U see it?")
                    global xn = x 
                    println("Then fucking save it")
                end

                if t%(100) < 1e-5
                    u,D = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>(D+h-3.0),"h"=>h])
                    println("Saved")
                    println("done $t/$Tend")
                    Di = interpolate_everywhere((D-3+h), P(0.0))
                    push!(prbDa, Di.(probe))
                end
                
            end
            prbDa = Matrix(prbDa)
            writedlm(namecsv, prbDa, ',')
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,name)) do pvd
            for (x,t) in x
                # u,D = x
                # Di = interpolate_everywhere((D-3+h), P(0.0))
                # push!(prbDa, Di.(probe)) 
                if abs(t - Tend) < 1e-5
                    global xn = x 
                    println("Saved in xn")
                    println("Spin up done")
                end
                if t%(100) < 1e-5
                    
                    u,D = x
                    Di = interpolate_everywhere((D-3+h), P(0.0))
                    push!(prbDa, Di.(probe))
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>(D+h-3.0),"h"=>h])
                    # println("Saved")
                    # println("done $t/$Tend")
                end 
            end
            prbDa = Matrix(prbDa)
            writedlm(namecsv, prbDa, ',')
        end
    end
end