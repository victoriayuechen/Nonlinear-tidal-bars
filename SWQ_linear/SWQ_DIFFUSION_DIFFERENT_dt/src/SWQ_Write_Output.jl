using WriteVTK
using Gridap
using DataFrames
using CSV
using DelimitedFiles

##''''''''''''''Save the beauty''''''''''''''##
function writing_output(dir, x, Ω, Tend, P)
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
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            for (x,t) in x
                u,ζ = x
                if t%(400) ==0
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                    println("Saved")
                    println("done $t/$Tend")
                end
                Di = interpolate_everywhere(ζ, P(0.0))
                push!(prbDa, Di.(probe)) 
            end
            prbDa = Matrix(prbDa)
            writedlm("zeta_diff.csv", prbDa, ',')
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output")) do pvd
            for (x,t) in x
                u,ζ = x
                if t%(400) ==0
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                    println("Saved")
                    println("done $t/$Tend")
                end
                Di = interpolate_everywhere(ζ, P(0.0))
                push!(prbDa, Di.(probe)) 

            end
            prbDa = Matrix(prbDa)
            writedlm("zeta_diff.csv", prbDa, ',')
        end
    end
end
