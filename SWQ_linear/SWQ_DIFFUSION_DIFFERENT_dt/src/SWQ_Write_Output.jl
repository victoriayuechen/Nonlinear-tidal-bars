using WriteVTK
using Gridap


##''''''''''''''Save function''''''''''''''##
function writing_output(dir, name, x, Ω, Tend, P, Param, spinup_save, h, save_CSV)
    @unpack T_save, B, L, nx_start, ny_start, nx, ny, CSVname = Param
    if save_CSV
        x_hep = []
        y_hep = []
        i_grid = nx_start
        j_grid = ny_start
        while i_grid <= L
            append!(x_hep, i_grid)
            i_grid += B/nx
        end
        while j_grid <= B
            append!(y_hep, j_grid)
            j_grid += L/ny
        end
        probe = [Point(i, j) for i in x_hep, j in y_hep]
        
        lDa = zeros(Float64, 1, length(probe))
        prbDa = DataFrame(lDa, :auto)
    end

    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,name))do pvd
            for (x,t) in x
                if t==Tend && spinup_save
                    global spinup_solution = x
                end
                if t%(T_save) ==0
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ])
                    println("Done $t/$Tend")
                    if save_CSV
                        Di = interpolate_everywhere(ζ, P(0.0))
                        push!(prbDa, Di.(probe))
                    end
                end
                
            end
            if save_CSV
                prbDa = Matrix(prbDa)
                writedlm(CSVname, prbDa, ',')
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,name)) do pvd
            for (x,t) in x
                if t==Tend && spinup_save
                    global spinup_solution = x 
                end
                
                if t%(T_save) ==0
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ])
                    println("Done $t/$Tend")
                    if save_CSV
                        Di = interpolate_everywhere(ζ, P(0.0))
                        push!(prbDa, Di.(probe))
                    end
                end
            end
            if save_CSV
                prbDa = Matrix(prbDa)
                writedlm(CSVname, prbDa, ',')
            end
        end
    end
end