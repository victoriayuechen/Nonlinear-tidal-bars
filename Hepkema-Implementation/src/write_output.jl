using HDF5, Dates, Printf, Parameters, DataFrames, CSV

includet("read_output.jl")

"""
looks at last logbook entry and returns the next
"""
function next_output_folder(logbookFile)
    df = DataFrame(CSV.File(logbookFile,header=1,delim=" ; "))
    
    time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")   # get current time
    date       = Dates.format(now(), "yyyy_mm_dd")      # get current date
    lb_enteries_today = filter(x->occursin(date,x), df.output_folder)     # find logbook enteries from today
    if isempty(lb_enteries_today)
        last_nr = 0
    else
        last_nr = parse(Int32, lb_enteries_today[end][end-1:end])   # get the last number
    end
    new_nr = @sprintf("%02.0f",last_nr+1)  # add new number and make string
    output_folder = date * "__" * new_nr   # make new folder name like "2020_08_25__03"
    return output_folder
end


"""
write: time ; output_folder ; description to logbook.csv
return output_folder
"""
function write_to_logbook(output_folder,logbookFile,description_logbook)
    
    time = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")   # get current time
    new_line = DataFrame(time = [time], folder = [output_folder], description = [description_logbook])
    CSV.write(logbookFile, new_line, append = true, delim=" ; ")
end

"""
make folder to store output and copy the current state of the code
returns output_path (to put in Params struct)
"""
function make_output_dir(;
    test=false, 
    out_path = "./output/",
    scr_path="./", 
    description_long = "no description", 
    description_logbook=" ",
    logbook = false,
    logbookFile = "logbook.csv")
    
    #make out_path if it does not exist
    if !isdir(out_path)
        mkpath(out_path)
    end

    if test
        println("test not in logbook")
        output_folder = "test"
    else
        output_folder = next_output_folder(logbookFile)
    end

    #write run in logbook.csv
    if logbook
        write_to_logbook(output_folder,logbookFile,description_logbook)
    end

    #make output folder with date (or test)
    output_path_date = mkpath(out_path * output_folder) * "/"

    #make file with description of the run
    open(output_path_date * "description.txt", "w") do io
        write(io, description_long)
        write(io,"\n")
    end
    
    #copy current state of program to output_path/src
    scr_path_cp = mkpath(output_path_date * "scr") * "/"
    files = filter(x->occursin(".jl",x), readdir(scr_path))
    files = [files; ["Project.toml"]]  
    for file in files
        cp(scr_path * file, scr_path_cp * file, force=true)
    end

    #make folder to store the output arrays and return it
    array_path = mkpath(output_path_date * "arrays") * "/"
    return array_path
end

"""
save grid to hfdf5 file named grid.hfd5
"""
function output_grid(Pars,Grid)
    @unpack nx,ny,Δx,Δy,output_path = Pars
    h5open(output_path * "grid.h5","cw") do fid
        fid["x"]  = Grid.x
        fid["y"]  = Grid.y
        fid["xu"] = Grid.xᵤ
        fid["xv"] = Grid.xᵥ
        fid["yu"] = Grid.yᵤ
        fid["yv"] = Grid.yᵥ
    end
end



"""
save ζ,u,v,c to hfd5 file named swec.hdf5
its organised as: /t/var
where t is a string with the current time and var ∈ {zeta,u,v,c}
"""
function output_swec(t,Ψ,A,Pars)
    @unpack output_path = Pars
    h5open(output_path * "swec.h5","cw") do fid
        # create new group named after the current simulation time
        fid_t = create_group(fid,string(t))

        # fill with datasets
        fid_t["zeta"] = Ψ[:,:,1]
        fid_t["u"]    = Ψ[:,:,2]
        fid_t["v"]    = Ψ[:,:,3]
        fid_t["c"]    = Ψ[:,:,4]
        fid_t["div_qs"] = A.div_q⃗s
        fid_t["div_qb"] = A.div_q⃗b
    end
end

"""
save residual ζ,u,v,c to hfd5 file named residual_swec.h5
all saved residuals are on C grid
"""
function output_residual_swec(sol,savedval,Pars)
    @unpack nx,ny,output_path = Pars

    ζ = @view sol[:,:,1,:] # x,y,var,t
    u = @view sol[:,:,2,:]
    v = @view sol[:,:,3,:]
    c = @view sol[:,:,4,:]

    
    
    ζ_avg  = reshape(mean(ζ,dims=3),(nx,ny))
    u_avgᵤ = reshape(mean(u,dims=3),(nx,ny))
    v_avgᵥ = reshape(mean(v,dims=3),(nx,ny))
    c_avg  = reshape(mean(c,dims=3),(nx,ny))

    u_avg = zeros(nx,ny)
    v_avg = zeros(nx,ny)
    interpolate_UC!(u_avg,u_avgᵤ)
    interpolate_VC!(v_avg,v_avgᵥ)
    
    qu_avgᵤ = zeros(nx,ny)
    qv_avgᵥ = zeros(nx,ny)
    for tn in 1:length(savedval)
        qu_avgᵤ += savedval[tn][1]
        qv_avgᵥ += savedval[tn][2]
    end
    qu_avgᵤ /= length(savedval)
    qv_avgᵥ /= length(savedval)
    qu_avg = zeros(nx,ny)
    qv_avg = zeros(nx,ny)
    interpolate_UC!(qu_avg,qu_avgᵤ)
    interpolate_VC!(qv_avg,qv_avgᵥ)

    h5open(output_path * "residual_swec.h5","cw") do fid
        # fill with datasets
        fid["residual_zeta"]  = ζ_avg
        fid["residual_u"]     = u_avg
        fid["residual_v"]     = v_avg
        fid["residual_c"]     = c_avg
        fid["residual_qu"]    = qu_avg
        fid["residual_qv"]    = qv_avg
    end


end


"""
save h to hfd5 file named bottomheight.hdf5
its organised as: /t/h
where t is a string with the current time in years
"""
function output_bottomheight(τ,A,Pars)
    @unpack output_path,ε,sec_in_year = Pars
    h5open(output_path * "bottomheight.h5","cw") do fid
        # create new group named after the current simulation time
        fid_t = create_group(fid,string(τ/(ε*sec_in_year)))

        # fill with datasets
        fid_t["h"]          = A.h
    end
end

