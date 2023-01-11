using HDF5, Parameters

"""
read grid file
returns x,y,xᵤ,xᵥ,yᵤ,yᵥ
"""
function read_grid(path)
    gridfile = h5open(path * "grid.h5","r")
    x  = read(gridfile,"x");
    y  = read(gridfile,"y");
    xᵤ = read(gridfile,"xu");
    xᵥ = read(gridfile,"xv");
    yᵤ = read(gridfile,"yu");
    yᵥ = read(gridfile,"yv");
    close(gridfile)
    return x,y,xᵤ,xᵥ,yᵤ,yᵥ
end


"""
read output swec.h5
returns t, t, div_qs, div_qb
The variables are 3d arrays (x,y,t)
t is in hours
"""
function read_divq(path ; physical_domain=true)
    # read grid
    x,y,xᵤ,xᵥ,yᵤ,yᵥ = read_grid(path);
    nx = length(x)
    ny = length(y) + 1

    # read swec
    dset = h5open(path * "swec.h5","r")
    str_t = keys(dset) # list of timesteps (as strings)
    Nt = length(str_t)  # number of timesteps
    
    

    if physical_domain
        # leave h useless ny grid point at C and U grid
        div_qs = zeros(nx,ny-1,Nt);
        div_qb = zeros(nx,ny-1,Nt);
        for n in eachindex(str_t)
            div_qs[:,:,n] .= read(dset, str_t[n] * "/div_qs")[:,1:end-1]
            div_qb[:,:,n] .= read(dset, str_t[n] * "/div_qb")[:,1:end-1]
        end
    else
        div_qs = zeros(nx,ny,Nt);
        div_qb = zeros(nx,ny,Nt);
        
        for n in eachindex(str_t)
            div_qs[:,:,n] .= read(dset, str_t[n] * "/div_qs")
            div_qb[:,:,n] .= read(dset, str_t[n] * "/div_qb")
            
        end
    end


    close(dset)

    # turn strings of str_t in floats
    t = map(x->parse(Float64,x),str_t)
    # get pertubation to sort t
    p = sortperm(t)

    # sort t and the rest.
    t = t[p] ./ 3600;
    div_qs .= div_qs[:,:,p];
    div_qb .= div_qb[:,:,p];
    

    return t, div_qs, div_qb
end

"""
read output swec.h5
returns t, ζ, u, v, c. 
The variables are 3d arrays (x,y,t)
t is in hours
"""
function read_swec(path ; physical_domain=true)
    # read grid
    x,y,xᵤ,xᵥ,yᵤ,yᵥ = read_grid(path);
    nx = length(x)
    ny = length(y) + 1

    # read swec
    dset = h5open(path * "swec.h5","r")
    str_t = keys(dset) # list of timesteps (as strings)
    Nt = length(str_t)  # number of timesteps
    
    

    if physical_domain
        # leave h useless ny grid point at C and U grid
        ζ = zeros(nx,ny-1,Nt);
        u = zeros(nx,ny-1,Nt);
        v = zeros(nx,ny  ,Nt);
        c = zeros(nx,ny-1,Nt);
        for n in eachindex(str_t)
            ζ[:,:,n] .= read(dset, str_t[n] * "/zeta")[:,1:end-1]
            u[:,:,n] .= read(dset, str_t[n] * "/u")[:,1:end-1]
            v[:,:,n] .= read(dset, str_t[n] * "/v")
            c[:,:,n] .= read(dset, str_t[n] * "/c")[:,1:end-1]
        end
    else
        ζ = zeros(nx,ny,Nt);
        u = zeros(nx,ny,Nt);
        v = zeros(nx,ny,Nt);
        c = zeros(nx,ny,Nt);
        for n in eachindex(str_t)
            ζ[:,:,n] .= read(dset, str_t[n] * "/zeta")
            u[:,:,n] .= read(dset, str_t[n] * "/u")
            v[:,:,n] .= read(dset, str_t[n] * "/v")
            c[:,:,n] .= read(dset, str_t[n] * "/c")
        end
    end


    close(dset)

    # turn strings of str_t in floats
    t = map(x->parse(Float64,x),str_t)
    # get pertubation to sort t
    p = sortperm(t)

    # sort t and the rest.
    t = t[p] ./ 3600;
    ζ .= ζ[:,:,p];
    u .= u[:,:,p];
    v .= v[:,:,p];
    c .= c[:,:,p];

    return t, ζ, u, v, c
end

"""
read residual_swec files
returns ζ_avg, u_avg, v_avg, c_avg on C grid
"""
function read_residual_swec(path)
    # read residual file
    residualfile = h5open(path * "residual_swec.h5","r")
    ζ_avg = read(residualfile, "residual_zeta")[:,1:end-1]
    u_avg = read(residualfile, "residual_u")[:,1:end-1]
    v_avg = read(residualfile, "residual_v")[:,1:end-1]
    c_avg = read(residualfile, "residual_c")[:,1:end-1]
    qu_avg = read(residualfile, "residual_qu")[:,1:end-1]
    qv_avg = read(residualfile, "residual_qv")[:,1:end-1]
    close(residualfile)
    return ζ_avg, u_avg, v_avg, c_avg, qu_avg, qv_avg
end

"""
read hput bottomheight.h5
returns t, h. 
h is a 3d arrays (x,y,t)
t is in years
"""
function read_bottomheight(path; physical_domain=true)    
    # read swec
    dset = h5open(path * "bottomheight.h5","r")
    str_t = keys(dset) # list of timesteps (as strings)
    Nt = length(str_t)  # number of timesteps
    nx,ny = size(read(dset, str_t[1] * "/h"))

    if physical_domain
        # leave out useless ny grid point at C and U grid
        h = zeros(nx,ny-1,Nt);
        for n in eachindex(str_t)
            h[:,:,n] .= read(dset, str_t[n] * "/h")[:,1:end-1]
        end
    else
        h = zeros(nx,ny,Nt);
        for n in eachindex(str_t)
            h[:,:,n] .= read(dset, str_t[n] * "/h")
        end
    end
    close(dset)

    # turn strings of str_t in floats
    t = map(x->parse(Float64,x),str_t)
    # get pertubation that sorts τ
    p = sortperm(t)

    # sort τ, calculate t = ετ and get h.
    t = t[p] # in years
    h .= h[:,:,p];

    return t, h
end

"""
sync hput of gemini to local
returns path to directory with hput of the given version (name of hput folder)
"""
function retrieve_output_gemini(version;old=false)
    output_path_local = "/users/tjebbe/Documents/Morpho_output_julia/"
    if old
        path_gemini= "/nethome/3678636/morpho/output/" * version
    else
        path_gemini= "/nethome/3678636/morpho2/output/" * version
    end
    cmd = `rsync -avz uu:"$path_gemini" "$output_path_local"`
    run(cmd)
    return output_path_local * version * "/"
end


"""
calculate hrms vs t
"""
function calculate_hrms(h)
    nx,ny,nt = size(h)
    hrms = zeros(nt)
    for n in 1:nt
        for j in 1:ny-1
            for i in 2:nx-1
                hrms[n] += h[i,j,n]^2
            end
        end
        hrms[n] /= (nx-2)*(ny-1)
    end
    return sqrt.(hrms)
end

