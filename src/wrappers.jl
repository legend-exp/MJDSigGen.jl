# This file is a part of MJDSigGen, licensed under the MIT License (MIT).

function read_config!(setup::SigGenSetup, config_filename::AbstractString)
    result = ccall(
        @sgsym(:read_config), Cint,
        (Cstring, Ptr{SigGenSetup}),
        config_filename, Ref(setup)
    )
    
    result != 0 && error("read_config failed.")
    return setup
end

read_config(config_filename::AbstractString) =
    read_config!(SigGenSetup(), config_filename)


function signal_calc_init!(setup::SigGenSetup, config_filename::AbstractString)
    result = ccall(
        @sgsym(:signal_calc_init), Cint,
        (Cstring, Ptr{SigGenSetup}),
        config_filename, Ref(setup)
    )
    result != 0 && error("signal_calc_init failed.")
    finalizer(signal_calc_finalize!, setup)
    return setup
end

signal_calc_init(config_filename::AbstractString) =
    signal_calc_init!(SigGenSetup(), config_filename)


function signal_calc_finalize!(setup::SigGenSetup)
    result = ccall(
        @sgsym(:signal_calc_finalize), Cint,
        (Ptr{SigGenSetup},),
        Ref(setup)
    )
    result != 0 && error("signal_calc_finalize failed.")
    return setup
end


function get_signal!(signal::DenseArray{Float32, 1}, setup::SigGenSetup, location::NTuple{3})
    sl = length(signal)
    nsteps = setup.ntsteps_out
    sl != nsteps && throw(
        BoundsError("signal length ($ls) must equal number of time steps ($nsteps)"))

    pt = Struct_point(location[1], location[2], location[3])

    result = ccall(
        @sgsym(:get_signal), Cint,
        (Struct_point, Ptr{Float32}, Ptr{SigGenSetup}),
        pt, signal, Ref(setup)
    )
    result < 0 && error("Point not in crystal or has no field: $pt")

    return signal
end


get_signal!(setup::SigGenSetup, location::NTuple{3}) =
    get_signal!(zeros(Float32, setup.ntsteps_out), setup, location)


function _drift_path_ptr(setup::SigGenSetup, t::Symbol)
    if t == :e
        setup.dpath_e
    elseif t == :h
        setup.dpath_h
    else
        error("Charge carrier type must be :e or :h")
    end
end


function drift_path_len(setup::SigGenSetup, t::Symbol)
    path_ptr = _drift_path_ptr(setup, t)
    n = setup.time_steps_calc
    @inbounds for i in 1:n
        pt = unsafe_load(path_ptr, i)
        if pt == Struct_point()
            return i - 1
        end
    end
    return n
end


function drift_path!(path::DenseArray{Float32, 2}, setup::SigGenSetup, t::Symbol)
    (size(path, 2) < 2) && throw(BoundsError())

    path_ptr = _drift_path_ptr(setup, t)
    n = min(size(path, 1), setup.time_steps_calc)
    path_idxs = axes(path, 1)
    @inbounds for i in 1:n
        pt = unsafe_load(path_ptr, i)
        j = path_idxs[i]
        path[j, 1] = pt.x
        path[j, 2] = pt.y
        path[j, 3] = pt.z
    end

    path
end

drift_path(setup, t::Symbol) =
    drift_path!(zeros(Float32, drift_path_len(setup, t), 3), setup, t)


function _instant_vel_ptr(setup::SigGenSetup, t::Symbol)
    if t == :e
		setup.instant_vel_e
    elseif t == :h
        setup.instant_vel_h
    else
        error("Charge carrier type must be :e or :h")
    end
end


function instant_vel_len(setup::SigGenSetup, t::Symbol)
    vel_ptr = _instant_vel_ptr(setup, t)
    n = setup.time_steps_calc
    @inbounds for i in 1:n
        pt = unsafe_load(vel_ptr, i)
        if pt == Struct_point()
            return i - 1
        end
    end
    return n
end


function instant_vel!(vel, setup, t::Symbol)
    (size(vel, 2) < 2) && throw(BoundsError())

    vel_ptr = _instant_vel_ptr(setup, t)
    n = min(size(vel, 1), setup.time_steps_calc)
    vel_idxs = axes(vel, 1)
    @inbounds for i in 1:n
        pt = unsafe_load(vel_ptr, i)
        j = vel_idxs[i]
        vel[j, 1] = pt.x
        vel[j, 2] = pt.y
        vel[j, 3] = pt.z
    end

	vel 
end

instant_vel(setup, t::Symbol) =
    instant_vel!(zeros(Float32, instant_vel_len(setup, t), 3), setup, t)


function _instant_charge_size_ptr(setup, t::Symbol)
    if t == :e
		setup.instant_charge_size_e
    elseif t == :h
        setup.instant_charge_size_h
    else
        error("Charge carrier type must be :e or :h")
    end
end


function instant_charge_size_len(setup, t::Symbol)
    size_ptr = _instant_charge_size_ptr(setup, t)
    n = setup.time_steps_calc
    @inbounds for i in 1:n
        pt = unsafe_load(size_ptr, i)
    end
    return n
end


function instant_charge_size!(cloudSize, setup, t::Symbol)
	(size(cloudSize, 1) < 2) && throw(BoundsError())

    size_ptr = _instant_charge_size_ptr(setup, t)
    n = min(size(cloudSize,1), setup.time_steps_calc)
    size_idxs = axes(cloudSize, 1)
    @inbounds for i in 1:n
        pt = unsafe_load(size_ptr, i)
        j = size_idxs[i]
        cloudSize[j] = pt
    end

	cloudSize
end

instant_charge_size(setup, t::Symbol) =
    instant_charge_size!(zeros(Float32, instant_charge_size_len(setup, t)), setup, t)

function outside_detector(setup, location::NTuple{3})
    pt = Struct_point(location[1], location[2], location[3])

    r = ccall(
        @sgsym(:outside_detector), Cint,
        (Struct_point, Ptr{SigGenSetup}),
        pt, Ref(setup)
    )

    return r != 0
end


"""
    nearest_field_grid_index(setup::SigGenSetup, location::NTuple{3})

Returns:
* (:outside, i, j), if outside crystal or too far from a valid grid point
* (:interpol, i, j, if interpolation is okay
* (:extrapol, i, j), if we can find a point but extrapolation is needed
"""
function nearest_field_grid_index(setup::SigGenSetup, location::NTuple{3})
    r, ϕ, z = cart2cyl(location ...)
    cyl_pos_val = Struct_cyl_pt(r, ϕ, z)
    cyl_idx_ref = Ref(Struct_cyl_int_pt(0, 0, 0))

    retcode = ccall(
        @sgsym(:nearest_field_grid_index), Cint,
        (Struct_cyl_pt, Ptr{Struct_cyl_int_pt}, Ptr{SigGenSetup}, ),
        cyl_pos_val, cyl_idx_ref, Ref(setup)
    )

    cyl_idx = cyl_idx_ref.x

    i,j = Int(cyl_idx.z + 1), Int(cyl_idx.r + 1)

    if retcode < 0
        (:outside, 0, 0)
    elseif retcode == 0
        (:interpol, i, j)
    elseif retcode == 1
        (:extrapol, i, j)
    else
        error("Don't know how to interpret return code $retcode of nearest_field_grid_index unknown")
    end
end


const fieldgen_exe = joinpath(dirname(@__FILE__), "..", "deps", "usr", "bin", "mjd_fieldgen")

function fieldgen(config_filename::AbstractString)
    setup = read_config(config_filename)

    mkpath(dirname(joinpath(dirname(config_filename), setup.field_name)))
    mkpath(dirname(joinpath(dirname(config_filename), setup.wp_name)))

    run(`$fieldgen_exe -c $config_filename`)
end


function get_drift_velocity(setup::SigGenSetup, location::NTuple{3}, t::Symbol)
	if t==:e
		q = -1;
	elseif t==:h
		q = +1;
	else
		error("Charge carrier type must be :e or :h")
	end

    pt	= Struct_point(location[1], location[2], location[3])
	vel	= Ref(Struct_point(0, 0, 0))

    ccall(
        @sgsym(:drift_velocity), Cint,
        (Struct_point, Float32, Ptr{Struct_point}, Ptr{SigGenSetup}),
        pt, q, vel, Ref(setup)
    ) < 0 && error("Point not in crystal or has no field: $pt")

	#asd = unsafe_load(vel);
    return vel.x;
	
	
end

function get_drift_velocity_w_Eadd(setup::SigGenSetup, location::NTuple{3}, t::Symbol, Eadd_cart::NTuple{3})
	if t==:e
		q = -1;
	elseif t==:h
		q = +1;
	else
		error("Charge carrier type must be :e or :h")
	end
	
    pt	= Struct_point(location[1], location[2], location[3])
	vel	= Ref(Struct_point(0, 0, 0))
	Er, Eϕ, Ez = cart2cyl(Eadd_cart ...);
	if(Eadd_cart[1]<0)
		Er =-Er;
		Eϕ = 0;
	end	
	Eadd = Struct_cyl_pt(Er, Eϕ, Ez);

    ccall(
        @sgsym(:drift_velocity_w_Eadd), Cint,
        (Struct_point, Float32, Ptr{Struct_point}, Ptr{SigGenSetup}, Struct_cyl_pt),
        pt, q, vel, Ref(setup), Eadd
    ) < 0 && error("Point not in crystal or has no field: $pt")

	#asd = unsafe_load(vel);
    return vel.x;
	
end


function get_drift_velocity_from_Efield(setup::SigGenSetup, location::NTuple{3}, t::Symbol, Efield_cart::NTuple{3})
	if t==:e
		q = -1;
	elseif t==:h
		q = +1;
	else
		error("Charge carrier type must be :e or :h")
	end
	
    pt	= Struct_point(location[1], location[2], location[3])
	vel	= Ref(Struct_point(0, 0, 0))
	Er, Eϕ, Ez = cart2cyl(Efield_cart ...);
	if(Efield_cart[1]<0)
		Er =-Er;
		Eϕ = 0;
	end	
	Efield = Struct_cyl_pt(Er, Eϕ, Ez);

    ccall(
        @sgsym(:drift_velocity_from_Efield), Cint,
        (Struct_point, Float32, Ptr{Struct_point}, Ptr{SigGenSetup}, Struct_cyl_pt),
        pt, q, vel, Ref(setup), Efield
    ) < 0 && error("Point not in crystal or has no field: $pt")

    return Tuple(vel.x);
	
end
