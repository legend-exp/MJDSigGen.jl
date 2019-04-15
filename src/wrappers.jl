# This file is a part of MJDSigGen, licensed under the MIT License (MIT).

function read_config!(setup::Struct_MJD_Siggen_Setup, config_filename::AbstractString)
    ccall(
        @sgsym(:read_config), Cint,
        (Cstring, Ptr{Struct_MJD_Siggen_Setup}),
        config_filename, Ref(setup)
    ) != 0 && error("read_config failed.")
    setup
end

read_config(config_filename::AbstractString) =
    read_config!(Struct_MJD_Siggen_Setup(), config_filename)


function signal_calc_init!(setup::Struct_MJD_Siggen_Setup, config_filename::AbstractString)
    ccall(
        @sgsym(:signal_calc_init), Cint,
        (Cstring, Ptr{Struct_MJD_Siggen_Setup}),
        config_filename, Ref(setup)
    ) != 0 && error("signal_calc_init failed.")
    finalizer(signal_calc_finalize!, setup)
    setup
end

signal_calc_init(config_filename::AbstractString) =
    signal_calc_init!(Struct_MJD_Siggen_Setup(), config_filename)


function signal_calc_finalize!(setup::Struct_MJD_Siggen_Setup)
    ccall(
        @sgsym(:signal_calc_finalize), Cint,
        (Ptr{Struct_MJD_Siggen_Setup},),
        Ref(setup)
    ) != 0 && error("signal_calc_finalize failed.")
    setup
end


function get_signal!(signal::DenseArray{Float32, 1}, setup::Struct_MJD_Siggen_Setup, location::NTuple{3})
    (length(eachindex(signal)) < setup.ntsteps_out) && throw(BoundsError())

    pt = Struct_point(location[1], location[2], location[3])

    ccall(
        @sgsym(:get_signal), Cint,
        (Struct_point, Ptr{Float32}, Ptr{Struct_MJD_Siggen_Setup}),
        pt, signal, Ref(setup)
    ) < 0 && error("Point not in crystal or has no field: $pt")

    signal
end


get_signal!(setup::Struct_MJD_Siggen_Setup, location::NTuple{3}) =
    get_signal!(zeros(Float32, setup.ntsteps_out), setup, location)


function _drift_path_ptr(setup::Struct_MJD_Siggen_Setup, t::Symbol)
    if t == :e
        setup.dpath_e
    elseif t == :h
        setup.dpath_h
    else
        error("Charge carrier type must be :e or :h")
    end
end


function drift_path_len(setup::Struct_MJD_Siggen_Setup, t::Symbol)
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


function drift_path!(path::DenseArray{Float32, 2}, setup::Struct_MJD_Siggen_Setup, t::Symbol)
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

drift_path(setup::Struct_MJD_Siggen_Setup, t::Symbol) =
    drift_path!(zeros(Float32, drift_path_len(setup, t), 3), setup, t)


function _instant_vel_ptr(setup::Struct_MJD_Siggen_Setup, t::Symbol)
    if t == :e
		setup.instant_vel_h
    elseif t == :h
        setup.instant_vel_h
    else
        error("Charge carrier type must be :e or :h")
    end
end


function instant_vel_len(setup::Struct_MJD_Siggen_Setup, t::Symbol)
    vel_ptr = _instant_vel_ptr(setup, t)
    n = setup.time_steps_calc
    @inbounds for i in 1:n
        pt = unsafe_load(path_ptr, i)
        if pt == Struct_point()
            return i - 1
        end
    end
    return n
end


function instant_vel!(vel::DenseArray{Float32, 2}, setup::Struct_MJD_Siggen_Setup, t::Symbol)
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

instant_vel(setup::Struct_MJD_Siggen_Setup, t::Symbol) =
    instant_vel!(zeros(Float32, instant_vel_len(setup, t), 3), setup, t)

function outside_detector(setup::Struct_MJD_Siggen_Setup, location::NTuple{3})
    pt = Struct_point(location[1], location[2], location[3])

    r = ccall(
        @sgsym(:outside_detector), Cint,
        (Struct_point, Ptr{Struct_MJD_Siggen_Setup}),
        pt, Ref(setup)
    )

    r != 0
end


"""
    nearest_field_grid_index(setup::Struct_MJD_Siggen_Setup, location::NTuple{3})

Returns:
* (:outside, i, j), if outside crystal or too far from a valid grid point
* (:interpol, i, j, if interpolation is okay
* (:extrapol, i, j), if we can find a point but extrapolation is needed
"""
function nearest_field_grid_index(setup::Struct_MJD_Siggen_Setup, location::NTuple{3})
    r, ϕ, z = cart2cyl(location ...)
    cyl_pos_val = Struct_cyl_pt(r, ϕ, z)
    cyl_idx_ref = Ref(Struct_cyl_int_pt(0, 0, 0))

    retcode = ccall(
        @sgsym(:nearest_field_grid_index), Cint,
        (Struct_cyl_pt, Ptr{Struct_cyl_int_pt}, Ptr{Struct_MJD_Siggen_Setup}, ),
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
    setup = Struct_MJD_Siggen_Setup()
    read_config!(setup, config_filename)

    mkpath(dirname(field_file_name(setup)))
    mkpath(dirname(wpot_file_name(setup)))

    run(`$fieldgen_exe -c $config_filename`)
end
