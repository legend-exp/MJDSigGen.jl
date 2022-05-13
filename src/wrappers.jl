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


function get_signal!(signal::DenseArray{Float32, 1}, setup::SigGenSetup, location::NTuple{3, Real})
    sl = length(signal)
    nsteps = setup.ntsteps_out
    sl != nsteps && throw(
        ArgumentError("signal length ($sl) must equal number of time steps ($nsteps)"))

    pt = CartPoint{Cfloat}(location[1], location[2], location[3])

    result = ccall(
        @sgsym(:get_signal), Cint,
        (CartPoint{Cfloat}, Ptr{Float32}, Ptr{SigGenSetup}),
        pt, signal, Ref(setup)
    )
    result < 0 && throw(ArgumentError("point not in crystal or has no field: $pt"))

    return signal
end


get_signal!(setup::SigGenSetup, location::NTuple{3, Real}) =
    get_signal!(zeros(Float32, setup.ntsteps_out), setup, location)


function outside_detector(setup, location::NTuple{3, Real})
    pt = CartPoint{Cfloat}(location[1], location[2], location[3])

    r = ccall(
        @sgsym(:outside_detector), Cint,
        (CartPoint{Cfloat}, Ptr{SigGenSetup}),
        pt, Ref(setup)
    )

    return r != 0
end


"""
    nearest_field_grid_index(setup::SigGenSetup, location::NTuple{3, Real})

Returns:
* (:outside, i, j), if outside crystal or too far from a valid grid point
* (:interpol, i, j, if interpolation is okay
* (:extrapol, i, j), if we can find a point but extrapolation is needed
"""
function nearest_field_grid_index(setup::SigGenSetup, location::NTuple{3, Real})
    r, ϕ, z = cart2cyl(location ...)
    cyl_pos_val = CylPoint{Cfloat}(r, ϕ, z)
    cyl_idx_ref = Ref(CylPoint{Cint}(0, 0, 0))

    retcode = ccall(
        @sgsym(:nearest_field_grid_index), Cint,
        (CylPoint{Cfloat}, Ptr{CylPoint{Cint}}, Ptr{SigGenSetup}, ),
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
        throw(ArgumentError(
            "Don't know how to interpret return code $retcode of
            nearest_field_grid_index unknown"
        ))
    end
end


const fieldgen_exe = joinpath(dirname(@__FILE__), "..", "deps", "usr", "bin", "mjd_fieldgen")

function fieldgen(config_filename::AbstractString)
    setup = read_config(config_filename)

    mkpath(dirname(joinpath(dirname(config_filename), setup.field_name)))
    mkpath(dirname(joinpath(dirname(config_filename), setup.wp_name)))

    run(`$fieldgen_exe $config_filename`)
end


function get_drift_velocity(setup::SigGenSetup, location::NTuple{3, Real}, t::Symbol)
	if t==:e
		q = -1;
	elseif t==:h
		q = +1;
	else
		error("Charge carrier type must be :e or :h")
	end

    pt	= CartPoint{Cfloat}(location[1], location[2], location[3])
	vel	= Ref(CartPoint{Cfloat}(0, 0, 0))

    ccall(
        @sgsym(:drift_velocity), Cint,
        (CartPoint{Cfloat}, Float32, Ptr{CartPoint{Cfloat}}, Ptr{SigGenSetup}),
        pt, q, vel, Ref(setup)
    ) < 0 && error("Point not in crystal or has no field: $pt")

	#asd = unsafe_load(vel);
    return vel.x;
end

function get_drift_velocity_w_Eadd(setup::SigGenSetup, location::NTuple{3, Real}, t::Symbol, Eadd_cart::NTuple{3, Real})
	if t==:e
		q = -1;
	elseif t==:h
		q = +1;
	else
		error("Charge carrier type must be :e or :h")
	end
	
    pt	= CartPoint{Cfloat}(location[1], location[2], location[3])
	vel	= Ref(CartPoint{Cfloat}(0, 0, 0))
	Er, Eϕ, Ez = cart2cyl(Eadd_cart ...);
	if(Eadd_cart[1]<0)
		Er =-Er;
		Eϕ = 0;
	end	
	Eadd = CylPoint{Cfloat}(Er, Eϕ, Ez);

    ccall(
        @sgsym(:drift_velocity_w_Eadd), Cint,
        (CartPoint{Cfloat}, Float32, Ptr{CartPoint{Cfloat}},
        Ptr{SigGenSetup}, CylPoint{Cfloat}),
        pt, q, vel, Ref(setup), Eadd
    ) < 0 && error("Point not in crystal or has no field: $pt")

	#asd = unsafe_load(vel);
    return vel.x;
end


function get_drift_velocity_from_Efield(setup::SigGenSetup, location::NTuple{3, Real}, t::Symbol, Efield_cart::NTuple{3, Real})
	if t==:e
		q = -1;
	elseif t==:h
		q = +1;
	else
		error("Charge carrier type must be :e or :h")
	end
	
    pt	= CartPoint{Cfloat}(location[1], location[2], location[3])
	vel	= Ref(CartPoint{Cfloat}(0, 0, 0))
	Er, Eϕ, Ez = cart2cyl(Efield_cart...);
	if(Efield_cart[1]<0)
		Er =-Er;
		Eϕ = 0;
	end	
	Efield = CylPoint{Cfloat}(Er, Eϕ, Ez);

    ccall(
        @sgsym(:drift_velocity_from_Efield), Cint,
        (CartPoint{Cfloat}, Float32, Ptr{CartPoint{Cfloat}}, Ptr{SigGenSetup}, CylPoint{Cfloat}),
        pt, q, vel, Ref(setup), Efield
    ) < 0 && error("Point not in crystal or has no field: $pt")

    return Tuple(vel.x);
end
