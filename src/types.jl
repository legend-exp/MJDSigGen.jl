# This file is a part of MJDSigGen, licensed under the MIT License (MIT).

struct CartPoint{T<:Real}
    x::T
    y::T
    z::T
end

CartPoint{T}() where {T} = CartPoint{T}(0, 0, 0)


struct CylPoint{T<:Real}
    r::T
    phi::T
    z::T
end

CylPoint{T}() where {T} = CylPoint{T}(0, 0, 0)


mutable struct VelocityLookup
    e::Cfloat
    e100::Cfloat
    e110::Cfloat
    e111::Cfloat
    h100::Cfloat
    h110::Cfloat
    h111::Cfloat
    ea::Cfloat
    eb::Cfloat
    ec::Cfloat
    ebp::Cfloat
    ecp::Cfloat
    ha::Cfloat
    hb::Cfloat
    hc::Cfloat
    hbp::Cfloat
    hcp::Cfloat
    hcorr::Cfloat
    ecorr::Cfloat

    VelocityLookup() = new()
end


mutable struct SigGenSetup
    # general
    verbosity::Cint

    # geometry
    xtal_length::Cfloat
    xtal_radius::Cfloat
    pc_length::Cfloat
    pc_radius::Cfloat
    wrap_around_radius::Cfloat
    ditch_depth::Cfloat
    ditch_thickness::Cfloat
    bottom_taper_length::Cfloat
    hole_length::Cfloat
    hole_radius::Cfloat
    outer_taper_length::Cfloat
    taper_angle::Cfloat
    inner_taper_length::Cfloat
    outer_taper_width::Cfloat
    inner_taper_width::Cfloat
    top_bullet_radius::Cfloat
    bottom_bullet_radius::Cfloat
    hole_bullet_radius::Cfloat
    Li_thickness::Cfloat
    vacuum_gap::Cfloat

    # electric fields & weighing potentials
    xtal_grid::Cfloat
    impurity_z0::Cfloat
    impurity_gradient::Cfloat
    impurity_quadratic::Cfloat
    impurity_surface::Cfloat
    impurity_radial_add::Cfloat
    impurity_radial_mult::Cfloat
    impurity_rpower::Cfloat
    xtal_HV::Cfloat
    capacitance::Cfloat
    max_iterations::Cint
    write_field::Cint
    write_WP::Cint
    bulletize_PC::Cint

    # file names
    drift_name::NTuple{256,Cchar}
    field_name::NTuple{256,Cchar}
    wp_name::NTuple{256,Cchar}
    
    # signal calculation 
    xtal_temp::Cfloat
    preamp_tau::Cfloat
    time_steps_calc::Cint
    step_time_calc::Cfloat
    step_time_out::Cfloat
    charge_cloud_size::Cfloat
    use_acceleration::Cint
    use_repulsion::Cint
    use_diffusion::Cint
    energy::Cfloat
    
    coord_type::Cint
    ntsteps_out::Cint
    
    # data for fields.c
    rmin::Cfloat
    rmax::Cfloat
    rstep::Cfloat
    zmin::Cfloat
    zmax::Cfloat
    zstep::Cfloat
    rlen::Cint
    zlen::Cint
    v_lookup_len::Cint
    v_lookup::Ptr{VelocityLookup}

    # for fieldgen:
    v::NTuple{2,Ptr{Ptr{Cdouble}}}
    eps::Ptr{Ptr{Cdouble}}
    eps_dr::Ptr{Ptr{Cdouble}}
    eps_dz::Ptr{Ptr{Cdouble}}
    impurity::Ptr{Ptr{Cdouble}}
    vfraction::Ptr{Ptr{Cdouble}}
    s1::Ptr{Ptr{Cdouble}}
    s2::Ptr{Ptr{Cdouble}}
    vsave::Ptr{Ptr{Cdouble}}
    point_type::Ptr{Ptr{Cchar}}
    undepleted::Ptr{Ptr{Cchar}}
    fully_depleted::Cint
    bubble_volts::Cfloat
    Emin::Cfloat
    rho_z_spe::NTuple{1024,Cfloat}
    dr::NTuple{2,Ptr{Ptr{Cdouble}}}
    dz::NTuple{2,Ptr{Ptr{Cdouble}}}


    # for siggen
    efld::Ptr{Ptr{CylPoint{Cfloat}}}
    wpot::Ptr{Ptr{Cfloat}}
    
    config_file_name::NTuple{256,Cchar}
    
    # data for calc_signal.c
    dpath_e::Ptr{CartPoint{Cfloat}}
    dpath_h::Ptr{CartPoint{Cfloat}}
    surface_drift_vel_factor::Cfloat
    instant_vel_e::Ptr{CartPoint{Cfloat}}
    instant_vel_h::Ptr{CartPoint{Cfloat}}
    instant_charge_size_e::Ptr{Cfloat}
    instant_charge_size_h::Ptr{Cfloat}
    initial_vel::Cfloat
    final_vel::Cfloat
    dv_dE::Cfloat
    v_over_E::Cfloat
    final_charge_size::Cdouble

    function SigGenSetup()
        s = new()
        fill!(unsafe_wrap(Array, Ptr{UInt8}(pointer_from_objref(s)), sizeof(s)), 0)
        return s
    end
end

SigGenSetup(config_filename::AbstractString) = signal_calc_init(config_filename)

function tup2str(tup::NTuple{N,C}) where {N, C<:Union{AbstractChar, Cchar}}
    out = Vector{UInt8}(undef, N)

    for i in 1:N
        c = @inbounds tup[i]
        c == 0 && (resize!(out, i-1); break)
        @inbounds out[i] = tup[i]
    end

    return String(out)
end

function Base.getproperty(setup::SigGenSetup, sym::Symbol)
    if sym in (:config_file_name, :drift_name, :field_name, :wp_name)
        return tup2str(getfield(setup, sym))
    elseif sym in (:dpath_e, :dpath_h, :instant_vel_e, :instant_vel_h)
        ptr = getfield(setup, sym)
        N = setup.time_steps_calc
        v = unsafe_wrap(Vector{CartPoint{Cfloat}}, ptr, N)
        return transpose(reshape(reinterpret(Cfloat, v), (3, N)))
    elseif sym in (:instant_charge_size_e, :instant_charge_size_h)
        ptr = getfield(setup, sym)
        return unsafe_wrap(Vector{Cfloat}, ptr, setup.time_steps_calc)
    else
        return getfield(setup, sym)
    end
end
