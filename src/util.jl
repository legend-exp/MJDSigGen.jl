# This file is a part of MJDSigGen, licensed under the MIT License (MIT).

function read_fields(setup::SigGenSetup)
    config_dir = dirname(setup.config_name)
    field_data = readdlm(joinpath(config_dir, setup.field_name), comments=true)
    wpot_data = readdlm(joinpath(config_dir, setup.wp_name), comments=true)

    Nr::Int = setup.xtal_radius / setup.xtal_grid + 1
    Nz::Int = setup.xtal_length / setup.xtal_grid + 1
    N = Nz * Nr

    @assert size(field_data, 1) == size(wpot_data, 1) == N

    E_pot = copy(reshape(view(field_data, :, 3), (Nz, Nr)))
    W_pot = copy(reshape(view(wpot_data, :, 3), (Nz, Nr)))

    E_abs = copy(reshape(view(field_data, :, 4), (Nz, Nr)))
    E_r = copy(reshape(view(field_data, :, 5), (Nz, Nr)))
    E_z = copy(reshape(view(field_data, :, 6), (Nz, Nr)))

    return E_pot, W_pot, E_abs, E_r, E_z
end


cyl2cart(r, ϕ, z) = (r * cos(ϕ), r * sin(ϕ), z)


function cart2cyl(x, y, z)
    r = √(x^2 + y^2)
    if x == 0
        ϕ = (y > 0) ? π/2 : -π/2
    else
        ϕ = (x < 0) ? atan(y/x) + π : atan(y/x)
    end
    return r, ϕ, z
end


function coll_effects_off!(setup::SigGenSetup)
    setup.energy            = 0
	setup.charge_cloud_size = 0
	setup.use_diffusion     = 0
	setup.use_acceleration  = 0
	setup.use_repulsion     = 0
	return setup
end

function with_coll_effects!(f::Function, setup::SigGenSetup, E::Real, ch_cld_size::Real)
	try
		setup.energy            = E
		setup.charge_cloud_size = ch_cld_size
		setup.use_diffusion     = 1
		setup.use_acceleration  = 1
		setup.use_repulsion     = 1

		return f()
	finally
		coll_effects_off!(setup)
	end
end


getδτ(charge_size::Real, speed::Real) = charge_size / speed

getδτ(setup::SigGenSetup) = getδτ(setup.final_charge_size, setup.final_vel)

function getδτ!(setup::SigGenSetup, pos::NTuple{3, Real})
    get_signal!(setup, pos)
    return getδτ(setup)
end
