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
