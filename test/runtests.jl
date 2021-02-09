using Test, MJDSigGen, DelimitedFiles

@testset "MJDSigGen.jl" begin
    mktempdir() do scratchdir
        cp(joinpath("..", "examples"), joinpath(scratchdir, "examples"))

        scratchdir = abspath(joinpath(dirname(pathof(MJDSigGen)), "../"))

        config_file_rel_path = joinpath("examples", "config", "example.config")
        config_file = joinpath(scratchdir, config_file_rel_path)

        # Run fieldgen in scratchdir, as it writes "undepleted.txt" into current directory
        # cd(scratchdir) do
        #     fieldgen(config_file_rel_path)
        # end

        @testset "read_config" begin
            let setup1 = SigGenSetup(), setup2 = read_config(config_file)
                read_config!(setup1, config_file)

                @test setup1 !== setup2
                
                # Rquire all fields set to 0 initially?

                # test some fields
                for setup in (setup1, setup2)
                    setup isa SigGenSetup

                    @test setup.xtal_length         == 80
                    @test setup.xtal_radius         == 35
                    @test setup.bottom_taper_length == 0
                    @test setup.outer_taper_length  == 60

                    @test setup.xtal_grid      == 0.1f0
                    @test setup.xtal_HV        == 3500
                    @test setup.max_iterations == 30000

                    @test setup.config_name == config_file
                    @test setup.drift_name  == "./drift_vel_tcorr.tab"
                    @test setup.field_name  == "../fields/ev_example.dat"
                    @test setup.wp_name     == "../fields/wp_example.dat"

                    @test setup.xtal_temp      == 90
                    @test setup.step_time_calc == 1.0
                    @test setup.use_diffusion  == 0
                end
            end
        end

        @testset "fieldgen and read_fields" begin
            field_data, wpot_data = cd(joinpath(scratchdir, "examples", "fields")) do 
                readdlm("ev_example.dat", comments=true),
                readdlm("wp_example.dat", comments=true)
            end

            setup = read_config(config_file)

            # Number of points
            Nr::Int = setup.xtal_radius / setup.xtal_grid + 1
            Nz::Int = setup.xtal_length / setup.xtal_grid + 1

            @test size(field_data, 1) == size(wpot_data, 1) == Nr * Nz
            @test size(field_data, 2) == 6
            @test size(wpot_data, 2)  == 3

            # 0 <= radius <= max radius
            @test all(r -> 0 <= r <= setup.xtal_radius, view(field_data, :, 1))
            @test all(r -> 0 <= r <= setup.xtal_radius, view(wpot_data, :, 1))

            # 0 <= z <= length
            @test all(z -> 0 <= z <= setup.xtal_length, view(field_data, :, 2))
            @test all(z -> 0 <= z <= setup.xtal_length, view(wpot_data, :, 2))

            # 0 <= wp <= 1
            @test all(p -> 0 <= p <= 1, view(wpot_data, :, 3))

            # E^2 = E_z^2 + E_r^2 (within 5%)
            @test all(eachrow(view(field_data, :, 4:6))) do v 
                E, E_r, E_z = v
                d = abs((E^2 - E_r^2 - E_z^2))
                return E == 0 ? d ≈ 0 : sqrt(d / E^2) < 0.05
            end

            # read_fields
            E_pot, W_pot, E_abs, E_r, E_z = read_fields(setup)
            @test E_pot == reshape(field_data[:, 3], (Nz, Nr))
            @test W_pot == reshape(wpot_data[:, 3], (Nz, Nr))
            @test E_abs == reshape(field_data[:, 4], (Nz, Nr))
            @test E_r   == reshape(field_data[:, 5], (Nz, Nr))
            @test E_z   == reshape(field_data[:, 6], (Nz, Nr))
        end

        @testset "cyl2cart and cart2cyl" begin
            carts = [(1.,  1., 1.), (1., 0., -5.), (0.,   -2., 3.)]
            cyls  = [(√2, π/4, 1.), (1., 0., -5.), (2., -π/2., 3.)]

            for ((x, y, z1), (r, ϕ, z2)) in zip(carts, cyls)
                @test collect(MJDSigGen.cart2cyl(x, y, z1)) ≈ [r, ϕ, z2]
                @test collect(MJDSigGen.cyl2cart(r, ϕ, z2)) ≈ [x, y, z1]
            end
        end

        # @testset "geometry" begin
        #     @test MJDSigGen.outside_detector(setup, (-10, 0, -10)) == true
        #     @test MJDSigGen.outside_detector(setup, (10, 0, 10)) == false

        #     grid = setup.xtal_grid
        #     @test MJDSigGen.nearest_field_grid_index(setup, (8, 0, 12)) == (:interpol, round.(Int, [12, 8] / grid .+ 1)...)
        #     @test MJDSigGen.nearest_field_grid_index(setup, (10, 0, -10)) == (:outside,0,0)
        #     @test MJDSigGen.nearest_field_grid_index(setup, (1.0 * setup.xtal_radius, 0.0, 10.0)) == (:extrapol, round.(Int, [10.0, 1.0 * setup.xtal_radius] / grid + [1, 0])...)
        # end


        # @testset "siggen" begin
        #     charge_signal = Array{Float32, 1}()
        #     @test begin
        #         charge_signal = MJDSigGen.get_signal!(setup, (10.1, 10.0, 10.0))::typeof(charge_signal)
        #         size(charge_signal) == (setup.ntsteps_out,)
        #     end

        #     path_e = Array{Float32}(undef, 0,0)
        #     @test begin
        #         path_e = MJDSigGen.drift_path(setup, :e)::typeof(path_e)
        #         (size(path_e, 2) == 3) && (0 < size(path_e, 1) <= setup.time_steps_calc)
        #     end

        #     path_h = Array{Float32}(undef, 0,0)
        #     @test begin
        #         path_h = MJDSigGen.drift_path(setup, :h)::typeof(path_h)
        #         (size(path_h, 2) == 3) && (0 < size(path_h, 1) <= setup.time_steps_calc)
        #     end

        #     @test 0 <= charge_signal[1] < 0.01
        #     n_steps_drift_out = ceil(Int, (max(size(path_e, 1), size(path_h, 1)) / setup.time_steps_calc * setup.ntsteps_out)) + 1
        #     @test charge_signal[n_steps_drift_out] ≈ 1

        #     @test begin
        #         MJDSigGen.signal_calc_finalize!(setup)
        #         result = (setup.wpot == C_NULL && setup.dpath_h == C_NULL)
        #         MJDSigGen.signal_calc_finalize!(setup)
        #         result
        #     end
        # end
    end
end
