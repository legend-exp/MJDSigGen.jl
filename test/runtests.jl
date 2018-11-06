using MJDSigGen
using Test

@testset "Package MJDSigGen" begin
    mktempdir() do scratchdir
        cp(joinpath("..", "examples"), joinpath(scratchdir, "examples"))

        config_file_rel_path = joinpath("examples", "config", "example.config")
        config_file = joinpath(scratchdir, config_file_rel_path)

        setup = MJDSigGen.Struct_MJD_Siggen_Setup()


        @testset "config" begin
            @test (setup = MJDSigGen.read_config(config_file); true)
        end


        @testset "util" begin
            @testset "file names" begin
                @test MJDSigGen.config_file_name(setup) == config_file
                @test length(MJDSigGen.field_file_name(setup)) > 0
                @test length(MJDSigGen.wpot_file_name(setup)) > 0
            end

            @testset "cyl2cart and cart2cyl" begin
                @test [MJDSigGen.cyl2cart(MJDSigGen.cart2cyl(10.1, 12.1, 2.1)...)...] ≈ [10.1, 12.1, 2.1]
                @test [MJDSigGen.cyl2cart(MJDSigGen.cart2cyl(-10.2, 12.2, 2.2)...)...] ≈ [-10.2, 12.2, 2.2]
                @test [MJDSigGen.cyl2cart(MJDSigGen.cart2cyl(-10.3, -12.3, 2.3)...)...] ≈ [-10.3, -12.3, 2.3]
                @test [MJDSigGen.cyl2cart(MJDSigGen.cart2cyl(10.4, -12.4, 2.4)...)...] ≈ [10.4, -12.4, 2.4]
            end
        end


        @testset "fieldgen" begin
            # Run fieldgen in scratchdir, as it writes "undepleted.txt" into current directory:
            cd(scratchdir) do
                @test (MJDSigGen.fieldgen(config_file_rel_path); true)
            end

            @test (setup = MJDSigGen.signal_calc_init(config_file); true)

            E_pot, W_pot, E_abs, E_r, E_z = MJDSigGen.read_fields(setup)

            @test (
                (setup.zlen, setup.rlen) == size(E_pot) == size(W_pot)
                == size(E_abs) == size(E_r) == size(E_z)
            )
        end


        @testset "geometry" begin
            @test MJDSigGen.outside_detector(setup, (-10, 0, -10)) == true
            @test MJDSigGen.outside_detector(setup, (10, 0, 10)) == false

            grid = setup.xtal_grid
            @test MJDSigGen.nearest_field_grid_index(setup, (8, 0, 12)) == (:interpol, round.(Int, [12, 8] / grid .+ 1)...)
            @test MJDSigGen.nearest_field_grid_index(setup, (10, 0, -10)) == (:outside,0,0)
            @test MJDSigGen.nearest_field_grid_index(setup, (1.0 * setup.xtal_radius, 0.0, 10.0)) == (:extrapol, round.(Int, [10.0, 1.0 * setup.xtal_radius] / grid + [1, 0])...)
        end


        @testset "siggen" begin
            charge_signal = Array{Float32, 1}()
            @test begin
                charge_signal = MJDSigGen.get_signal!(setup, (10.1, 10.0, 10.0))::typeof(charge_signal)
                size(charge_signal) == (setup.ntsteps_out,)
            end

            path_e = Array{Float32}(undef, 0,0)
            @test begin
                path_e = MJDSigGen.drift_path(setup, :e)::typeof(path_e)
                (size(path_e, 2) == 3) && (0 < size(path_e, 1) <= setup.time_steps_calc)
            end

            path_h = Array{Float32}(undef, 0,0)
            @test begin
                path_h = MJDSigGen.drift_path(setup, :h)::typeof(path_h)
                (size(path_h, 2) == 3) && (0 < size(path_h, 1) <= setup.time_steps_calc)
            end

            @test 0 <= charge_signal[1] < 0.01
            n_steps_drift_out = ceil(Int, (max(size(path_e, 1), size(path_h, 1)) / setup.time_steps_calc * setup.ntsteps_out)) + 1
            @test charge_signal[n_steps_drift_out] ≈ 1

            @test begin
                MJDSigGen.signal_calc_finalize!(setup)
                result = (setup.wpot == C_NULL && setup.dpath_h == C_NULL)
                MJDSigGen.signal_calc_finalize!(setup)
                result
            end
        end
    end
end
