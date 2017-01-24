# This file is a part of MJDSigGen, licensed under the MIT License (MIT).

export signal_calc_init
export get_signal, get_signal!


function signal_calc_init(config_filename::AbstractString)
     setup = Struct_MJD_Siggen_Setup()

    r = ccall(
        @sgsym(:signal_calc_init),
        Cint,
        (Cstring, Ptr{Struct_MJD_Siggen_Setup}),
        pointer(config_filename), Ref(setup)
    )
    r != 0 && error("signal_calc_init failed.")
    setup
end


function get_signal!(signal, setup::Struct_MJD_Siggen_Setup, location::NTuple{3})
    (length(linearindices(signal)) < setup.ntsteps_out) && throw(BoundsError())

    pt = Struct_point(location[1], location[2], location[3])

    r = ccall(
        @sgsym(:get_signal),
        Cint,
        (Struct_point, Ptr{Float32}, Ptr{Struct_MJD_Siggen_Setup}),
        pt, signal, Ref(setup)
    )
    r < 0 && error("Point not in crystal or has no field: $pt")

    signal
end


get_signal(setup::Struct_MJD_Siggen_Setup, location::NTuple{3}) =
    get_signal!(zeros(Float32, setup.ntsteps_out), setup, location)
