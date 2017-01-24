# This file is a part of MJDSigGen, licensed under the MIT License (MIT).

function tuplestr{N}(s::NTuple{N,Cchar})
    a = [c % UInt8 for c in s]
    unsafe_string(pointer(a))
end
