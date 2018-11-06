using BinDeps
using Base.Sys: isapple,isbsd,islinux,isunix,iswindows

@BinDeps.setup

libmjd_siggen = library_dependency("libmjd_siggen", aliases=[])

prefix = joinpath(BinDeps.depsdir(libmjd_siggen), "usr")
srcdir = joinpath(BinDeps.srcdir(libmjd_siggen), "mjd_siggen")
builddir = joinpath(BinDeps.builddir(libmjd_siggen), "mjd_siggen")

# BSD systems (other than macOS) use BSD Make rather than GNU Make by default
# We need GNU Make, and on such systems GNU make is invoked as `gmake`
make = isbsd() && !isapple() ? "gmake" : "make"

provides(SimpleBuild,
    (@build_steps begin
        if isdir(builddir)
            rm(builddir, recursive = true)
        end
        CreateDirectory(builddir)
        @build_steps begin
            ChangeDirectory(builddir)
            `cmake -DCMAKE_INSTALL_PREFIX="$prefix" $srcdir`
            `$make -j$(min(Sys.CPU_THREADS, 8)) install`
        end
    end), [libmjd_siggen], os = :Unix)

@BinDeps.install Dict(:libmjd_siggen => :libmjd_siggen)
