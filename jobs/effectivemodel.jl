module effectivemodel
#using LinearAlgebra
using ..LQMRunner
using LatticeQM
using Plots

function getfigS2(micromodel, macromodel, outputdir)

    outputpath = joinpath(outputdir, "pllfit.png")

	H_micro = gethamiltonian(micromodel)
	lat_micro = getlattice(micromodel)

	H_macro = gethamiltonian(macromodel)

    ks = kpath(lat_micro, ["γ", "κ", "μ", "κ'", "γ"]; num_points = 100)
    bands_micro = getbands(H_micro, ks, format = :sparse, num_bands = 14)
    tp = getparams(micromodel)[:tperp]
    tps = getparams(micromodel)[:tperp_scaled]
    bands_micro.bands = (tp / tps) * bands_micro.bands
    bands_macro = getbands(H_macro, ks)

    plt = plot(bands_micro)
    plot!(bands_macro; markercolor=:red, linestyle=:dash)
    savefig(plt, outputpath)
end
end
