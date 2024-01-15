module effectivebias
using ..LQMRunner
using LatticeQM
using Plots

function getfigS3(micromodel, macromodel, outputdir)

    outputpath = joinpath(outputdir, "effbias.png")
	
	H_micro = gethamiltonian(micromodel)
	lat_micro = getlattice(micromodel)

	H_macro = gethamiltonian(macromodel)
	lat_macro = getlattice(macromodel)

    ks = kpath(lat_micro, ["γ", "κ", "μ", "κ'", "γ"]; num_points = 100)

    valley = Operators.valley(lat_micro; spinhalf=false)
    bands_micro = getbands(H_micro, ks, [valley]; format = :sparse, num_bands = 16)
    tp = getparams(micromodel)[:tperp]
    tps = getparams(micromodel)[:tperp_scaled]
    bands_micro.bands = (tp / tps) * bands_micro.bands

    valley = Operators.spin(lat_macro, "sz")
    bands_macro = getbands(H_macro, ks, [valley])

	plt = plot(bands_micro; ylims=(-0.005, 0.005), colorbar=true, colorbar_title="valley")
	plot!(bands_macro)
	
	savefig(plt, outputpath)
end
end
