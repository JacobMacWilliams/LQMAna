module GetFig2Wolf2019
using ..LQMRunner
using LatticeQM: getbands, kpath, Operators, Lattices
using LatticeQM.Structure: regulargrid
using LinearAlgebra
using Plots

function getfig2(model, outputdir)
	
	name = getname(model)
    outputpatha = joinpath(outputdir, name * "_bands.png")
    outputpathb = joinpath(outputdir, name * "_spinmap.png")
	
	# Calculating the chemical potential at PLL-1 half filling
	klin = getparams(model)[:klin]
	ks = regulargrid(nk=klin)
	mu = getchemicalpotential(model, -6, ks)
    
	plta = getfig2a(model, mu)
    pltb = getfig2b(model, mu, ks)
    savefig(plta, outputpatha)
    savefig(pltb, outputpathb)
end

function getfig2a(model, mu)
	
	H = gethamiltonian(model)

	Operators.addchemicalpotential!(H, -mu)

	# Retrieving band structure along high symmetry path
	lat = getlattice(model)
    ks = kpath(lat, ["γ", "κ", "μ", "κ'", "γ"]; num_points = 100)
	sz = Operators.spin(lat, "sz")
	bands = getbands(H, ks, [sz]; format=:sparse, num_bands = 24)
	tp = getparams(model)[:tperp]
	tps = getparams(model)[:tperp_scaled]
    bands.bands = (tp / tps) .* bands.bands
    
	# Plotting band structure
	e1 = min(bands.bands...)
    e2 = max(bands.bands...)
    plt = plot(
        bands;
        size = (300, 250),
        ylim = (e1, e2),
        markercolor = :diverging_bkr_55_10_c35_n256,
        colorbar = true,
    )
    return plt
end

function getfig2b(model, mu, ks)
	
	lat = getlattice(model)
	ρ = getdensitymatrix(model, mu, ks; multimode=:distributed)
    
	basis = Lattices.getA(lat)
    sitenumber = Lattices.countorbitals(lat)
    latpos = Lattices.positions(lat)
    sz = Operators.spin(lat, "sz")

    siteprojector(latsize, i) = kron(Diagonal([j == i for j in 1:latsize]), Matrix(1.0I, 2, 2))
	
    expspinzs = []
    for i in 1:sitenumber
        proj_i = siteprojector(sitenumber, i)
		expspinz = real(Operators.trace(ρ, sz * proj_i))
        push!(expspinzs, expspinz)
    end
	
	translations = [zeros(3), basis[:, 1], basis[:, 2], basis[:, 2] - basis[:, 1]]
	latpos = hcat([latpos .+ b for b in translations]...)
	expspinzs = vcat([expspinzs for i in 1:length(translations)]...)
	plt = plot(scatter(latpos[1, :], latpos[2, :]; zcolor=expspinzs, color=:roma))
    return plt
end

end
