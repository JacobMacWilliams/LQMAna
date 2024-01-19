module GetFig2Wolf2019
using ..LQMRunner
using LatticeQM: getbands, kpath, Operators, Lattices
using LatticeQM.Structure: regulargrid
using LinearAlgebra: I, kron, Diagonal
using Distributed
using Plots
using ProgressMeter

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

function localspin_distributed(ρ, sz, sites::Vector{Int}, latsize::Int, results::RemoteChannel, pool::RemoteChannel)
	spinzs = []
	for i in sites
		proj_i = kron(Diagonal([j == i for j in 1:latsize]), Matrix(1.0I, 2, 2))
		spinz = real(Operators.trace(ρ, sz * proj_i))
		push!(spinzs, spinz)
	end
	put!(results, (sites, spinzs))
	put!(pool, myid())
end

function getfig2b(model, mu, ks)
	
	n = nworkers()
	if !(n > 1)
		error("This code must be run distributed.")
	end
	
	lat = getlattice(model)
    sitenumber = Lattices.countorbitals(lat)
	basis = Lattices.getA(lat)
    latpos = Lattices.positions(lat)

	# we will be translating the spin map by bravais vectors to see the emergent spin
	# texture on the full lattice.
	translations = [zeros(3), basis[:, 1], basis[:, 2], basis[:, 2] - basis[:, 1]]
	latpos = hcat([latpos .+ t for t in translations]...)

	# Recalculating density matrix from Hamiltonian
	ρ = getdensitymatrix(model, mu, ks; multimode=:distributed)
    sz = Operators.spin(lat, "sz")
	
	# Parsing latticeidx chunks for scheduling
	chunksize = ceil(sitenumber / n)
	sites = [i for i in 1:sitenumber]
	chunks = []
	for i in 1:n
		startidx = Int(chunksize * (i - 1) + 1)
		endidx = Int(chunksize * i)
		if endidx > sitenumber
			endidx = sitenumber
		end
		push!(chunks, sites[startidx:endidx])
	end
	
	# asynchronously schedule matrix multiplication on the workers and greedily scatter
	# results as as soon as availabe in the results channel. The worker pool guarantees that
	# the individual workers aren't "overscheduled/overworked".
	results = RemoteChannel()
	pool  = RemoteChannel(() -> Channel{Int}(n))
	for i in workers()
		put!(pool, i)
	end

	plt = plot()
	@sync begin

		# Greedy processing
		@async @showprogress "Plotting spin map..." for _ in eachindex(chunks)
			siteidxs, spins = take!(results)
			println(size(spins))
			# copying results for translated unit cells
			siteidxs = vcat([siteidxs .+ (i - 1) * sitenumber for i in eachindex(translations)]...)
			spins = vcat([spins for i in eachindex(translations)]...)
			scatter!(latpos[1,siteidxs], latpos[2, siteidxs]; markersize=2, zcolor=spins, color=:roma)
		end
		
		# Scheduling
		for chunk in chunks
			worker = take!(pool)
			remote_do(localspin_distributed, worker, ρ, sz, chunk, sitenumber, results, pool)
		end
	end
	close(results)
    return plt
end

end
