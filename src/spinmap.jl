include("setup.jl")
using ProgressMeter
@everywhere using LinearAlgebra: kron, Diagonal, I

@everywhere function localspin_distributed(ρ, sz, sites::Vector{Int}, latsize::Int, results::RemoteChannel, pool::RemoteChannel)
	spinzs = []
	for i in sites
		proj_i = kron(Diagonal([j == i for j in 1:latsize]), Matrix(1.0I, 2, 2))
		spinz = real(Operators.trace(ρ, sz * proj_i))
		push!(spinzs, spinz)
	end
	put!(results, (sites, spinzs))
	put!(pool, myid())
end

@info "Calculating the spatially resolved magnetization of the ground state..."
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
chunksize = ceil(sitenumber / N)
sites = [i for i in 1:sitenumber]
chunks = []
for i in 1:N
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
pool  = RemoteChannel(() -> Channel{Int}(N))
for i in workers()
	put!(pool, i)
end

pltb = plot()
@sync begin
	# Greedy processing
	@async @showprogress "Plotting spin map..." for _ in eachindex(chunks)
		siteidxs, spins = take!(results)
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

fname = joinpath(OUTPUTDIR, "spinmap.png")
savefig(pltb, fname)
