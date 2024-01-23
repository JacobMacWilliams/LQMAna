using Distributed
using Plots
using JLD2
using ProgressMeter

const N = nworkers()
if !(N > 1)
	error("This code must be run distributed.")
end

@everywhere using LatticeQM
@everywhere using LQMRunner
@everywhere import LinearAlgebra.BLAS
@everywhere using LinearAlgebra: kron, Diagonal, I
@everywhere BLAS.set_num_threads(1)

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

function get_ldos(H, klin, bandwidth)
	get_ldos(H, klin, -bandwidth / 2, bandwidth / 2)
end

function get_ldos(H, klin, e_low, e_high)
	nbins = klin ^ 2
	#Γ = abs(e_low - e_high) ^ 2
	Γ = abs(e_low - e_high) / 2 / nbins / sqrt(exp(1) - 1)
	get_ldos(H, klin, e_low, e_high, nbins, Γ)
end

function get_ldos(H, klin, e_low, e_high, nbins, Γ)
	# Calculating the ldos for the appropriate pll band
	ks = Structure.regulargrid(nk=klin^2)
	frequencies = [e_low + e_high * n / nbins for n in 0:nbins]
	#display(frequencies) # Sanity check
	ldos = Spectrum.ldos(H, ks, frequencies; Γ=Γ, num_bands=12)
	return ldos
end
	
	
function latheatmap(lat, colors, title, outpath)
	# Plotting the ldos for multiple unit cells to better display structure
	# resulting from patterns on the edge.
	positions = Lattices.positions(lat)
	basis = Lattices.getA(lat)
	translations = [zeros(size(basis[:, 1])), basis[:, 1], basis[:, 2], basis[:, 2] - basis[:, 1]]
	
	positions = hcat([positions .+ b for b in translations]...)
	colors = vcat([colors for i in 1:length(translations)]...)
	
	plot(title=title)
	#scatter!(positions[1, :], positions[2, :]; markersize=2, zcolor=colors, color=:thermal, clims=(0.0, 0.05))
	scatter!(positions[1, :], positions[2, :]; markersize=2, zcolor=colors, color=:thermal)
	savefig(outpath)
end

const MODELNAME = ARGS[1]

@info "Loading model " * MODELNAME * "..."
const OUTPUTDIR = joinpath("output", MODELNAME * "_ana")
if !isdir(OUTPUTDIR)
	mkpath(OUTPUTDIR)
end

#const PROJECTROOT = pkgdir(LQMRunner)
const MODELPATH = joinpath(pkgdir(LQMRunner), "models", MODELNAME * ".jld2")
if !isfile(MODELPATH)
	error("File not found")
end

model = load(MODELPATH)["model"];
lat = getlattice(model)
H = gethamiltonian(model)
#klin = getparams(model)[:klin]
klin = 2
T = getparams(model)[:T]
tp = getparams(model)[:tperp]
tps = getparams(model)[:tperp_scaled]
#doping = getparams(model)[:doping]
doping = -6

@info "Calculating chemical potential at " * string(doping) * " charge doping..."
ks = Structure.regulargrid(;nk=klin^2)
mu = getchemicalpotential(model, doping, ks; multimode=:distributed)

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

fname = joinpath(OUTPUTDIR, "distributedoptspinmap.png")
savefig(pltb, fname)

@info "Calculating the magnetization resolved bandstructure at the fermi surface..."
ks = kpath(lat, ["γ", "κ", "μ", "κ'", "γ"]; num_points = 100)
Operators.addchemicalpotential!(H, -mu)
bands = getbands(H, ks, [sz]; format=:sparse, num_bands = 24)
bands.bands = (tp / tps) .* bands.bands

# Plotting band structure
e1 = min(bands.bands...)
e2 = max(bands.bands...)
plta = plot(
    bands;
    size = (300, 250),
    ylim = (e1, e2),
    markercolor = :diverging_bkr_55_10_c35_n256,
    colorbar = true,
)

fname = joinpath(OUTPUTDIR, "distributedoptbands.png")
savefig(plta, fname)

@info "Calculating the ldos..."
#ldos = get_ldos(H, klin, -0.01, 0.01, 100, 0.01)
ldos = get_ldos(H, klin, T)
fname = "test_ldos.png"
latheatmap(lat, ldos, "ldos plot", joinpath(OUTPUTDIR, fname))

sleep(2) # Allow workers to close
println(mu)
