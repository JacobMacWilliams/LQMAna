using Distributed
using JLD2
using Plots

const N = nworkers()
if !(N > 1)
	error("This code must be run distributed.")
end

@everywhere using LatticeQM
@everywhere using LQMRunner
@everywhere import LinearAlgebra.BLAS
@everywhere BLAS.set_num_threads(1)

const MODELNAME = ARGS[1]

@info "Loading model " * MODELNAME * "..."
const OUTPUTDIR = joinpath(@__DIR__, "..", "output", MODELNAME * "_ana")
if !isdir(OUTPUTDIR)
	mkpath(OUTPUTDIR)
end

#const PROJECTROOT = pkgdir(LQMRunner)
const MODELPATH = joinpath(pkgdir(LQMRunner), "models", MODELNAME * ".jld2")
if !isfile(MODELPATH)
	error("File not found")
end


const model = load(MODELPATH)["model"];
const lat = getlattice(model)
const H = gethamiltonian(model)
const klin = getparams(model)[:klin]
const tp = getparams(model)[:tperp]
const tps = getparams(model)[:tperp_scaled]
const ks = Structure.regulargrid(;nk=klin^2)
const T = getparams(model)[:T]
#const doping = getparams(model)[:doping]
const doping = -6
const num_bands = 24
const num_bins = 100
const gamma = 0.005

@info "Calculating chemical potential at " * string(doping) * " charge doping..."
const mu = getchemicalpotential(model, doping, ks; multimode=:distributed)

failed = 0
for f in ARGS[2:end]
	if !isfile(joinpath(@__DIR__, f))
		f = f * ".jl"
		if !isfile(joinpath(@__DIR__, f))
			global failed += 1
			continue
		end
	end
	include(f)
end

if failed == 0
	println("Jobs ran.")
else
	println(string(failed) * " jobs failed to run.")
end
