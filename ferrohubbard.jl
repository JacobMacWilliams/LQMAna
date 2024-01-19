using Distributed
using Plots: savefig
if !(nprocs() > 1)
	error("This script should run as distributed")
end

@everywhere using LatticeQM
@everywhere import LinearAlgebra.BLAS
#@everywhere using LinearAlgebra: Diagonal
@everywhere using LQMRunner
@everywhere include("jobs/GetFig2Wolf2019.jl")
@everywhere using .GetFig2Wolf2019
BLAS.set_num_threads(1)
using JLD2

const PROJECTROOT = pkgdir(LQMRunner)
const modelpath = joinpath(PROJECTROOT, "models", ARGS[1] * ".jld2")
if !isfile(modelpath)
	error("File not found")
end

model = load(modelpath)["model"];
klin = getparams(model)[:klin]
#doping = getparams(model)[:doping]
doping = -6


ks = Structure.regulargrid(;nk=klin^2)

@info "Calculating chemical potential at " * string(doping) * " charge doping..."
mu = getchemicalpotential(model, doping, ks; multimode=:distributed)

@info "Calculating the magnetization resolved bandstructure at the fermi surface..."
@time plta = GetFig2Wolf2019.getfig2a(model, mu)
savefig(plta, "distributedoptbands.png")

@info "Calculating the spatially resolved magnetization of the ground state..."
@time pltb = GetFig2Wolf2019.getfig2b(model, mu, ks)
savefig(pltb, "distributedoptspinmap.png")
println(mu)

