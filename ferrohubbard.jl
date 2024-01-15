using Distributed
using Plots: savefig
if !(nprocs() > 1)
	error("This script should run as distributed")
end

@everywhere using LatticeQM
using LQMRunner
using JLD2

const PROJECTROOT = pkgdir(LQMRunner)
const modelpath = joinpath(PROJECTROOT, "models", ARGS[1] * ".jld2")
if !isfile(modelpath)
	error("File not found")
end

model = load(modelpath)["model"];

include("jobs/GetFig2Wolf2019.jl")
using .GetFig2Wolf2019

ks = Structure.regulargrid(;nk=81)
mu = getchemicalpotential(model, -6, ks)

@time plta = GetFig2Wolf2019.getfig2a(model, mu)
@time pltb = GetFig2Wolf2019.getfig2b(model, mu, ks)
savefig(plta, "distributeoptbands.png")
savefig(pltb, "distributedoptspinmap.png")
println(mu)

