module FlatBand
using FromFile
#using Distributed
#nprocs() < 2 ? addprocs(3) : nothing
using LatticeQM: kpath, getbands, Lattices.twistangle as twistangle
@from "../../utils/BandPlots.jl" using BandPlots
using ..LQMRunner
#import Statistics
#import Optim
#using Plots
#export getunscaledtwist

getunscaledtwist(theta::Float64, tperp::Float64) = asind(sind(theta / 2) * 0.12 / tperp) * 2
getunscaledtwist(n::Int64, tperp::Float64) =
    getunscaledtwist(twistangle(n) * 180 / pi, tperp)

function getfig1df(model, outputdir::AbstractString)
	n = getparams(model)[:n]
	tp = getparams(model)[:tperp]
	tps = getparams(model)[:tperp_scaled]
	
	alpha = getunscaledtwist(n, tps)
	astring = string(round(alpha; digits = 3))
	filename = "twistedbands_" * astring * ".png"
	filepath = joinpath(outputdir, filename)
	
	H = gethamiltonian(model)
	lat = getlattice(model)
	ks = kpath(lat, ["γ", "κ", "μ", "κ'", "γ"]; num_points = 100)
	bands = getbands(H, ks; format=:sparse, num_bands = 14, multimode=:distributed)
	bands.bands = (tp / tps) .* bands.bands
	bandplotsave(bands, filepath)
end

#=
function findflatbands()
    #t_perp = 0.37
    n = 14
    tperp = Optim.optimize(t -> flatloss(n, t), 0.2, 0.6)
    return tperp
end

function flatloss(n, t_perp)
    bands = twistedbands(n, t_perp)
    pll1 = bands.bands[5, :]
    pll2 = bands.bands[7, :]

    mean1 = Statistics.mean(pll1)
    score1 = sum(broadcast(abs, pll1 .- mean1))
    mean2 = Statistics.mean(pll2)
    score2 = sum(broadcast(abs, pll2 .- mean2))
    score = score1 + score2
    return score
end

toptobj = findflatbands()
topt = Optim.minimizer(toptobj)
Optim.summary(toptobj)
bands = twistedbands(14, topt)
p1 = plot(bands; size=(300,250), ylim=(-0.02, 0.02), markercolor=:diverging_bkr_55_10_c35_n256, colorbar=true)

!isdir("output/tsbg") ? mkpath("output") : nothing
filename = "twistedbands_" * string(round(getunscaledtwist(14, topt); digits=3)) * ".png"
savefig("output/"*filename)
=#

#getfig1df()
end
