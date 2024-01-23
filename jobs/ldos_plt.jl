using Distributed
@everywhere using LatticeQM
using Plots

# Plotting bands to ensure that the correct band lies in the energy
# window for the ldos calculation
function showband(H, lat, bandwidth)
	showband(H, lat, -bandwidth / 2, bandwidth / 2)
end

function showband(H, lat, e_low, e_high)
	ks = kpath(lat, ["γ", "κ", "μ", "κ'", "γ"]; num_points=100)
	bands = getbands(H, ks; format=:sparse, num_bands=12)
	#bands.bands = 0.12 / abs(tz) .* bands.bands
	plot(ylim=(e_low, e_high), title="Band Plot")
	#plot(ylim=(-energy, energy), title="Band Plot")
	scatter!(transpose(bands.bands))
	savefig("bandplot.png")
end

function get_ldos(H, klin, bandwidth)
	get_ldos(H, klin, -bandwidth / 2, bandwidth / 2)
end


function get_ldos(H, klin, e_low, e_high)
	nbins = klin ^ 2
	Γ = abs(e_low - e_high) ^ 2
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
	scatter!(positions[1, :], positions[2, :]; markersize=2, zcolor=colors, color=:thermal, clims=(0.0, 0.05))
	savefig(outpath)
end
