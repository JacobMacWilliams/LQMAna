include("setup.jl")
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


@info "Calculating the ldos..."
#ldos = get_ldos(H, klin, -T / 2, T / 2, 100, 0.01)

frequencies = LinRange(-T / 2, T / 2, num_bins)
Operators.addchemicalpotential!(H, -mu)
ldos = Spectrum.ldos(H, ks, collect(frequencies); Γ=gamma, format=:sparse, num_bands=num_bands)

# See comment in dos.jl
Operators.addchemicalpotential!(H, mu)

fname = "ldos.png"
latheatmap(lat, ldos, "ldos plot", joinpath(OUTPUTDIR, fname))
