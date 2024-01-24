include("setup.jl")
@info "Calculating the dos..."
#avgspacing = (e2 - e1) * tps / tp / klin^2 / num_bands / 5

#gamma = avgspacing / 2 / sqrt(exp(1)-1)

Operators.addchemicalpotential!(H, -mu)
frequencies, dos = Spectrum.getdos(H, -T / 2, T / 2, num_bins; Γ=gamma, klin=klin, format=:sparse, num_bands=num_bands)
plt = plot(title="Density of States")
bar!(frequencies, dos)
xlabel!("Energy")
ylabel!("DOS")
fname = "DOS.png"
savefig(plt, joinpath(OUTPUTDIR, fname))
