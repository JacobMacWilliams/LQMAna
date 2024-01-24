include("setup.jl")
@info "Calculating the dos..."
#avgspacing = (e2 - e1) * tps / tp / klin^2 / num_bands / 5

#gamma = avgspacing / 2 / sqrt(exp(1)-1)

Operators.addchemicalpotential!(H, -mu)
frequencies, dos = Spectrum.getdos(H, -T / 2, T / 2, num_bins; Γ=gamma, klin=klin, format=:sparse, num_bands=num_bands)

# This operation has to be undone here be spinmap.jl
# expects the chemical potential hasn't been
# embedded within the hamiltonian. This behaviour can't
# be changed without changing LatticeQM because
# densitymatrix_distributed is called in spinmap and
# it is this function that ultimately makes this
# assumption. Another work around would be to add the
# chemical potential and then have spinmap assume the
# chemical potential is zero. But that manipulates state
# throughout the length of the program and I prefer this
# because it doesn't.
Operators.addchemicalpotential!(H, mu)

plt = plot(title="Density of States")
bar!(frequencies, dos)
xlabel!("Energy")
ylabel!("DOS")
fname = "DOS.png"
savefig(plt, joinpath(OUTPUTDIR, fname))
