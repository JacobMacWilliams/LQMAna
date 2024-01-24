include("setup.jl")
@info "Calculating the magnetization resolved bandstructure at the fermi surface..."
path = kpath(lat, ["γ", "κ", "μ", "κ'", "γ"]; num_points = 100)
Operators.addchemicalpotential!(H, -mu)
sz = Operators.spin(lat, "sz")
bands = getbands(H, path, [sz]; format=:sparse, num_bands = num_bands)
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

fname = joinpath(OUTPUTDIR, "bands.png")
savefig(plta, fname)
