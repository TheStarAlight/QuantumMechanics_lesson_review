using Plots
pyplot()

include("TDSE_1D.jl")

V(x) = ifelse(abs(x)<=10.0, 0.0, 0.4)
Nx = 4096
dx = 0.05
dt = 0.02
job_im = TDSEJob_1D_SOF_ImTime(Nx=Nx,dx=dx,dt=dt, V0=V, N_ex=6)
prop!(job_im)

ψ_scale_shift(job, i, scale) = real(job.ψ[i]) .* scale .+ job.E0[i]
scale = 0.05

plot(job_im.x, ψ_scale_shift(job_im, 1, scale))
plot!(job_im.x, ψ_scale_shift(job_im, 2, scale))
plot!(job_im.x, ψ_scale_shift(job_im, 3, scale))
plot!(job_im.x, ψ_scale_shift(job_im, 4, scale))
plot!(job_im.x, ψ_scale_shift(job_im, 5, scale))
plot!(job_im.x, ψ_scale_shift(job_im, 6, scale))
plot!(job_im.x, V, linecolor=:black)
annotate!(-10.0, -0.01, text(raw"$V(x)$", :left, :top))
annotate!(-14.0, job_im.E0[1], text(raw"$ψ_1(x)$", :left, :bottom))
annotate!(-14.0, job_im.E0[2], text(raw"$ψ_2(x)$", :left, :bottom))
annotate!(-14.0, job_im.E0[3], text(raw"$ψ_3(x)$", :left, :bottom))
annotate!(-14.0, job_im.E0[4], text(raw"$ψ_4(x)$", :left, :bottom))
annotate!(-14.0, job_im.E0[5], text(raw"$ψ_5(x)$", :left, :bottom))
annotate!(-14.0, job_im.E0[6], text(raw"$ψ_6(x)$", :left, :bottom))
plot!([-15.0, 15.0], [job_im.E0[1], job_im.E0[1]], linecolor=:black, alpha=0.2, linestyle=:dash)
plot!([-15.0, 15.0], [job_im.E0[2], job_im.E0[2]], linecolor=:black, alpha=0.2, linestyle=:dash)
plot!([-15.0, 15.0], [job_im.E0[3], job_im.E0[3]], linecolor=:black, alpha=0.2, linestyle=:dash)
plot!([-15.0, 15.0], [job_im.E0[4], job_im.E0[4]], linecolor=:black, alpha=0.2, linestyle=:dash)
plot!([-15.0, 15.0], [job_im.E0[5], job_im.E0[5]], linecolor=:black, alpha=0.2, linestyle=:dash)
plot!([-15.0, 15.0], [job_im.E0[6], job_im.E0[6]], linecolor=:black, alpha=0.2, linestyle=:dash)
plot!(xlim=(-15,15), ylim=(-0.05,0.45), legend=:none, framestyle=:box, xlabel=raw"position $x$", ylabel=raw"energy $E$")

savefig("figures/sec3_figures/finite_well_eigenstates.pdf")
