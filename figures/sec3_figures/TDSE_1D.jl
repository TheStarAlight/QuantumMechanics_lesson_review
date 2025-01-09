using FFTW
using ProgressMeter
using Printf

abstract type TDSEJob end
abstract type TDSEJob_1D end

"One-dimensional TDSE job which uses imaginary-time evolution to obtain the ground state using split-operator-Fourier method."
mutable struct TDSEJob_1D_SOF_ImTime
    # Parameters
    Nx::Int             # number of grid points
    dx::Float64         # coord grid size
    dt::Float64         # time step
    conv_ΔE::Float64    # criteria of energy convergence
    N_ex::Int           # number of excited states

    # Calculated Parameters
    x1::Float64         # begin of coord
    x2::Float64         # end of coord
    E0::Vector{Float64} # energies of gnd and excited states

    # Memory
    ψ::Vector{Vector{ComplexF64}}  # wavefunction
    ψ_pur::Vector{Vector{ComplexF64}}  # wavefunction after purification [after executing purify!]
    x::Vector{Float64}  # coords
    p::Vector{Float64}  # momentum coords
    V0::Vector{Float64} # evaluated atomic potential in coordinate space
    T0::Vector{Float64} # kinetic energy in momentum space
    mask::Vector{Float64}   # absorbing boundary mask function
    UVt::Vector{ComplexF64} # [mutable] exp(-iV(x)Δt/2)
    UTt::Vector{ComplexF64} # [mutable] exp(-iT(k)Δt)
    pF::FFTW.Plan  # FFT plan
    pI::FFTW.Plan  # IFFT plan
end

function TDSEJob_1D_SOF_ImTime(;Nx,dx,dt, guess_wfn=(x-> exp(-(x-1)^2)), V0, conv_ΔE=1e-16, N_ex=1)
    @assert N_ex >= 1 "Number of excited states must be at least 1."
    x = [(i-0.5(Nx+1))*dx for i=1:Nx]
    p = 2π*fftfreq(Nx,1.0/dx)

    x1 = x[1]
    x2 = x[end]
    E0 = Float64[]

    V0_ = if typeof(V0) <: Function
        V0.(x)
    else # Vector
        V0 .|> Float64
    end
    T0 = @. p^2/2
    UVt = @. ComplexF64(exp(-0.5*V0_*dt))
    UTt = @. ComplexF64(exp(-T0*dt))
    mask = gen_mask_1D(Nx, 0.8)

    ψ0 = if typeof(guess_wfn) <: Function
        guess_wfn.(x) .|> ComplexF64
    else # Vector
        guess_wfn .|> ComplexF64
    end
    ψ = [copy(ψ0) for i=1:N_ex]
    ψ_pur = []

    pF = FFTW.plan_fft!(ψ[1])
    pI = FFTW.plan_ifft!(ψ[1])

    return TDSEJob_1D_SOF_ImTime(Nx, dx, dt, conv_ΔE, N_ex, x1, x2, E0, ψ, ψ_pur, x, p, V0_, T0, mask, UVt, UTt, pF, pI)
end

function gen_mask_1D(Nx, edge_ratio)
    mask = ones(Nx)
    @inbounds for i in 1:Nx
        u = abs(i-Nx/2)/(Nx/2)
        u0 = edge_ratio
        u>=u0 && (mask[i] = cos((u-u0)/(1-u0+0.01) * 0.5π)^(1/6))
    end
    return mask
end

function prop!(job::TDSEJob_1D_SOF_ImTime)
    x = job.x
    V0 = job.V0
    T0 = job.T0
    dx = job.dx
    ψ = job.ψ
    N_ex = job.N_ex
    UTt = job.UTt
    UVt = job.UVt
    pF = job.pF
    pI = job.pI
    kfac = dx/sqrt(2π)
    dp = job.p[2] - job.p[1]
    E_old = -0.0
    E_new = -1.0
    job.E0 = Float64[]
    prog = ProgressUnknown(desc="Running imag time prop for gnd state...")
    while abs(E_old-E_new) > job.conv_ΔE
        ψ[1] .*= UVt
        pF * ψ[1]
        ψ[1] .*= UTt
        pI * ψ[1]
        ψ[1] .*= UVt
        ψ[1] ./= sqrt((ψ[1] .|> abs2 |> sum)*dx)
        E_old = E_new
        pF * ψ[1]
        E_new = dp*kfac^2*sum(T0.*abs2.(ψ[1]))
        pI * ψ[1]
        E_new += dx*sum(V0.*abs2.(ψ[1]))
        update!(prog, desc=@sprintf("Running imag time prop for gnd state... E=%.10f ΔE=%.2e", E_new, abs(E_new-E_old)), spinner=raw"\-/|")
    end
    finish!(prog)
    println()
    push!(job.E0, E_new)
    if N_ex > 1
        for i in 2:N_ex
            E_old = -0.0
            E_new = -1.0
            @info "Running imag time prop for excited state #$i..."
            prog = ProgressUnknown(desc="Running imag time prop for excited state #$i...")
            while abs(E_old-E_new) > job.conv_ΔE
                for j in 1:i-1
                    ψ[i] .-= sum(ψ[j] .* ψ[i]) * dx .* ψ[j]     # proj to ψ[j] and remove it
                    # ψ[i] ./= sqrt((ψ[i] .|> abs2 |> sum)*dx)    # normalize
                end
                # continue imag time propagation
                ψ[i] .*= UVt
                pF * ψ[i]
                ψ[i] .*= UTt
                pI * ψ[i]
                ψ[i] .*= UVt
                ψ[i] ./= sqrt((ψ[i] .|> abs2 |> sum)*dx)
                E_old = E_new
                pF * ψ[i]
                E_new = dp*kfac^2*sum(T0.*abs2.(ψ[i]))
                pI * ψ[i]
                E_new += dx*sum(V0.*abs2.(ψ[i]))
                update!(prog, desc=@sprintf("Running imag time prop for excited state #%d... E=%.10f ΔE=%.2e", i, E_new, abs(E_new-E_old)), spinner=raw"\-/|")
            end
            finish!(prog)
            println(join([i, E_old, E_new], ", "))
            push!(job.E0, E_new)
        end
    end
    return
end

function purify!(job::TDSEJob_1D_SOF_ImTime, steps=100000)
    dt = job.dt
    dx = job.dx
    ψ = job.ψ
    ψ_pur = job.ψ_pur
    mask = job.mask
    UVt = job.UVt
    UTt = job.UTt
    pF = job.pF
    pI = job.pI
    @. UVt = exp(-0.5im * job.V0 * dt)
    @. UTt = exp(-1im * job.T0 * dt)
    t = 0.0
    for i in 1:job.N_ex
        prog = Progress(steps; desc="Purifying state $i...", dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true)
        push!(ψ_pur, similar(ψ[1]))
        ψ_pur[i] .= 0.0
        for it in 1:steps
            t += dt
            ψ[i] .*= UVt
            pF * ψ[i]
            ψ[i] .*= UTt
            pI * ψ[i]
            ψ[i] .*= UVt
            ψ[i] .*= mask
            ψ_pur[i] .+= ψ[i] * (1-cos(2π/steps*it)) * cis(job.E0[i]*it*dt)
            next!(prog)
        end
        ψ_pur[i] ./= sqrt((ψ_pur[i] .|> abs2 |> sum)*dx)
        finish!(prog)
    end
    return
end
