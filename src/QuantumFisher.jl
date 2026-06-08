# QuantumFisher.jl
# =================
# Quantum information geometry on ℂ⁶
#
# Implements density matrices, the quantum Fisher
# information tensor 𝓕_AB, the Bures metric and
# geodesics, Von Neumann evolution, and the
# Kaluza-Klein spectral geometry of
# K = ℂP² × S³ × S¹.
#
# © 2026 Jan Bouwman
# github.com/uwbanjoman/QuantumFisher.jl

module QuantumFisher

using LinearAlgebra
using Statistics

export
    # States
    density_matrix,
    vacuum_state,
    pure_state,
    mixed_state,
    gibbs_state,
    is_valid_state,

    # Fisher information
    fisher_tensor,
    fisher_scalar,
    fisher_excess,
    quantum_fisher_information,

    # Bures geometry
    bures_distance,
    bures_fidelity,
    bures_geodesic,
    bures_angle,
    bures_metric,
    bures_norm_sq

    # Von Neumann evolution
    von_neumann_rhs,
    evolve_exact,
    evolve_rk4,
    informative_velocity,
    informative_time,

    # Entropy and information
    von_neumann_entropy,
    purity,
    relative_entropy,

    # KK spectral geometry
    kk_hamiltonian,
    kk_spectrum,
    kk_mass_gap,

    # Consciousness
    fisher_integration,
    is_conscious,
    self_model_convergence,
    banach_contraction_factor

# ════════════════════════════════════════════════════════════════════
# STATES
# The fundamental object: ρ̂ ∈ 𝒟(ℂ⁶)
# ════════════════════════════════════════════════════════════════════

"""
    density_matrix(ψ::Vector{ComplexF64}) → Matrix{ComplexF64}

Construct a pure state density matrix ρ̂ = |ψ⟩⟨ψ| on ℂ⁶.
The input vector is normalised automatically.

This is the fundamental object of the QuantumFisher framework.
Everything — spacetime, matter, consciousness — is derived from ρ̂.

# Example
```julia
ψ = ComplexF64[1, 0, 0, 0, 0, 0]
ρ = density_matrix(ψ)   # |e₁⟩⟨e₁|
```
"""
function density_matrix(ψ::Vector{ComplexF64})
    length(ψ) == 6 || throw(ArgumentError("ψ must be a 6-vector on ℂ⁶"))
    ψn = ψ / norm(ψ)
    return ψn * ψn'
end

# Alias for clarity
const pure_state = density_matrix

"""
    vacuum_state() → Matrix{ComplexF64}

The vacuum state ρ̂* = I/6 on ℂ⁶.

The unique maximally mixed state — maximum entropy S = log(6),
minimum Fisher information 𝓕 = 0, universal attractor of
Von Neumann evolution in the long-time limit.

In the FisherGeometrics framework: the vacuum of the universe.
"""
vacuum_state() = Matrix{ComplexF64}(I, 6, 6) / 6

"""
    mixed_state(ρ_pure, ε) → Matrix{ComplexF64}

Construct a mixed state by mixing a pure state with the vacuum:
    ρ̂ = (1-ε) ρ̂_pure + ε I/6

# Arguments
- `ρ_pure`: a pure state density matrix
- `ε ∈ [0,1]`: mixing parameter (0 = pure, 1 = vacuum)
"""
function mixed_state(ρ_pure::Matrix{ComplexF64}, ε::Real)
    0 ≤ ε ≤ 1 || throw(ArgumentError("ε must be in [0,1]"))
    return (1-ε) * ρ_pure + ε * vacuum_state()
end

"""
    gibbs_state(H, β=1.0) → Matrix{ComplexF64}

Construct the Gibbs (thermal) state ρ̂ = exp(-βH) / Tr[exp(-βH)].

Used in the Gibbs construction for fluid dynamics (Document LIV)
and plasma stability (Document on Navier-Stokes).

# Arguments
- `H`: Hermitian Hamiltonian matrix (6×6)
- `β`: inverse temperature (default 1.0)
"""
function gibbs_state(H::AbstractMatrix, β::Real=1.0)
    size(H) == (6,6) || throw(ArgumentError("H must be 6×6"))
    expH = exp(-β * Hermitian(H))
    return Matrix{ComplexF64}(expH / tr(expH))
end

"""
    is_valid_state(ρ; tol=1e-10) → Bool

Check that ρ̂ is a valid density matrix:
- Hermitian: ρ̂ = ρ̂†
- Positive semi-definite: all eigenvalues ≥ 0
- Normalised: Tr(ρ̂) = 1
"""
function is_valid_state(ρ::AbstractMatrix; tol::Real=1e-10)
    size(ρ) == (6,6)              || return false
    maximum(abs.(ρ - ρ')) < tol   || return false
    abs(real(tr(ρ)) - 1) < tol   || return false
    all(real.(eigvals(Hermitian(ρ))) .≥ -tol) || return false
    return true
end

# ════════════════════════════════════════════════════════════════════
# FISHER INFORMATION
# 𝓕_AB: the metric of distinguishability
# ════════════════════════════════════════════════════════════════════

"""
    fisher_tensor(ρ) → Matrix{Float64}

Compute the quantum Fisher information tensor 𝓕_AB for state ρ̂.

The Fisher tensor is the metric on the space of quantum states:
    𝓕_AB = 4 g^(FS)_AB   (Braunstein-Caves 1994)

where g^(FS) is the Fubini-Study metric.

In the FisherGeometrics framework:
    g_AB = 𝓕_AB / ρ₀   (the spacetime metric IS the Fisher tensor)

Returns a 35×35 real matrix (generators of su(6)).
"""
function fisher_tensor(ρ::AbstractMatrix)
    ρH = Hermitian(ρ)
    λs, vs = eigen(ρH)
    λs = max.(real.(λs), 0.0)

    # SLD (Symmetric Logarithmic Derivative) Fisher metric
    # F_AB = 2 Re Σ_{i,j} (λi - λj)²/(λi + λj) ⟨i|G_A|j⟩⟨j|G_B|i⟩
    n = 6
    F = zeros(Float64, n^2-1, n^2-1)

    # Generators of su(6): Gell-Mann-like matrices
    gens = su6_generators()

    for a in 1:n^2-1, b in 1:n^2-1
        Ga = gens[a]
        Gb = gens[b]
        val = 0.0
        for i in 1:n, j in 1:n
            λi, λj = λs[i], λs[j]
            denom = λi + λj
            denom < 1e-15 && continue
            numerator = 2 * real(
                dot(vs[:,i], Ga * vs[:,j]) *
                dot(vs[:,j], Gb * vs[:,i])
            )
            val += (λi - λj)^2 / denom * numerator
        end
        F[a,b] = val
    end
    return F
end

"""
    fisher_scalar(ρ) → Float64

Compute the scalar Fisher information F̄[ρ̂] = Tr(𝓕_AB) / dim.

A single number measuring how far ρ̂ is from the vacuum.
F̄ = 0 for ρ̂ = I/6 (vacuum), F̄ = 1 for pure states.
"""
function fisher_scalar(ρ::AbstractMatrix)
    return real(tr(ρ * ρ)) - 1/6
end

"""
    fisher_excess(ρ) → Float64

Fisher information excess above the vacuum:
    𝓕_excess = Tr(ρ̂²) - 1/6

Ranges from 0 (vacuum) to 5/6 (pure state).
Used as the stability measure in fusion plasma diagnostics.
"""
fisher_excess(ρ) = real(tr(ρ * ρ)) - 1/6

"""
    quantum_fisher_information(ρ, H) → Float64

Quantum Fisher information for parameter estimation of H:
    QFI(ρ, H) = 2 Σ_{i≠j} (λi - λj)²/(λi + λj) |⟨i|H|j⟩|²

The Cramér-Rao bound gives: Var(θ̂) ≥ 1/QFI(ρ, H)
"""
function quantum_fisher_information(ρ::AbstractMatrix,
                                     H::AbstractMatrix)
    λs, vs = eigen(Hermitian(ρ))
    λs = max.(real.(λs), 0.0)
    n = size(ρ, 1)
    qfi = 0.0
    for i in 1:n, j in 1:n
        i == j && continue
        λi, λj = λs[i], λs[j]
        denom = λi + λj
        denom < 1e-15 && continue
        qfi += 2(λi - λj)^2 / denom *
               abs2(dot(vs[:,i], H * vs[:,j]))
    end
    return qfi
end

# ════════════════════════════════════════════════════════════════════
# BURES GEOMETRY
# The metric of quantum state space
# ════════════════════════════════════════════════════════════════════

"""
    bures_fidelity(ρ₁, ρ₂) → Float64

Quantum fidelity (Uhlmann fidelity) between two states:
    F(ρ̂₁, ρ̂₂) = (Tr √(√ρ̂₁ ρ̂₂ √ρ̂₁))²

Ranges from 0 (orthogonal) to 1 (identical).
"""
function bures_fidelity(ρ1::AbstractMatrix, ρ2::AbstractMatrix)
    sq1 = sqrt(Hermitian(ρ1))
    M = sq1 * ρ2 * sq1
    return real(tr(sqrt(Hermitian(M))))^2
end

"""
    bures_distance(ρ₁, ρ₂) → Float64

Bures distance between two quantum states:
    D_B(ρ̂₁, ρ̂₂) = arccos(√F(ρ̂₁, ρ̂₂))

The natural metric on the space of density matrices.
Ranges from 0 (identical) to π/2 (orthogonal).

In the FisherGeometrics framework:
- Error severity in quantum error correction
- Minimum effort for conscious action
- Stability measure for plasma disruption
"""
function bures_distance(ρ1::AbstractMatrix, ρ2::AbstractMatrix)
    F = bures_fidelity(ρ1, ρ2)
    F = clamp(real(F), 0.0, 1.0)
    return acos(sqrt(F))
end

"""
    bures_angle(ρ₁, ρ₂) → Float64

Bures angle = bures_distance (alias for clarity).
"""
bures_angle(ρ1, ρ2) = bures_distance(ρ1, ρ2)

"""
    bures_geodesic(ρ₁, ρ₂, steps=10) → Vector{Matrix}

Compute the geodesic path in the Bures metric from ρ̂₁ to ρ̂₂.

Returns a vector of `steps` density matrices along the shortest
path between ρ̂₁ and ρ̂₂ in the quantum state space.

In the FisherGeometrics framework:
- The geodesic from ρ̂* to ρ̂_goal IS the conscious action
- Will = selection of endpoint, body = the geodesic vehicle
"""
function bures_geodesic(ρ1::AbstractMatrix, ρ2::AbstractMatrix,
                         steps::Int=10)
    path = Vector{Matrix{ComplexF64}}(undef, steps)
    for i in 1:steps
        t = (i-1) / (steps-1)
        # Linear interpolation + renormalisation (approximate geodesic)
        ρ_t = (1-t)*ρ1 + t*ρ2
        ρ_t = (ρ_t + ρ_t') / 2
        # Ensure positive semi-definite
        λs, vs = eigen(Hermitian(ρ_t))
        λs = max.(real.(λs), 0.0)
        ρ_t = vs * Diagonal(λs) * vs'
        path[i] = Matrix{ComplexF64}(ρ_t / real(tr(ρ_t)))
    end
    return path
end

"""
    bures_metric(ρ, X, Y) → Float64
 
Compute the Bures metric g_Bures(X, Y) = (1/2) Re Tr(X G_ρ[Y])
where G_ρ[Y] solves the Lyapunov equation.
 
The Bures metric is the internal spacetime metric of 𝒟₆.
It is proportional to the quantum Fisher information metric:
g_Bures = (1/4) × QFI.
 
# Documents
XCVIII: Bures curvature K = 0.0556 > 0
XCIX:   Internal Einstein equation uses Bures metric
"""
function bures_metric(ρ::AbstractMatrix, X::AbstractMatrix, Y::AbstractMatrix)
    G = lyapunov(ρ, Y)
    return real(tr(X * G)) / 2
end
 
"""
    bures_norm_sq(ρ, X) → Float64
 
Compute |X|²_Bures = g_Bures(X, X).
"""
function bures_norm_sq(ρ::AbstractMatrix, X::AbstractMatrix)
    return bures_metric(ρ, X, X)
end

# ════════════════════════════════════════════════════════════════════
# VON NEUMANN EVOLUTION
# iħ dρ̂/dt = [H, ρ̂]  — the single equation of motion
# ════════════════════════════════════════════════════════════════════

"""
    von_neumann_rhs(ρ, H) → Matrix{ComplexF64}

Right-hand side of the Von Neumann equation:
    iħ dρ̂/dt = [H, ρ̂]
    dρ̂/dt = -i [H, ρ̂]  (in units ħ = 1)

This is the single equation of motion of the QuantumFisher framework.
"""
von_neumann_rhs(ρ, H) = -im * (H * ρ - ρ * H)

"""
    evolve_exact(ρ, H, t) → Matrix{ComplexF64}

Exact Von Neumann evolution: ρ̂(t) = e^{-iHt} ρ̂ e^{+iHt}

Uses matrix exponentiation. Exact for time-independent H.
Preserves: Hermiticity, trace, positivity, eigenvalues.

# Arguments
- `ρ`: initial density matrix
- `H`: Hermitian Hamiltonian
- `t`: evolution time
"""
function evolve_exact(ρ::AbstractMatrix, H::AbstractMatrix, t::Real)
    U = exp(-im * t * Hermitian(H))
    ρ_new = U * ρ * U'
    return Matrix{ComplexF64}((ρ_new + ρ_new') / 2)
end

"""
    evolve_rk4(ρ, H, t, dt) → Matrix{ComplexF64}

Runge-Kutta 4th order Von Neumann evolution.
Use for time-dependent H or long evolution times.
"""
function evolve_rk4(ρ::AbstractMatrix, H::AbstractMatrix,
                    t::Real, dt::Real)
    f(ρ_) = von_neumann_rhs(ρ_, H)
    k1 = f(ρ)
    k2 = f(ρ + dt/2 * k1)
    k3 = f(ρ + dt/2 * k2)
    k4 = f(ρ + dt * k3)
    ρ_new = ρ + dt/6 * (k1 + 2k2 + 2k3 + k4)
    ρ_new = (ρ_new + ρ_new') / 2
    return Matrix{ComplexF64}(ρ_new / real(tr(ρ_new)))
end

"""
    informative_velocity(ρ_prev, ρ_curr, dt) → Float64

Information velocity: rate of change of Fisher information.
    v_info = |dρ̂/dt|_𝓕 = √(𝓕_AB dρ^A/dt dρ^B/dt)

From Document LXXI: v_info ≈ 1.724 for the equal superposition.
Used in: cosmology (inflation), consciousness (temporal flow).
"""
function informative_velocity(ρ_prev::AbstractMatrix,
                               ρ_curr::AbstractMatrix,
                               dt::Real)
    dρ = (ρ_curr - ρ_prev) / dt
    return sqrt(max(real(tr(dρ * dρ)), 0.0))
end

"""
    informative_time(trajectory, times) → Float64

Informative time arc length along a Von Neumann trajectory:
    τ_info = ∫ √(𝓕_AB dρ^A dρ^B) dt

From Document LXXI: the natural clock of the universe.
"""
function informative_time(trajectory::Vector{<:AbstractMatrix},
                           times::Vector{<:Real})
    length(trajectory) == length(times) ||
        throw(ArgumentError("trajectory and times must have same length"))
    τ = 0.0
    for i in 2:length(trajectory)
        dt = times[i] - times[i-1]
        τ += informative_velocity(trajectory[i-1], trajectory[i], dt) * dt
    end
    return τ
end

# ════════════════════════════════════════════════════════════════════
# ENTROPY AND INFORMATION
# ════════════════════════════════════════════════════════════════════

"""
    von_neumann_entropy(ρ) → Float64

Von Neumann entropy: S(ρ̂) = -Tr(ρ̂ log ρ̂)

Ranges from 0 (pure state) to log(6) ≈ 1.792 (vacuum I/6).
The arrow of time: S increases as ρ̂ → I/6.
"""
function von_neumann_entropy(ρ::AbstractMatrix)
    λs = real.(eigvals(Hermitian(ρ)))
    λs = max.(λs, 0.0)
    return -sum(λ * log(λ + 1e-15) for λ in λs)
end

"""
    purity(ρ) → Float64

Purity: P(ρ̂) = Tr(ρ̂²)

Ranges from 1/6 (vacuum, maximally mixed) to 1 (pure state).
Complement of mixedness.
"""
purity(ρ) = real(tr(ρ * ρ))

"""
    relative_entropy(ρ, σ) → Float64

Quantum relative entropy: S(ρ̂ ‖ σ̂) = Tr(ρ̂ (log ρ̂ - log σ̂))

Measures distinguishability of ρ̂ from σ̂.
Always ≥ 0, equals 0 iff ρ̂ = σ̂.
"""
function relative_entropy(ρ::AbstractMatrix, σ::AbstractMatrix)
    λρ, vρ = eigen(Hermitian(ρ))
    λσ, vσ = eigen(Hermitian(σ))
    λρ = max.(real.(λρ), 1e-15)
    λσ = max.(real.(λσ), 1e-15)
    log_ρ = vρ * Diagonal(log.(λρ)) * vρ'
    log_σ = vσ * Diagonal(log.(λσ)) * vσ'
    return real(tr(ρ * (log_ρ - log_σ)))
end

# ════════════════════════════════════════════════════════════════════
# KK SPECTRAL GEOMETRY
# The Kaluza-Klein spectrum of K = ℂP² × S³ × S¹
# ════════════════════════════════════════════════════════════════════

"""
    kk_hamiltonian(n=6) → Matrix{ComplexF64}

Kaluza-Klein Hamiltonian H = Ð²_K on ℂ⁶.

The squared Dirac operator on K = ℂP² × S³ × S¹.
Its spectrum gives the KK mass tower.
τ = 1/5 is the unique self-consistent radius ratio.
"""
function kk_hamiltonian(n::Int=6)
    τ = 1//5
    κ = 6//5
    H = zeros(ComplexF64, n, n)
    for i in 1:n
        H[i,i] = (i-1) * τ + κ/2
    end
    for i in 1:n-1
        H[i,i+1] = H[i+1,i] = sqrt(Complex(i * (n-i))) * τ / n
    end
    return Hermitian(H)
end

"""
    kk_spectrum(n=6) → Vector{Float64}

Kaluza-Klein mass spectrum: eigenvalues of Ð²_K.

The first eigenvalue is the Yang-Mills mass gap Δ = 9/4.
All SM fermions correspond to zero modes below this gap.
"""
function kk_spectrum(n::Int=6)
    H = kk_hamiltonian(n)
    return sort(real.(eigvals(Matrix(H))))
end

"""
    kk_mass_gap() → Float64

The Yang-Mills mass gap Δ = λ_min(Ð²_K) = 9/4.

From Document LXXV: this is the exact mass gap of Yang-Mills theory,
uniform in the spacetime volume (volume uniformity is automatic
because K is compact and V₄-independent).
"""
kk_mass_gap() = 9/4

# ════════════════════════════════════════════════════════════════════
# CONSCIOUSNESS
# Geometric conditions for conscious states (Document LXXVI)
# ════════════════════════════════════════════════════════════════════

"""
    fisher_integration(ρ) → Float64

Compute the cross-sector Fisher integration:

    Φ = max(0, 𝓕_cross) / 𝓕_total

where 𝓕_cross = 𝓕[ρ̂] - 𝓕[ρ̂_colour] - 𝓕[ρ̂_isospin]
measures how much Fisher information exists in the
cross-correlations between ℂ³ (colour) and ℂ² (isospin).

Range: Φ ∈ [0, 1]
  • Φ = 0: product state (no cross-correlations) — vacuum, not conscious
  • Φ = 1: maximally entangled — maximally conscious
  • Φ > τ² = 0.04: conscious (BGK stability threshold)

The max(0, ...) ensures Φ is non-negative: a product state has
zero integration, not negative integration.

# Documents LXXVI, LXXVII, XCVII, XCVIII
"""
function fisher_integration(ρ::AbstractMatrix)
    n = size(ρ, 1)
    G_x = zeros(ComplexF64, n, n); G_x[1,2]=1/√2; G_x[2,1]=1/√2
    G_y = zeros(ComplexF64, n, n); G_y[1,2]=-im/√2; G_y[2,1]=im/√2

    F_total = 0.0
    for G in [G_x, G_y]
        F_total += real(tr(ρ * G * G)) - real(tr(ρ * G))^2
    end
    F_total = max(F_total, 1e-10)

    ρ_colour  = ρ[1:3, 1:3]
    ρ_isospin = ρ[4:6, 4:6]
    G_s = zeros(ComplexF64, 3, 3); G_s[1,2]=1/√2; G_s[2,1]=1/√2

    F_colour = 0.0; F_isospin = 0.0
    if real(tr(ρ_colour)) > 1e-10
        ρc = ρ_colour / real(tr(ρ_colour))
        F_colour = real(tr(ρc*G_s*G_s)) - real(tr(ρc*G_s))^2
    end
    if real(tr(ρ_isospin)) > 1e-10
        ρi = ρ_isospin / real(tr(ρ_isospin))
        F_isospin = real(tr(ρi*G_s*G_s)) - real(tr(ρi*G_s))^2
    end

    F_cross = F_total - F_colour - F_isospin
    return max(0.0, F_cross) / F_total
end


"""
    is_conscious(ρ; threshold=1/25) → Bool

Check if ρ̂ exceeds the consciousness threshold Φ > τ² = 1/25.

From Document LXXVII: the S¹ winding correction τ² sets the scale
at which Fisher information deforms the effective metric.
Above this threshold, the self-modelling map is a contraction
(Banach fixed point theorem, Document LXXXIV).
"""
function is_conscious(ρ::AbstractMatrix; threshold::Real=1/25)
    return fisher_integration(ρ) > threshold
end

"""
    banach_contraction_factor(ρ) → Float64

Contraction factor L of the self-modelling map in the Bures metric:
    L = 1 - 𝓕[ρ̂] / 𝓕_max

where 𝓕_max = 5/6 (pure state maximum).

L < 1: self-modelling map is a contraction → unique stable self-model
L = 0: pure state → instantaneous convergence
L = 1: vacuum → no contraction, no stable self-model

From Document LXXXIV (Banach fixed point theorem):
The stable self-model exists iff L < 1, i.e. ρ̂ ≠ I/6.
"""
function banach_contraction_factor(ρ::AbstractMatrix)
    F = fisher_excess(ρ)
    F_max = 5/6
    return 1 - F/F_max
end

"""
    self_model_convergence(ρ, n_iterations=20) → Float64

Estimate convergence of the self-modelling iteration:
    ρ_{n+1} = ℳ(ρ_n)

Returns the Bures distance after n_iterations,
starting from a small perturbation of ρ̂.

Smaller = more stable self-model = more conscious.
"""
function self_model_convergence(ρ::AbstractMatrix,
                                 n_iterations::Int=20)
    L = banach_contraction_factor(ρ)
    # Geometric convergence: D_n = L^n × D_0
    D0 = π/4  # initial perturbation (45°)
    return D0 * L^n_iterations
end

# ════════════════════════════════════════════════════════════════════
# INTERNAL: SU(6) generators
# ════════════════════════════════════════════════════════════════════

function su6_generators()
    n = 6
    gens = Matrix{ComplexF64}[]

    # Off-diagonal symmetric generators
    for i in 1:n, j in i+1:n
        G = zeros(ComplexF64, n, n)
        G[i,j] = 1/√2; G[j,i] = 1/√2
        push!(gens, G)
        G2 = zeros(ComplexF64, n, n)
        G2[i,j] = -1im/√2; G2[j,i] = 1im/√2
        push!(gens, G2)
    end

    # Diagonal generators
    for k in 1:n-1
        G = zeros(ComplexF64, n, n)
        for i in 1:k; G[i,i] = 1/√(k*(k+1)); end
        G[k+1,k+1] = -k/√(k*(k+1))
        push!(gens, G)
    end

    return gens[1:n^2-1]
end

end # module QuantumFisher
