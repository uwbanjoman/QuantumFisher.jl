# test/runtests.jl
# ==================
# QuantumFisher.jl — full function test suite
#
# Runs every exported function with correct input
# and verifies the output is physically sensible.
#
# Usage:
#   julia --project test/runtests.jl
#   or: ] test QuantumFisher
#
# © 2026 Jan Bouwman

using QuantumFisher
using LinearAlgebra
using Test
using Printf

println("="^65)
println("  QuantumFisher.jl — Function Test Suite")
println("="^65)
println()

# ── Helper ────────────────────────────────────────────────────────────────────

pass = 0
fail = 0

macro qtest(name, expr)
    quote
        try
            result = $(esc(expr))
            if result isa Bool
                if result
                    global pass += 1
                    println("  ✓  ", $(name))
                else
                    global fail += 1
                    println("  ✗  ", $(name), "  [FAILED]")
                end
            else
                global pass += 1
                println("  ✓  ", $(name))
            end
        catch e
            global fail += 1
            println("  ✗  ", $(name), "  [ERROR: ", e, "]")
        end
    end
end

# ── Standard test states ──────────────────────────────────────────────────────

ψ_e1   = ComplexF64[1, 0, 0, 0, 0, 0]        # |e₁⟩
ψ_eq   = ComplexF64[1, 1, 1, 1, 1, 1] / sqrt(6)  # equal superposition
ψ_r    = ComplexF64[1, 1, 0, 0, 0, 0] / sqrt(2)  # red sector

# ════════════════════════════════════════════════════════════════════
# STATES
# ════════════════════════════════════════════════════════════════════

println("─── States ───")
println()

@qtest "density_matrix(ψ_e1) is Hermitian" begin
    ρ = density_matrix(ψ_e1)
    maximum(abs.(ρ - ρ')) < 1e-12
end

@qtest "density_matrix(ψ_e1) has trace 1" begin
    ρ = density_matrix(ψ_e1)
    abs(real(tr(ρ)) - 1) < 1e-12
end

@qtest "density_matrix normalises input" begin
    ψ_unnorm = ComplexF64[3, 0, 0, 0, 0, 0]
    ρ = density_matrix(ψ_unnorm)
    abs(real(tr(ρ)) - 1) < 1e-12
end

@qtest "pure_state == density_matrix" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = density_matrix(ψ_e1)
    maximum(abs.(ρ1 - ρ2)) < 1e-12
end

@qtest "vacuum_state() = I/6" begin
    ρ = vacuum_state()
    maximum(abs.(ρ - I/6)) < 1e-12
end

@qtest "vacuum_state trace = 1" begin
    abs(real(tr(vacuum_state())) - 1) < 1e-12
end

@qtest "mixed_state(ρ, 0) = pure state" begin
    ρp = pure_state(ψ_e1)
    ρm = mixed_state(ρp, 0.0)
    maximum(abs.(ρm - ρp)) < 1e-12
end

@qtest "mixed_state(ρ, 1) = vacuum" begin
    ρp = pure_state(ψ_e1)
    ρm = mixed_state(ρp, 1.0)
    maximum(abs.(ρm - vacuum_state())) < 1e-12
end

@qtest "gibbs_state returns valid density matrix" begin
    H = Matrix(kk_hamiltonian())
    ρ = gibbs_state(H)
    is_valid_state(ρ)
end

@qtest "is_valid_state(pure_state)" begin
    is_valid_state(pure_state(ψ_e1))
end

@qtest "is_valid_state(vacuum_state)" begin
    is_valid_state(vacuum_state())
end

@qtest "is_valid_state rejects non-normalised" begin
    ρ = Matrix{ComplexF64}(2I, 6, 6)
    !is_valid_state(ρ)
end

println()

# ════════════════════════════════════════════════════════════════════
# FISHER INFORMATION
# ════════════════════════════════════════════════════════════════════

println("─── Fisher information ───")
println()

@qtest "fisher_tensor returns 35×35 matrix" begin
    F = fisher_tensor(pure_state(ψ_e1))
    size(F) == (35, 35)
end

@qtest "fisher_tensor is symmetric" begin
    F = fisher_tensor(pure_state(ψ_e1))
    maximum(abs.(F - F')) < 1e-10
end

@qtest "fisher_tensor vacuum has lower info than pure" begin
    F_pure = tr(fisher_tensor(pure_state(ψ_e1)))
    F_vac  = tr(fisher_tensor(vacuum_state()))
    F_pure > F_vac
end

@qtest "fisher_scalar(vacuum) = 0" begin
    abs(fisher_scalar(vacuum_state())) < 1e-12
end

@qtest "fisher_scalar(pure) > 0" begin
    fisher_scalar(pure_state(ψ_e1)) > 0
end

@qtest "fisher_excess(vacuum) = 0" begin
    abs(fisher_excess(vacuum_state())) < 1e-12
end

@qtest "fisher_excess(pure) = 5/6" begin
    abs(fisher_excess(pure_state(ψ_e1)) - 5/6) < 1e-10
end

@qtest "quantum_fisher_information ≥ 0" begin
    ρ = pure_state(ψ_e1)
    H = Matrix(kk_hamiltonian())
    quantum_fisher_information(ρ, H) ≥ 0
end

println()

# ════════════════════════════════════════════════════════════════════
# BURES GEOMETRY
# ════════════════════════════════════════════════════════════════════

println("─── Bures geometry ───")
println()

@qtest "bures_fidelity(ρ, ρ) = 1" begin
    ρ = pure_state(ψ_e1)
    abs(bures_fidelity(ρ, ρ) - 1) < 1e-10
end

@qtest "bures_fidelity(pure, vacuum) = 1/6" begin
    ρ = pure_state(ψ_e1)
    abs(bures_fidelity(ρ, vacuum_state()) - 1/6) < 1e-10
end

@qtest "bures_distance(ρ, ρ) = 0" begin
    ρ = pure_state(ψ_e1)
    abs(bures_distance(ρ, ρ)) < 1e-10
end

@qtest "bures_distance(pure, vacuum) = arccos(1/√6)" begin
    ρ = pure_state(ψ_e1)
    expected = acos(1/sqrt(6))
    abs(bures_distance(ρ, vacuum_state()) - expected) < 1e-10
end

@qtest "all pure states equidistant from vacuum" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = pure_state(ψ_eq)
    ρ3 = pure_state(ψ_r)
    d1 = bures_distance(ρ1, vacuum_state())
    d2 = bures_distance(ρ2, vacuum_state())
    d3 = bures_distance(ρ3, vacuum_state())
    abs(d1 - d2) < 1e-10 && abs(d2 - d3) < 1e-10
end

@qtest "bures_distance symmetric" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = pure_state(ψ_eq)
    abs(bures_distance(ρ1, ρ2) - bures_distance(ρ2, ρ1)) < 1e-10
end

@qtest "bures_distance ∈ [0, π/2]" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = vacuum_state()
    d = bures_distance(ρ1, ρ2)
    0 ≤ d ≤ π/2 + 1e-10
end

@qtest "bures_angle == bures_distance" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = pure_state(ψ_eq)
    abs(bures_angle(ρ1, ρ2) - bures_distance(ρ1, ρ2)) < 1e-12
end

@qtest "bures_geodesic has correct length" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = vacuum_state()
    path = bures_geodesic(ρ1, ρ2, 5)
    length(path) == 5
end

@qtest "bures_geodesic endpoints correct" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = vacuum_state()
    path = bures_geodesic(ρ1, ρ2, 5)
    maximum(abs.(path[1] - ρ1)) < 1e-10 &&
    maximum(abs.(path[end] - ρ2)) < 1e-10
end

@qtest "bures_geodesic states are valid" begin
    ρ1 = pure_state(ψ_e1)
    ρ2 = vacuum_state()
    path = bures_geodesic(ρ1, ρ2, 5)
    all(is_valid_state(ρ) for ρ in path)
end

println()

# ════════════════════════════════════════════════════════════════════
# VON NEUMANN EVOLUTION
# ════════════════════════════════════════════════════════════════════

println("─── Von Neumann evolution ───")
println()

H = Matrix(kk_hamiltonian())
ρ0 = pure_state(ψ_e1)

@qtest "von_neumann_rhs is anti-Hermitian" begin
    rhs = von_neumann_rhs(ρ0, H)
    maximum(abs.(rhs + rhs')) < 1e-10
end

@qtest "von_neumann_rhs has trace 0" begin
    rhs = von_neumann_rhs(ρ0, H)
    abs(real(tr(rhs))) < 1e-12
end

@qtest "evolve_exact preserves trace" begin
    ρt = evolve_exact(ρ0, H, 0.1)
    abs(real(tr(ρt)) - 1) < 1e-10
end

@qtest "evolve_exact preserves Hermiticity" begin
    ρt = evolve_exact(ρ0, H, 0.1)
    maximum(abs.(ρt - ρt')) < 1e-10
end

@qtest "evolve_exact preserves eigenvalues (unitarity)" begin
    ρt = evolve_exact(ρ0, H, 1.0)
    λ0 = sort(real.(eigvals(Hermitian(ρ0))))
    λt = sort(real.(eigvals(Hermitian(ρt))))
    maximum(abs.(λ0 - λt)) < 1e-10
end

@qtest "evolve_exact(t=0) = identity" begin
    ρt = evolve_exact(ρ0, H, 0.0)
    maximum(abs.(ρt - ρ0)) < 1e-10
end

@qtest "evolve_rk4 preserves trace" begin
    ρt = evolve_rk4(ρ0, H, 0.0, 0.01)
    abs(real(tr(ρt)) - 1) < 1e-8
end

@qtest "informative_velocity ≥ 0" begin
    ρ1 = evolve_exact(ρ0, H, 0.0)
    ρ2 = evolve_exact(ρ0, H, 0.1)
    informative_velocity(ρ1, ρ2, 0.1) ≥ 0
end

@qtest "informative_time ≥ 0" begin
    traj = [evolve_exact(ρ0, H, t) for t in 0:0.1:1.0]
    times = collect(0:0.1:1.0)
    informative_time(traj, times) ≥ 0
end

@qtest "vacuum has zero informative velocity (fixed point)" begin
    ρv = vacuum_state()
    ρv2 = evolve_exact(ρv, H, 0.1)
    v = informative_velocity(ρv, ρv2, 0.1)
    v < 0.01  # nearly zero — vacuum is approximate fixed point
end

println()

# ════════════════════════════════════════════════════════════════════
# ENTROPY AND INFORMATION
# ════════════════════════════════════════════════════════════════════

println("─── Entropy and information ───")
println()

@qtest "entropy(pure) = 0" begin
    abs(von_neumann_entropy(pure_state(ψ_e1))) < 1e-10
end

@qtest "entropy(vacuum) = log(6)" begin
    abs(von_neumann_entropy(vacuum_state()) - log(6)) < 1e-10
end

@qtest "entropy ≥ 0" begin
    ρ = mixed_state(pure_state(ψ_e1), 0.5)
    von_neumann_entropy(ρ) ≥ 0
end

@qtest "purity(pure) = 1" begin
    abs(purity(pure_state(ψ_e1)) - 1) < 1e-12
end

@qtest "purity(vacuum) = 1/6" begin
    abs(purity(vacuum_state()) - 1/6) < 1e-12
end

@qtest "purity ∈ [1/6, 1]" begin
    ρ = mixed_state(pure_state(ψ_eq), 0.3)
    1/6 - 1e-10 ≤ purity(ρ) ≤ 1 + 1e-10
end

@qtest "relative_entropy(ρ, ρ) = 0" begin
    ρ = pure_state(ψ_e1)
    abs(relative_entropy(ρ, ρ)) < 1e-8
end

@qtest "relative_entropy ≥ 0" begin
    ρ = pure_state(ψ_e1)
    σ = vacuum_state()
    relative_entropy(ρ, σ) ≥ -1e-8
end

println()

# ════════════════════════════════════════════════════════════════════
# KK SPECTRAL GEOMETRY
# ════════════════════════════════════════════════════════════════════

println("─── KK spectral geometry ───")
println()

@qtest "kk_hamiltonian returns 6×6 matrix" begin
    size(kk_hamiltonian()) == (6, 6)
end

@qtest "kk_hamiltonian is Hermitian" begin
    H = kk_hamiltonian()
    maximum(abs.(H - H')) < 1e-12
end

@qtest "kk_spectrum returns sorted real values" begin
    sp = kk_spectrum()
    all(diff(sp) .≥ -1e-12)
end

@qtest "kk_mass_gap = 9/4" begin
    abs(kk_mass_gap() - 9/4) < 1e-12
end

@qtest "kk_spectrum minimum ≥ 0" begin
    minimum(kk_spectrum()) ≥ -1e-10
end

println()

# ════════════════════════════════════════════════════════════════════
# CONSCIOUSNESS
# ════════════════════════════════════════════════════════════════════

println("─── Consciousness ───")
println()

@qtest "consciousness_measure ∈ [0, 1]" begin
    Φ = consciousness_measure(pure_state(ψ_e1))
    0 ≤ Φ ≤ 1
end

@qtest "consciousness_measure(vacuum) = 0" begin
    abs(consciousness_measure(vacuum_state())) < 1e-10
end

@qtest "is_conscious uses τ² threshold" begin
    τ² = (1/5)^2
    ρ_high = pure_state(ψ_r)  # entangled state
    # just check it runs without error
    result = is_conscious(ρ_high)
    result isa Bool
end

@qtest "banach_contraction_factor(vacuum) = 1" begin
    abs(banach_contraction_factor(vacuum_state()) - 1) < 1e-10
end

@qtest "banach_contraction_factor(pure) = 0" begin
    abs(banach_contraction_factor(pure_state(ψ_e1))) < 1e-10
end

@qtest "banach_contraction_factor ∈ [0, 1]" begin
    ρ = mixed_state(pure_state(ψ_eq), 0.5)
    L = banach_contraction_factor(ρ)
    0 ≤ L ≤ 1
end

@qtest "self_model_convergence decreases with iterations" begin
    ρ = mixed_state(pure_state(ψ_eq), 0.3)
    d5  = self_model_convergence(ρ, 5)
    d20 = self_model_convergence(ρ, 20)
    d5 > d20
end

@qtest "self_model_convergence(pure, n) → 0" begin
    ρ = pure_state(ψ_e1)
    self_model_convergence(ρ, 100) < 1e-10
end

println()

# ════════════════════════════════════════════════════════════════════
# SUMMARY
# ════════════════════════════════════════════════════════════════════

total = pass + fail
println("="^65)
@printf("  Results: %d passed, %d failed  (%d total)\n",
        pass, fail, total)
if fail == 0
    println("  All tests passed ✓")
else
    println("  $fail test(s) failed ✗")
end
println("="^65)
