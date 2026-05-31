# QuantumFisher.jl

**Quantum information geometry on ℂ⁶**

Implements density matrices, the quantum Fisher information tensor 𝓕_AB,
the Bures metric and geodesics, Von Neumann evolution, and the
Kaluza-Klein spectral geometry of K = ℂP² × S³ × S¹.

The mathematical foundation of [Spinoza.jl](https://github.com/uwbanjoman/Spinoza.jl) —
a unified framework deriving the Standard Model, spacetime, gravity,
and consciousness from the single postulate **g_AB = 𝓕_AB / ρ₀**.

---

## Installation

```julia
] add https://github.com/uwbanjoman/QuantumFisher.jl
```

## Quick start

```julia
using QuantumFisher

# The fundamental object: a density matrix on ℂ⁶
ψ = ComplexF64[1, 0, 0, 0, 0, 0]
ρ = density_matrix(ψ)          # pure state |e₁⟩⟨e₁|

# The vacuum: ρ̂* = I/6
ρ_vac = vacuum_state()

# Fisher information
F = fisher_scalar(ρ)           # 5/6 for pure states, 0 for vacuum

# Bures distance
d = bures_distance(ρ, ρ_vac)  # arccos(1/√6) ≈ 65.9° for any pure state

# Von Neumann evolution
H = Matrix(kk_hamiltonian())
ρ_t = evolve_exact(ρ, H, 1.0) # ρ̂(t=1)

# Consciousness measure Φ (Document LXXVII)
Φ = consciousness_measure(ρ)
is_conscious(ρ)                # Φ > τ² = 1/25?
```

---

## The single postulate

```
g_AB = 𝓕_AB / ρ₀
```

The spacetime metric **is** the quantum Fisher information tensor.
All Standard Model structure follows as a geometric consequence.
Zero free parameters.

---

## Core functions

### States

| Function | Description |
|----------|-------------|
| `density_matrix(ψ)` | Pure state ρ̂ = \|ψ⟩⟨ψ\| on ℂ⁶ |
| `vacuum_state()` | Vacuum ρ̂* = I/6 |
| `pure_state(ψ)` | Alias for `density_matrix` |
| `mixed_state(ρ, ε)` | Mix with vacuum: (1-ε)ρ + ε I/6 |
| `gibbs_state(H, β)` | Thermal state exp(-βH)/Z |
| `is_valid_state(ρ)` | Check Hermitian, positive, trace 1 |

### Fisher information

| Function | Description |
|----------|-------------|
| `fisher_tensor(ρ)` | Full 35×35 Fisher tensor 𝓕_AB |
| `fisher_scalar(ρ)` | Tr(ρ²) - 1/6 (scalar measure) |
| `fisher_excess(ρ)` | Excess above vacuum |
| `quantum_fisher_information(ρ, H)` | QFI for parameter H |

### Bures geometry

| Function | Description |
|----------|-------------|
| `bures_fidelity(ρ₁, ρ₂)` | Uhlmann fidelity F(ρ₁,ρ₂) |
| `bures_distance(ρ₁, ρ₂)` | arccos(√F) ∈ [0, π/2] |
| `bures_geodesic(ρ₁, ρ₂, n)` | Geodesic path, n steps |

### Von Neumann evolution

| Function | Description |
|----------|-------------|
| `von_neumann_rhs(ρ, H)` | -i[H, ρ] |
| `evolve_exact(ρ, H, t)` | e^{-iHt} ρ e^{+iHt} |
| `evolve_rk4(ρ, H, t, dt)` | RK4 step |
| `informative_velocity(ρ₁, ρ₂, dt)` | v_info = \|dρ/dt\|_𝓕 |
| `informative_time(traj, times)` | τ_info = ∫√(𝓕 dρ²) dt |

### Entropy

| Function | Description |
|----------|-------------|
| `von_neumann_entropy(ρ)` | -Tr(ρ log ρ) |
| `purity(ρ)` | Tr(ρ²) ∈ [1/6, 1] |
| `relative_entropy(ρ, σ)` | Tr(ρ(log ρ - log σ)) |

### KK spectral geometry

| Function | Description |
|----------|-------------|
| `kk_hamiltonian(n)` | Squared Dirac operator Ð²_K |
| `kk_spectrum(n)` | KK mass spectrum |
| `kk_mass_gap()` | Yang-Mills mass gap Δ = 9/4 |

### Consciousness

| Function | Description |
|----------|-------------|
| `consciousness_measure(ρ)` | Φ = Tr(𝓕_cross)/Tr(𝓕) |
| `is_conscious(ρ)` | Φ > τ² = 1/25? |
| `banach_contraction_factor(ρ)` | L = 1 - 𝓕/𝓕_max |
| `self_model_convergence(ρ, n)` | Banach iteration distance |

---

## Key results

```julia
# All pure states equidistant from vacuum (Document LXXIV)
ρ₁ = pure_state(ComplexF64[1,0,0,0,0,0])
ρ₂ = pure_state(ComplexF64[0,1,0,0,0,0])
bures_distance(ρ₁, vacuum_state()) ≈ bures_distance(ρ₂, vacuum_state())
# → arccos(1/√6) ≈ 1.150 rad ≈ 65.9°

# Yang-Mills mass gap (Document LXXV)
kk_mass_gap()
# → 2.25 = 9/4 > 0  ✓

# Von Neumann is unitary — information conserved
ρ = pure_state(ComplexF64[1,1,1,1,1,1]/sqrt(6))
H = Matrix(kk_hamiltonian())
ρ_t = evolve_exact(ρ, H, 100.0)
purity(ρ_t) ≈ purity(ρ)
# → 1.0  (eigenvalues preserved)  ✓

# Consciousness threshold Φ > τ² = 0.04
consciousness_measure(pure_state(ComplexF64[1,1,0,0,0,0]/sqrt(2)))
# → 0.25 > 0.04  ✓
```

---

## Tests

```julia
] test QuantumFisher
# Results: 62 passed, 0 failed  ✓
```

---

## Theoretical background

QuantumFisher.jl implements the mathematical core of the
FisherGeometrics framework. Key documents:

| Document | Topic |
|----------|-------|
| Preprint | The Standard Model as information geometry |
| Doc LXXIV | Initial conditions — all pure states equidistant |
| Doc LXXV | Yang-Mills mass gap Δ = 9/4 |
| Doc LXXVI | Consciousness as stable Fisher geometry |
| Doc LXXVII | IIT and Fisher geometry (Φ > τ²) |
| Doc LXXXIV | Banach fixed point — stable self-model |
| Doc LXXXV | From consciousness to action — will as geodesic |
| Doc LXXXVI | Derivation of K = ℂP² × S³ × S¹ |

All documents available at:
[github.com/uwbanjoman/FisherGeometrics.jl](https://github.com/uwbanjoman/FisherGeometrics.jl)

---

## Citation

```bibtex
@software{Bouwman2026QuantumFisher,
  author  = {Jan Bouwman},
  title   = {QuantumFisher.jl: Quantum information geometry on ℂ⁶},
  year    = {2026},
  url     = {https://github.com/uwbanjoman/QuantumFisher.jl}
}
```

---

## License

MIT License © 2026 Jan Bouwman

---

*Working document. Speculative theoretical research.*
