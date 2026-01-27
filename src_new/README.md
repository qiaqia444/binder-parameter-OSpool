# src_new: Restructured Quantum Simulation Library

This directory contains a refactored and modularized version of the quantum many-body simulation code, organized following best practices inspired by the `src_1` architecture.

## Architecture Overview

The codebase is organized into logical modules, each with a single responsibility:

```
src_new/
├── types.jl                   # State type definitions (PureStateMPS, DiagonalStateMPS, MixedStateMPS)
├── state_constructors.jl      # Initial state preparation (type-safe constructors)
├── measurements.jl            # Weak measurement operators (basic, for backward compatibility)
├── measurements_typed.jl      # Type-safe measurements with multiple dispatch
├── dynamics.jl                # Time evolution and quantum channels
├── binder.jl                  # Edwards-Anderson Binder parameter calculations
├── entropies.jl              # Entanglement entropy computations
├── example_typed.jl          # Examples using the typed framework
└── README.md                  # This file
```
types.jl
**Purpose**: Define state type wrappers and polymorphic operations

**Key Types**:
- `PureStateMPS` - Pure quantum states |ψ⟩ (for quantum trajectories)
- `DiagonalStateMPS` - Diagonal density matrices (classical mixtures)
- `MixedStateMPS` - General mixed states (doubled MPS representation)

**Key Functions**:
- `norm(stateBasic weak quantum measurement protocols (for backward compatibility)

**Key Functions**:
- `create_weak_measurement_operators(sites, λₓ, λ_zz)` - Create Kraus operators
- `sample_and_apply(ψ, K0, K1, sites)` - Probabilistically apply measurements
- `compute_correlators(ψ, sites)` - Calculate 2-point and 4-point correlations

**Physical Model**:
Weak measurements with strength 0 ≤ λ ≤ 1:
- K₀ = (I + λM) / √(2(1 + λ²))
- K₁ = (I - λM) / √(2(1 + λ²))

### measurements_typed.jl
**Purpose**: Type-safe measurements with multiple dispatch

**Key Functions**:
- `expval(state::StateMPS, M, pos)` - Expectation value (dispatches on state type)
- `measure(state::StateMPS, M, λ, pos, outcome)` - Apply measurement with fixed outcome
- `measure_with_outcome(state::StateMPS, M, λ, pos)` - Sample measurement outcome

**Multiple Dispatch**:
Each function has specialized implementations for:
- `PureStateMPS`: ⟨ψ|M|ψ⟩ calculations
- `DiagonalStateMPS`: Tr(ρM) with diagonal ρ
- `MixedStateMPS`: Tr(ρM) with doubled MPS representation)` - General product state
- `zero_state(Type, L)` - Type-safe zero state constructor
- `ghz_state(Type, L)` - GHZ state (|00...0⟩ + |11...1⟩)/√2
- `neel_state(Type, L)` - Alternating up/down spins
- `create_up_state_mps(L)` - All spins up |↑↑...↑⟩
- `create_plus_state_mps(L)` - All spins in |+⟩ state
- `create_product_mps(L, state)` - General product state

### measurements.jl
**Purpose**: Implement weak quantum measurement protocols

**Key Functions**:
- `create_weak_measurement_operators(sites, λₓ, λ_zz)` - Create Kraus operators
- `sample_and_apply(ψ, K0, K1, sites)` - Probabilistically apply measurements
- `compute_correlators(ψ, sites)` - Calculate 2-point and 4-point correlations

**Physical Model**:
Weak measurements with strength 0 ≤ λ ≤ 1:
- K₀ = (I + λM) / √(2(1 + λ²))
- K₁ = (I - λM) / √(2(1 + λ²))

### dynamics.jl
**Purpose**: Time evolution under measurements and decoherence

**Key Functions**:
- `evolve_one_trial(L; λₓ, λ_zz)` - Standard measurement protocol
- `evolve_one_trial_dephasing(L; λₓ, λ_zz, Pₓ, P_zz)` - With dephasing channels
- `apply_single_site_gate(ψ, site, gate)` - Single-qubit operations
- `apply_two_site_gate(ψ, i, j, gate)` - Two-qubit operations
- `apply_x_dephasing(ψ, site, P)` - X-basis dephasing
- `apply_zz_dephasing(ψ, i, j, P)` - ZZ-basis dephasing

**Protocol**: T_max = 2L time steps (Backward Compatible)

```julia
using ITensors, ITensorMPS, ITensorCorrelators
include("src_new/binder.jl")

# Standard measurement protocol
result = ea_binder_mc(
    12;  # System size
    lambda_x = 0.5,
    lambda_zz = 0.5,
    ntrials = 1000,
    maxdim = 256,
    seed = 42
)

println("Binder parameter: ", result.B)
```

### Using the Typed Framework

```julia
include("src_new/types.jl")
include("src_new/state_constructors.jl")
include("src_new/measurements_typed.jl")

# Create a pure state
L = 10
state = zero_state(PureStateMPS, L)
println("Created: $state")

# Define measurement operator
X = [0 1; 1 0]  # Pauli X
λ = 0.5

# Compute expectation value
val = expval(state, X, 1)
println("⟨X₁⟩ = $val")

# Apply measurement with sampled outcome
state, outcome, val = measure_with_outcome(state, X, λ, 1)
println("Outcome: $outcome")

# The same interface works for all state types!
diag_state = zero_state(DiagonalStateMPS, L)
mixed_state = zero_state(MixedStateMPS, L)
```

### Multiple State Types

```julia
# Pure state trajectory
pure = zero_state(PureStateMPS, 8)
for i in 1:8
    pure, outcome, _ = measure_with_outcome(pure, X, 0.5, i)
end

# Diagonal density matrix (classical mixture)
diag = neel_state(DiagonalStateMPS, 8)
val = expval(diag, Z, 1)  # Tr(ρZ)

# Mixed density matrix (quantum mixture)
mixed = ghz_state(MixedStateMPS, 8)
val = expval(mixed, X, 1)  # Tr(ρX
M₂ = (1/L²) Σᵢⱼ |⟨ZᵢZⱼ⟩|²
M₄ = (1/L⁴) Σᵢⱼₖₗ |⟨ZᵢZⱼZₖZₗ⟩|²
```

### entropies.jl
**Purpose**: Compute entanglement measures

**Key Functions**:
- `calculate_bipartite_entropy(ψ, cut)` - von Neumann entropy across cut
- `weak_bipartite_entropy(L; λₓ, λ_zz, Pₓ, P_zz, ntrials)` - Average over trajectories

## Usage Examples

### Basic Binder Parameter Calculation

```julia
using ITensors, ITensorMPS, ITensorCorrelators
include("src_new/binder.jl")

# Standard measurement protocol
result = ea_binder_mc(
    12;  # System size
    lambda_x = 0.5,
    lambda_zz = 0.5,
    ntrials = 1000,
    maxdim = 256,
    seed = 42
)

println("Binder parameter: ", result.B)
```

### With Dephasing Channels

```julia
# Include dephasing noise
result = ea_binder_mc_dephasing(
    12;
    lambda_x = 0.5,
    lambda_zz = 0.5,
    P_x = 0.2,      # 20% X dephasing
    P_zz = 0.2,     # 20% ZZ dephasing
    ntrials = 1000
)
```

### Entanglement Entropy

```julia
include("src_new/entropies.jl")

entropy_result = weak_bipartite_entropy(
    12;
    lambda_x = 0.5,
    lambda_zz = 0.5,
    P_x = 0.0,
    P_zz = 0.0,
    ntrials = 100
)

println("Average entropy: ", entropy_result.entropy_mean)
```

## Design Principles

### Separation of Concerns
Each module handles a distinct aspect of the simulation:
- State preparation separate from dynamics
- Measurements separate from analysis
- Reusable components across different protocols

### Clear Dependencies
Module dependency flow:
```
state_constructors.jl
       ↓
measurements.jl
       ↓
dynamics.jl
       ↓
binder.jl / entropies.jl
```

### Consistency
- All functions have detailed docstrings
- Keyword arguments with sensible defaults
- Return values use named tuples for clarity
- Consistent parameter naming (lambda_x, lambda_zz, P_x, P_zz)

### Extensibility
Adding new protocols is straightforward:
1. Create new evolution function in `dynamics.jl`
2. Wrap it with analysis in `binder.jl` or similar
3. No need to modify existing modules

## Comparison with Original Code

### Before (src/)
- Monolithic modules: `BinderSim.jl` contained everything
- Code duplication between `BinderSim.jl` and `BinderSimDephasing.jl`
- Mixed concerns: state preparation, dynamics, and analysis together
- Difficult to test individual components

### After (src_new/)
- ✅ Modular design with single responsibilities
- ✅ Shared functions eliminate duplication
- ✅ Clear separation of concerns
- ✅ Easy to test and extend
- ✅ Better documentation
- ✅ Consistent naming conventions

## Migration Guide

To use the new structure in existing scripts:

**Old way:**
```julia
include("src/BinderSim.jl")
using .BinderSim
```

**New way:**
```julia
include("src_new/binder.jl")
# Functions are available directly
```

The new modules don't use module wrappers by default, making them simpler to include and use. Functions are exported and can be used directly after inclusion.

## Future Enhancements

Potential improvements:
- [ ] Add proper module wrapping if needed for package development
- [ ] Implement caching for measurement operators
- [ ] Add progress bars for long Monte Carlo runs
- [ ] Implement checkpointing for interrupted calculations
- [ ] Add more initial state constructors (Neel, GHZ, etc.)
- [ ] Support for periodic boundary conditions
- [ ] Generalized measurement protocols (non-uniform λ)

## Testing

Each module can be tested independently:

```julia
# Test state constructors
include("src_new/state_constructors.jl")
ψ, sites = create_up_state_mps(10)
@assert length(sites) == 10

# Test measurements
include("src_new/measurements.jl")
KX0, KX1, KZZ0, KZZ1 = create_weak_measurement_operators(sites, 0.5, 0.5)
@assert length(KX0) == 10

# Test dynamics (single trajectory)
include("src_new/dynamics.jl")
ψ_final, sites_final = evolve_one_trial(10; lambda_x=0.5, lambda_zz=0.5)
@assert length(ψ_final) == 10
```
