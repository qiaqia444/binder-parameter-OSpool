#!/usr/bin/env julia
"""
Example demonstrating the type-safe measurement framework.

This script shows how to use the typed state representations and
measurements with multiple dispatch.
"""

using ITensors, ITensorMPS
using Random, Statistics

# Include the typed modules
include("types.jl")
include("state_constructors.jl")
include("measurements_typed.jl")

using .Main: PureStateMPS, DiagonalStateMPS, MixedStateMPS
using .Main: zero_state, ghz_state, neel_state
using .Main: expval, measure, measure_with_outcome

function demo_pure_state()
    println("\n" * "="^70)
    println("DEMO: Pure State Measurements")
    println("="^70)
    
    L = 8
    state = zero_state(PureStateMPS, L)
    println("Created pure state: $state")
    
    # Define Pauli X operator
    X = [0 1; 1 0]
    
    # Measure at position 1 with λ=0.5
    λ = 0.5
    val = expval(state, X, 1)
    println("⟨X₁⟩ = $val")
    
    # Apply measurement with outcome 0
    state = measure(state, X, λ, 1, false)
    println("After measurement: $state")
    
    # Sample measurement outcome
    state, outcome, val = measure_with_outcome(state, X, λ, 2)
    println("Measured site 2: outcome=$outcome, ⟨X₂⟩=$val")
    
    return state
end

function demo_diagonal_state()
    println("\n" * "="^70)
    println("DEMO: Diagonal State Measurements")
    println("="^70)
    
    L = 6
    state = neel_state(DiagonalStateMPS, L)
    println("Created diagonal state: $state")
    
    Z = [1 0; 0 -1]
    
    # Expectation value on diagonal state
    val = expval(state, Z, 1)
    println("⟨Z₁⟩ = $val")
    
    # Apply measurement
    λ = 0.3
    state = measure(state, Z, λ, 1, true)
    println("After measurement: $state")
    
    return state
end

function demo_mixed_state()
    println("\n" * "="^70)
    println("DEMO: Mixed State Measurements")
    println("="^70)
    
    L = 4
    state = ghz_state(MixedStateMPS, L)
    println("Created mixed state: $state")
    
    X = [0 1; 1 0]
    
    # Expectation value on mixed state
    val = expval(state, X, 1)
    println("Tr(ρX₁) = $val")
    
    # Apply measurement
    λ = 0.5
    state = measure(state, X, λ, 1, false)
    println("After measurement: $state")
    
    return state
end

function demo_batch_measurements()
    println("\n" * "="^70)
    println("DEMO: Batch Measurements")
    println("="^70)
    
    L = 10
    state = zero_state(PureStateMPS, L)
    
    X = [0 1; 1 0]
    λ = 0.5
    
    # Measure all sites
    positions = 1:L
    state, outcomes, expvals = measure_with_outcome(state, X, λ, positions)
    
    println("Measured $L sites")
    println("Outcomes: $outcomes")
    println("Mean expectation value: $(mean(real.(expvals)))")
    
    return state
end

function compare_state_types()
    println("\n" * "="^70)
    println("COMPARISON: Different State Types")
    println("="^70)
    
    L = 6
    λ = 0.5
    X = [0 1; 1 0]
    
    # Pure state
    pure = zero_state(PureStateMPS, L)
    pure = measure(pure, X, λ, 1, false)
    println("Pure state after measurement: norm = $(norm(pure))")
    
    # Diagonal state
    diag = zero_state(DiagonalStateMPS, L)
    diag = measure(diag, X, λ, 1, false)
    println("Diagonal state after measurement: trace = $(norm(diag))")
    
    # Mixed state
    mixed = zero_state(MixedStateMPS, L)
    mixed = measure(mixed, X, λ, 1, false)
    println("Mixed state after measurement: trace = $(norm(mixed))")
    
    println("\nAll state types support the same measurement interface!")
end

function main()
    println("\n" * "="^70)
    println("TYPE-SAFE MEASUREMENT FRAMEWORK EXAMPLES")
    println("="^70)
    
    Random.seed!(42)
    
    demo_pure_state()
    demo_diagonal_state()
    demo_mixed_state()
    demo_batch_measurements()
    compare_state_types()
    
    println("\n" * "="^70)
    println("DEMOS COMPLETE")
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
