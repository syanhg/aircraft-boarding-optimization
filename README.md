# Optimizing Aircraft Boarding and Disembarkation

This repository contains the mathematical models and simulations for optimizing aircraft boarding and disembarkation processes through differential equations and numerical methods.

## Research Overview

The research applies fluid dynamics principles, differential equations, and numerical methods to model passenger movement during aircraft boarding and disembarkation. We focus on the Boeing 737-800 with a 3-3 seating configuration as our test case.

### Key Mathematical Approaches

1. **Fluid Dynamics Model**: We model passenger flow as a one-dimensional compressible fluid, using conservation laws to derive partial differential equations.

2. **Differential Equations**: The system is described by coupled PDEs representing passenger density and velocity evolution.

3. **Numerical Methods**:
   - Fourth-order Runge-Kutta method for high-accuracy solutions
   - Euler's method for validation and comparison
   - Lax-Friedrichs scheme for solving the advection-diffusion equation

4. **Bernoulli's Equation**: Applied to model flow acceleration in narrow sections of the aircraft aisle.

## Repository Contents

- **main.tex**: LaTeX source for the research paper
- **simulation.py**: Python implementation of the mathematical models
- Generated visualizations including:
  - Passenger density evolution plots
  - Strategy comparison charts
  - Numerical method comparisons

## Boarding Strategies Analyzed

1. **Random**: Passengers board in random order
2. **Back-to-Front**: Boarding from the back of the aircraft to the front
3. **Front-to-Back**: Boarding from the front of the aircraft to the back
4. **Window-Middle-Aisle (WMA)**: Passengers with window seats board first, followed by middle seats, then aisle seats
5. **Optimized**: Our proposed strategy combining WMA with back-to-front organization

## Running the Simulations

### Prerequisites

```
numpy
matplotlib
scipy
seaborn
```

### Execution

To run the simulations and generate all visualizations:

```
python simulation.py
```

This will:
- Compare all boarding strategies
- Generate density evolution plots for each strategy
- Produce numerical method comparison charts
- Analyze the Bernoulli effect in the aircraft aisle

## Key Findings

Our simulations demonstrate that:

1. The Window-Middle-Aisle strategy consistently outperforms traditional methods
2. Our optimized strategy provides a 22.5% improvement over random boarding
3. Front-to-back boarding is the least efficient method
4. The mathematical model accurately captures the queue formation and interference patterns observed in real boarding processes

## Paper Structure

The full mathematical development is presented in the LaTeX paper (main.tex), which includes:

1. Introduction and problem background
2. Mathematical formulation using fluid dynamics
3. Queueing theory integration
4. Differential equation system derivation
5. Numerical solution methods
6. Strategy evaluation and comparison
7. Results and analysis
8. Conclusions and practical implications

## Future Work

- Incorporating heterogeneous passenger behaviors
- Modeling multiple boarding doors
- Developing real-time adaptive strategies
- Extending the model to wide-body aircraft with multiple aisles