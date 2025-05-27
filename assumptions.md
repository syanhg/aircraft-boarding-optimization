# Extended Model Assumptions

This document outlines the complete set of assumptions used in our mathematical model of aircraft boarding and disembarkation processes.

## Physical Constraints

1. **Aircraft Configuration**
   - Boeing 737-800 with 3-3 seating configuration (single aisle)
   - Total of 33 rows (126 seats)
   - Prestige class in rows 7-9 (12 seats)
   - Economy class in remaining rows (114 seats)
   - Aisle width sufficient for single-file passenger movement only

2. **Spatial Parameters**
   - Aircraft cabin length: 30 meters
   - Row spacing: approximately 0.76 meters (30 inches)
   - Aisle width: 0.48 meters (19 inches)

## Passenger Movement Dynamics

1. **Unidirectional Flow**
   - All passenger movement occurs in a single direction
   - No overtaking or position swapping during movement
   - Movement toward the back during boarding, toward the exits during disembarkation
   - No path reversal or deviations to incorrect seats allowed

2. **Movement Characteristics**
   - All passengers move at a uniform base speed (1.2 m/s under free-flow conditions)
   - Speed reduces as a function of passenger density (Greenshields model)
   - No unnecessary stopping except for essential actions (stowing luggage, seating, etc.)

3. **Queue Formation**
   - Passengers form queues based on a first-in-first-out (FIFO) principle
   - Queue density at entrance follows a Poisson distribution
   - Queue formation within the aircraft follows conservation laws

## Mathematical Simplifications

1. **Continuous Flow Approximation**
   - The discrete process of individual passengers boarding is approximated as a continuous fluid flow
   - This allows application of fluid dynamics principles and PDEs
   - Passenger density (ρ) treated as a continuous variable over space

2. **Deterministic Base Parameters with Stochastic Variations**
   - Base seating time: mean μₛ = 15 seconds, variance σₛ² = 25 seconds²
   - Passenger attributes follow normal distribution N(μₐ, σₐ²)
   - Luggage stowing time proportional to number of carry-on items

3. **Conservation Principles**
   - Total passenger count remains constant throughout the process
   - Conservation equation: ∂ρ/∂t + ∂(ρv)/∂x = S(x,t)
   - Source/sink term S(x,t) represents passengers entering/leaving the aisle

## Boarding Process Specifics

1. **Boarding Strategies**
   - Random: passengers board in random order
   - Back-to-Front: boarding groups organized by row sections, from back to front
   - Front-to-Back: boarding groups organized by row sections, from front to back
   - Window-Middle-Aisle: passengers with window seats board first, followed by middle, then aisle
   - Optimized: combination of WMA with back-to-front organization

2. **Boarding Rate Parameters**
   - Passenger arrival rate at aircraft entrance: 0.5 passengers/second
   - Maximum processing rate at boarding gate: 0.6 passengers/second
   - Group size for zone-based boarding: 5 rows per group

3. **Interference Factors**
   - Aisle interference: passengers in aisle block others from reaching their seats
   - Seat interference: passengers must get up to let others reach window/middle seats
   - Luggage stowing interference: passengers stowing luggage temporarily block the aisle

## Disembarkation Process Specifics

1. **Initial Conditions**
   - All passengers initially seated
   - Exit rate limited by door capacity
   - Row-by-row emptying pattern

2. **Behavioral Parameters**
   - Standing and luggage retrieval time: mean 10 seconds, variance 16 seconds²
   - Exit door processing capacity: 0.6 passengers/second

## Simplifications and Limitations

1. **Excluded Factors**
   - No flight attendant movement or interference
   - No consideration of special boarding for passengers needing assistance
   - No bathroom usage during boarding/disembarkation
   - No psychological factors (stress, urgency perception)

2. **Steady-State Assumptions**
   - Constant environmental conditions
   - No external disruptions to the process
   - Consistent crew performance

These assumptions allow us to create a tractable mathematical model while maintaining sufficient realism to derive meaningful insights for practical aircraft boarding optimization.