import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import poisson
import seaborn as sns

class AircraftBoardingSimulation:
    """
    A comprehensive simulation of aircraft boarding processes using fluid dynamics and 
    differential equations.
    """
    
    def __init__(self, 
                 aircraft_length=30,       # Length of aircraft in meters
                 num_rows=33,              # Number of rows
                 seats_per_row=6,          # Seats per row (3-3 configuration)
                 free_flow_velocity=1.2,   # Maximum walking speed in m/s
                 jam_density=3.5,          # Maximum passenger density in pax/m
                 mean_seating_time=15,     # Average time to stow luggage and sit (seconds)
                 seating_time_variance=25, # Variance in seating time (seconds²)
                 dx=0.5,                   # Spatial discretization step (meters)
                 dt=0.1):                  # Time discretization step (seconds)
        
        # Initialize parameters
        self.L = aircraft_length
        self.num_rows = num_rows
        self.seats_per_row = seats_per_row
        self.v_free = free_flow_velocity
        self.rho_jam = jam_density
        self.mu_s = mean_seating_time
        self.sigma_s_squared = seating_time_variance
        self.dx = dx
        self.dt = dt
        
        # Setup spatial grid
        self.x = np.arange(0, self.L + self.dx, self.dx)
        self.nx = len(self.x)
        
        # Initialize density and velocity arrays
        self.rho = np.zeros(self.nx)
        self.v = np.zeros(self.nx)
        
        # Row positions (evenly spaced)
        self.row_positions = np.linspace(5, self.L - 5, self.num_rows)
        
        # Passenger assignments (seat number for each passenger)
        self.seat_assignments = np.zeros(self.num_rows * self.seats_per_row, dtype=int)
        
        # Time counter
        self.t = 0
        
        # Results storage
        self.time_history = []
        self.density_history = []
    
    def greenshields_velocity(self, rho):
        """
        Calculate velocity based on Greenshields model.
        
        Parameters:
        -----------
        rho : array_like
            Passenger density at each spatial point
            
        Returns:
        --------
        array_like
            Velocity at each spatial point
        """
        return self.v_free * (1 - np.minimum(rho / self.rho_jam, 1))
    
    def assign_boarding_strategy(self, strategy='random'):
        """
        Assign boarding order based on the selected strategy.
        
        Parameters:
        -----------
        strategy : str
            Boarding strategy: 'random', 'back_to_front', 'front_to_back', 'wma'
        """
        n_passengers = self.num_rows * self.seats_per_row
        
        # Create passenger indices
        passengers = np.arange(n_passengers)
        
        # Create row and seat indices for each passenger
        rows = np.repeat(np.arange(self.num_rows), self.seats_per_row)
        seats = np.tile(np.arange(self.seats_per_row), self.num_rows)
        
        # Map seat indices to seat types (0=window, 1=middle, 2=aisle)
        seat_types = np.zeros_like(seats)
        seat_types[seats == 0] = 0  # window (left)
        seat_types[seats == 5] = 0  # window (right)
        seat_types[seats == 1] = 1  # middle (left)
        seat_types[seats == 4] = 1  # middle (right)
        seat_types[seats == 2] = 2  # aisle (left)
        seat_types[seats == 3] = 2  # aisle (right)
        
        if strategy == 'random':
            # Random boarding
            np.random.shuffle(passengers)
            
        elif strategy == 'back_to_front':
            # Back to front (by groups of rows)
            group_size = 5  # Number of rows per group
            num_groups = int(np.ceil(self.num_rows / group_size))
            
            # Initialize ordering
            ordering = np.zeros_like(passengers)
            
            # Assign priorities by row groups (back to front)
            for g in range(num_groups):
                group_rows = np.arange(
                    max(0, self.num_rows - (g+1)*group_size), 
                    self.num_rows - g*group_size
                )
                
                # Find passengers in these rows
                group_mask = np.isin(rows, group_rows)
                ordering[group_mask] = g
            
            # Random ordering within each group
            for g in range(num_groups):
                group_indices = np.where(ordering == g)[0]
                np.random.shuffle(group_indices)
                passengers[group_indices] = passengers[group_indices]
                
        elif strategy == 'front_to_back':
            # Front to back (by groups of rows)
            group_size = 5  # Number of rows per group
            num_groups = int(np.ceil(self.num_rows / group_size))
            
            # Initialize ordering
            ordering = np.zeros_like(passengers)
            
            # Assign priorities by row groups (front to back)
            for g in range(num_groups):
                group_rows = np.arange(
                    g*group_size, 
                    min(self.num_rows, (g+1)*group_size)
                )
                
                # Find passengers in these rows
                group_mask = np.isin(rows, group_rows)
                ordering[group_mask] = g
            
            # Random ordering within each group
            for g in range(num_groups):
                group_indices = np.where(ordering == g)[0]
                np.random.shuffle(group_indices)
                passengers[group_indices] = passengers[group_indices]
                
        elif strategy == 'wma':
            # Window-Middle-Aisle
            
            # First sort by seat type
            sort_indices = np.lexsort((np.random.random(len(passengers)), seat_types))
            
            # Apply the sorting
            passengers = passengers[sort_indices]
        
        elif strategy == 'optimized':
            # Our proposed optimized strategy
            # Combination of back-to-front with window-middle-aisle priority
            
            # First sort by row (back to front)
            row_priority = self.num_rows - rows
            
            # Then by seat type
            sort_indices = np.lexsort((np.random.random(len(passengers)), seat_types, row_priority))
            
            # Apply the sorting
            passengers = passengers[sort_indices]
        
        else:
            raise ValueError(f"Unknown boarding strategy: {strategy}")
        
        # Assign the boarding order
        self.boarding_order = passengers
        
    def initialize_simulation(self, arrival_rate=0.5):
        """
        Initialize the simulation with passengers at the entrance.
        
        Parameters:
        -----------
        arrival_rate : float
            Rate at which passengers arrive at the aircraft entrance (passengers/second)
        """
        # Reset simulation variables
        self.rho = np.zeros(self.nx)
        self.v = np.zeros(self.nx)
        self.t = 0
        self.time_history = []
        self.density_history = []
        
        # Set initial condition (passengers at entrance)
        self.rho[0] = arrival_rate / self.v_free
        
        # Initialize velocity using Greenshields model
        self.v = self.greenshields_velocity(self.rho)
        
        # Store initial state
        self.time_history.append(self.t)
        self.density_history.append(self.rho.copy())
    
    def step(self):
        """
        Advance the simulation by one time step using Lax-Friedrichs method.
        """
        # Calculate velocity based on current density
        self.v = self.greenshields_velocity(self.rho)
        
        # Calculate flux
        flux = self.rho * self.v
        
        # Initialize new density array
        rho_new = np.zeros_like(self.rho)
        
        # Apply Lax-Friedrichs scheme for interior points
        for i in range(1, self.nx-1):
            # Advection term
            rho_new[i] = 0.5 * (self.rho[i+1] + self.rho[i-1]) - \
                         0.5 * self.dt/self.dx * (flux[i+1] - flux[i-1])
            
            # Source/sink term (passengers taking their seats)
            # Find closest row
            closest_row_idx = np.argmin(np.abs(self.x[i] - self.row_positions))
            closest_row_dist = np.min(np.abs(self.x[i] - self.row_positions))
            
            # If close enough to a row and there are passengers to be seated
            if closest_row_dist < self.dx and self.rho[i] > 0:
                # Probability of taking a seat in this row
                seat_prob = self.dt / self.mu_s
                
                # Apply seating (subtract from density)
                seating_rate = min(seat_prob * self.rho[i], self.rho[i])
                rho_new[i] -= seating_rate
        
        # Boundary conditions
        # Entrance (Dirichlet condition)
        rho_new[0] = self.rho[0]
        
        # Exit (Neumann condition)
        rho_new[-1] = rho_new[-2]
        
        # Update density
        self.rho = rho_new
        
        # Update time
        self.t += self.dt
        
        # Store state
        self.time_history.append(self.t)
        self.density_history.append(self.rho.copy())
    
    def run_simulation(self, strategy='random', max_time=2000):
        """
        Run the full simulation with the specified strategy.
        
        Parameters:
        -----------
        strategy : str
            Boarding strategy: 'random', 'back_to_front', 'front_to_back', 'wma', 'optimized'
        max_time : float
            Maximum simulation time (seconds)
            
        Returns:
        --------
        float
            Total boarding time
        """
        # Assign boarding strategy
        self.assign_boarding_strategy(strategy)
        
        # Initialize simulation
        self.initialize_simulation()
        
        # Run simulation until all passengers are seated or max time is reached
        while self.t < max_time and np.sum(self.rho) > 1e-6:
            self.step()
        
        # Return total boarding time
        return self.t
    
    def compare_strategies(self, n_simulations=10):
        """
        Compare different boarding strategies.
        
        Parameters:
        -----------
        n_simulations : int
            Number of simulations to run for each strategy
            
        Returns:
        --------
        dict
            Dictionary with statistics for each strategy
        """
        strategies = ['random', 'back_to_front', 'front_to_back', 'wma', 'optimized']
        results = {strategy: [] for strategy in strategies}
        
        for strategy in strategies:
            for _ in range(n_simulations):
                boarding_time = self.run_simulation(strategy)
                results[strategy].append(boarding_time)
            
            print(f"{strategy.capitalize()} strategy: "
                  f"Mean = {np.mean(results[strategy]):.2f}s, "
                  f"SD = {np.std(results[strategy]):.2f}s")
        
        return results
    
    def plot_density_evolution(self, strategy='random'):
        """
        Plot the evolution of passenger density over time for a given strategy.
        """
        # Run a simulation
        self.run_simulation(strategy)
        
        # Convert history to array
        density_history = np.array(self.density_history)
        time_history = np.array(self.time_history)
        
        # Create a time-space plot
        plt.figure(figsize=(10, 6))
        
        # Downsample for plotting
        stride = max(1, len(time_history) // 100)
        
        X, T = np.meshgrid(self.x, time_history[::stride])
        
        plt.pcolormesh(X, T, density_history[::stride], 
                       cmap='viridis', shading='auto')
        plt.colorbar(label='Passenger Density (pax/m)')
        plt.xlabel('Position in Aircraft (m)')
        plt.ylabel('Time (s)')
        plt.title(f'Passenger Density Evolution - {strategy.capitalize()} Strategy')
        
        # Mark row positions
        for row_pos in self.row_positions:
            plt.axvline(x=row_pos, color='r', linestyle='--', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'density_evolution_{strategy}.png', dpi=300)
        plt.close()
    
    def plot_strategy_comparison(self, results):
        """
        Plot comparison of boarding strategies.
        
        Parameters:
        -----------
        results : dict
            Dictionary with simulation results for each strategy
        """
        plt.figure(figsize=(10, 6))
        
        # Create violin plots
        sns.violinplot(data=[results[strategy] for strategy in results.keys()], 
                        inner='box')
        plt.xticks(range(len(results)), 
                  [s.replace('_', '-').capitalize() for s in results.keys()])
        plt.ylabel('Boarding Time (s)')
        plt.title('Comparison of Boarding Strategies')
        
        # Add mean values as text
        for i, strategy in enumerate(results.keys()):
            mean_val = np.mean(results[strategy])
            plt.text(i, mean_val + 50, f'{mean_val:.0f}s', 
                     ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('strategy_comparison.png', dpi=300)
        plt.close()
    
    def generate_runge_kutta_comparison(self):
        """
        Demonstrate the difference between Euler's method and Runge-Kutta for solving 
        a simplified version of our differential equations.
        """
        # Define a simplified ODE system for boarding
        def ode_system(t, y):
            rho, v = y
            # Simplified model (ignoring spatial dependence)
            drho_dt = -0.1 * rho  # Passengers exit the system at a rate proportional to density
            dv_dt = 0.05 * (self.v_free - v)  # Velocity approaches free-flow velocity
            return [drho_dt, dv_dt]
        
        # Initial conditions
        rho0 = 2.0  # Initial passenger density
        v0 = 0.5    # Initial velocity
        
        # Time points
        t_span = (0, 100)
        t_eval = np.linspace(0, 100, 1000)
        
        # Solve using RK45
        sol_rk45 = solve_ivp(ode_system, t_span, [rho0, v0], 
                              method='RK45', t_eval=t_eval)
        
        # Solve using Euler's method
        dt_euler = 1.0
        t_euler = np.arange(0, 100, dt_euler)
        sol_euler = np.zeros((2, len(t_euler)))
        sol_euler[0, 0] = rho0
        sol_euler[1, 0] = v0
        
        for i in range(1, len(t_euler)):
            derivatives = ode_system(t_euler[i-1], sol_euler[:, i-1])
            sol_euler[:, i] = sol_euler[:, i-1] + dt_euler * np.array(derivatives)
        
        # Plot comparison
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.plot(sol_rk45.t, sol_rk45.y[0], 'b-', label='RK45')
        plt.plot(t_euler, sol_euler[0], 'r--', label='Euler')
        plt.xlabel('Time (s)')
        plt.ylabel('Passenger Density')
        plt.legend()
        plt.title('Comparison of Numerical Methods')
        
        plt.subplot(1, 2, 2)
        plt.plot(sol_rk45.t, sol_rk45.y[1], 'b-', label='RK45')
        plt.plot(t_euler, sol_euler[1], 'r--', label='Euler')
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (m/s)')
        plt.legend()
        
        plt.tight_layout()
        plt.savefig('numerical_methods_comparison.png', dpi=300)
        plt.close()
    
    def analyze_bernoulli_effect(self):
        """
        Analyze the effect of applying Bernoulli's equation to model flow acceleration
        in narrow sections of the aisle.
        """
        # Define aisle widths (normalized)
        x = np.linspace(0, self.L, 100)
        width = np.ones_like(x)
        
        # Create a narrowing in the middle
        narrow_section = (x > 12) & (x < 18)
        width[narrow_section] = 0.8
        
        # Calculate velocity using Bernoulli's principle
        # Assuming constant flow rate, v ∝ 1/A
        v_normalized = 1 / width
        
        # Plot
        plt.figure(figsize=(10, 6))
        
        plt.subplot(2, 1, 1)
        plt.plot(x, width)
        plt.ylabel('Relative Aisle Width')
        plt.title('Bernoulli Effect in Aircraft Aisle')
        plt.grid(True)
        
        plt.subplot(2, 1, 2)
        plt.plot(x, v_normalized)
        plt.xlabel('Position in Aircraft (m)')
        plt.ylabel('Relative Velocity')
        plt.grid(True)
        
        plt.tight_layout()
        plt.savefig('bernoulli_effect.png', dpi=300)
        plt.close()

if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Create simulation
    sim = AircraftBoardingSimulation()
    
    # Compare different strategies
    results = sim.compare_strategies(n_simulations=20)
    
    # Plot strategy comparison
    sim.plot_strategy_comparison(results)
    
    # Plot density evolution for each strategy
    for strategy in ['random', 'back_to_front', 'front_to_back', 'wma', 'optimized']:
        sim.plot_density_evolution(strategy)
    
    # Generate numerical method comparison
    sim.generate_runge_kutta_comparison()
    
    # Analyze Bernoulli effect
    sim.analyze_bernoulli_effect()
    
    print("Simulation completed. Results and plots have been saved.")