import numpy as np

# Parameters for the pde model
R = 1                                    # growth rate
q = 0.01                                    # deleterious mutation
alpha = 0.5                                 # parameter alpha
D = np.sqrt(alpha*(q**3)/R)                # beneficial/deleterious mutations

Delta_x = 0.1                                # space discretisation step
Delta_t = 0.01                               # time discretisation step
Ld = -30                                     # left boundary
Lb = 30                                      # right boundary
K = 100                                      # Carrying capacity   
Tfinal = 10000                                  # Final time
N0 = 50                                      # Initial population size

#growth_indicator = 'exp'                    # No Carrying Capacity is included
growth_indicator = 'CC'                     # Model with a carrying capacity

# Parameters for the Gillespie
Deltat = 1.0/(10*R*180)                      # fixed timestep for Gillespie
nsimulations = 20#1000                           # Number of independent simulations   
h_compartment = 0.8                          # Compartment size for visualising the results   
threshold = 5                               # Threshold for death or not


filename = 'Results/' + 'alpha=' + str(alpha) + 'R=' + str(round(R,3)) + 'q=' + str(round(q,3)) + 'N0=' + str(N0) + 'K=' + str(K) 
