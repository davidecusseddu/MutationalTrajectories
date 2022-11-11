import numpy as np

def plot_distribution(X, h_compartment, min_X, max_X):
    N_compartments = int((max_X - min_X)/h_compartment)
    yvalue = np.zeros(N_compartments)
    xlocation = np.zeros(N_compartments)
    
    for i in range(N_compartments):
        xlocation[i] = min_X +  (i+0.5)*h_compartment
        for j in range(len(X)):
            if X[j] >= min_X + i*h_compartment and X[j] < min_X + (i+1)*h_compartment:
                yvalue[i] += 1
    return xlocation, yvalue

    
def average_distribution(final_populations, h_compartment):
    N_simulations = len(final_populations)
    min_final_pop = 1e10
    MAX_final_pop = -1e10
    
    for replicate in range(N_simulations):
        if len(final_populations[replicate])>0:
            if min(final_populations[replicate])< min_final_pop:
                min_final_pop = min(final_populations[replicate])
            if max(final_populations[replicate])> MAX_final_pop:
                MAX_final_pop = max(final_populations[replicate])
    print(min_final_pop)
    print(MAX_final_pop)
    N_compartments = int((MAX_final_pop - min_final_pop)/h_compartment)
    xlocation = np.arange(min_final_pop,MAX_final_pop,h_compartment)
    yvalues = np.zeros(len(xlocation))
    
    for replicate in range(N_simulations):

        X = final_populations[replicate]
        
        if len(X)>0:
            for j in range(len(X)):
                for i in range(N_compartments-1):
                    if X[j] >= xlocation[i] and X[j] < xlocation[i+1]:
                        yvalues[i] += 1
                        break
            
                    
    return xlocation, yvalues/N_simulations


def integral_representation(xmesh,y,xlocation,hcompartment):
    integraly_next_compartment = 0
    Ih_y = np.zeros(len(xlocation))

    for i in range(len(xlocation)):
        integraly_compartment = 0 + integraly_next_compartment
        integraly_next_compartment = 0
        
        A = xlocation[i]-0.5*hcompartment
        B = xlocation[i]+0.5*hcompartment
        
        for j in range(len(xmesh)-1):

            if xmesh[j]>=A and xmesh[j+1]<B:
                integraly_compartment += 0.5*(y[j]+y[j+1])*(xmesh[j+1]-xmesh[j])
            elif xmesh[j]>=A and xmesh[j]<B and xmesh[j+1]>B:
                yB = (y[j+1]-y[j])*(B-xmesh[j])/(xmesh[j+1]-xmesh[j]) + y[j]
                integraly_compartment += 0.5*(y[j]+yB)*(B-xmesh[j])
                integraly_next_compartment = 0.5*(y[j+1]+yB)*(xmesh[j+1]-B)
        
        Ih_y[i] = integraly_compartment
        
    return Ih_y


def mean_gillespie_surviving_replicates(totalpopulationsize, times, nsimulations,threshold):
    mean_tot_pop = np.zeros(len(times)) 
    n_survivors = 0 
    for simulation in range(nsimulations):
        if totalpopulationsize[-1, simulation] > threshold:   
            for t in range(len(times)):
                mean_tot_pop[t] = mean_tot_pop[t] + totalpopulationsize[t, simulation]
            n_survivors += 1
    if n_survivors > 0:        
        mean_tot_pop = mean_tot_pop/n_survivors
    return mean_tot_pop

