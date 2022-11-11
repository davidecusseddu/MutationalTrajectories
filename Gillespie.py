import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from GillespieFunctions import *


begin_index = 1
X0 = [i * 0.0 for i in list(range(1, N0+1))]  # List with all positions. Each component represents the position of a given individual

PDEmodel_tot_pop = np.loadtxt(filename +'-integraln.txt')
time_PDEmodel = np.loadtxt(filename + '-tmesh.txt')
Tfinal = len(PDEmodel_tot_pop)*Delta_t

Nt = int(Tfinal/Deltat)

mean_gaussian_mutation_distribution = 0.0

X = [i * 0.0 for i in list(range(1, N0+1))]
total_population = len(X) 
sqrt_2D_Deltat = np.sqrt(2*D*Deltat)
q_Deltat = q*Deltat
final_populations = []

saving_step = 5
population_history = []
times = Deltat*np.arange(0,Nt,saving_step)
totalpopulationsize = np.zeros([len(times)+1, nsimulations])
extinction = np.zeros(nsimulations, dtype = int)
extinction_time = Deltat*(Nt+1)*np.ones(nsimulations, dtype = int)
population_extension = np.zeros(nsimulations, dtype = int)
saving_indices = np.zeros(nsimulations, dtype=int)

for simulation in range(nsimulations):
    print('simulation ' + format(simulation, '04d') + ' out of ', nsimulations,  end ='\r')

    X = [i * 0.0 for i in list(range(1, N0+1))]
    total_population = len(X)
    totalpopulationsize[0, simulation] = total_population

    timestep = 0
    saving_index = 0
            
    extension_flag = 0
    
    Nt_simulation = Nt
    
    while timestep < Nt_simulation and total_population > 0:

        total_population_variation = 0  
        offspringX = []
        deadX = []
        
        ## CARRYING CAPACITY
        if growth_indicator == 'CC':
            crowding_percentage = total_population/K
            logistic = 1 - crowding_percentage
        ## PURELY EXPONENTIAL
        else:
            logistic = 1
            crowding_percentage = 1
        
        for individual in range(total_population-1,-1,-1):

            xi = np.random.normal(loc=mean_gaussian_mutation_distribution, scale=1.0, size=None)
            new_genotype_location = X[individual] + sqrt_2D_Deltat*xi - q_Deltat
            X[individual] = new_genotype_location

            r1 = np.random.uniform(0,1)

            if new_genotype_location>0: # growth
                event_rate = R*new_genotype_location*logistic*Deltat
            else: # death 
                event_rate = R*new_genotype_location*crowding_percentage*Deltat

            if r1 < abs(event_rate):
                if event_rate < 0:
                    deadX.append(individual)
                    total_population_variation -= 1
                else: 
                    offspring_location = new_genotype_location
                    offspringX.append(offspring_location)
                    total_population_variation += 1

        total_population = total_population + total_population_variation

        for dead_individual in deadX:
            del X[dead_individual]
            
        X = X + offspringX

        if timestep%saving_step == 0 and extension_flag == 0:
            saving_index += 1
            population_history.append(X)
            totalpopulationsize[saving_index, simulation] = total_population
        
        
        if len(X) == 0:
            extinction[simulation] = 1
            extinction_time[simulation] = timestep*Deltat
        
        timestep += 1

        if total_population < 3 and timestep == Nt_simulation-1:     # checks that if total population is less than 3 individuals, it either recovers or die out
            Nt_simulation +=1
            extension_flag = 1
            population_extension[simulation] += 1 
    
    
    final_populations.append(X)
    saving_indices[simulation] = saving_index


threshold = 1
mean_tot_pop = mean_gillespie_surviving_replicates(totalpopulationsize, times, nsimulations, threshold)

time_PDEmodel = time_PDEmodel[:len(PDEmodel_tot_pop)]

plt.figure('totalpopulations', (13,8))
for simulation in range(nsimulations):
    plt.plot(times,totalpopulationsize[:-1, simulation], alpha=0.2 )
plt.plot(times,mean_tot_pop, LineWidth= 4, label = 'mean with all survivors', color = 'red')
plt.plot(times,np.mean(totalpopulationsize[:-1, :],1), LineWidth= 2, label = 'mean between all replicates', color = 'blue' )
plt.plot(time_PDEmodel,PDEmodel_tot_pop,'r', LineWidth= 2,  label = 'pde model', color = 'green' )
plt.legend()
#plt.yscale('log')
#plt.ylim([0,150])
#plt.xlim([0,130])
plt.xlabel('time t',fontsize = 16)
plt.ylabel('total population size',fontsize = 16)
#plt.title( str(np.sum(extinction)) + ' extinctions out of ' + str(nsimulations) + ' simulations in time t = ' + str(round(Nt*Delta_t,2)) ,fontsize = 16)
plt.title( ' extinctions: ' + str(np.sum(extinction)) + '/' + str(nsimulations) + '. R = ' + str(round(R,4)) + ', q = ' + str(q) + ', D = ' + str(round(D,4))+ ', alpha = ' + str(alpha) ,fontsize = 16)
plt.savefig(filename + "-totalpopulations.png")

plt.figure('totalpopulations-zoom', (13,8))
for simulation in range(nsimulations):
    plt.plot(times,totalpopulationsize[:-1, simulation], alpha=0.2 )
plt.plot(times,mean_tot_pop, LineWidth= 4, label = 'mean with all survivors', color = 'red')
plt.plot(times,np.mean(totalpopulationsize[:-1, :],1), LineWidth= 2, label = 'mean between all replicates', color = 'blue' )
plt.plot(time_PDEmodel,PDEmodel_tot_pop,'r', LineWidth= 2,  label = 'pde model', color = 'green' )
plt.legend()
#plt.yscale('log')
#plt.ylim([0,150])
plt.xlim([0,130])
plt.xlabel('time t',fontsize = 16)
plt.ylabel('total population size',fontsize = 16)
#plt.title( str(np.sum(extinction)) + ' extinctions out of ' + str(nsimulations) + ' simulations in time t = ' + str(round(Nt*Delta_t,2)) ,fontsize = 16)
plt.title( ' extinctions: ' + str(np.sum(extinction)) + '/' + str(nsimulations) + '. R = ' + str(round(R,4)) + ', q = ' + str(q) + ', D = ' + str(round(D,4))+ ', alpha = ' + str(alpha) ,fontsize = 16)
plt.savefig(filename + "-totalpopulations-zoomed.png")


plt.figure('population distribution', (13,8))
PDEmodel_last_pop = np.loadtxt(filename + '-lastn.txt')
PDEmodel_xmesh = np.loadtxt(filename + '-xmesh.txt')
number_of_survivals = 0
h_compartment = 1
x_compartments,y_mean=average_distribution(final_populations, h_compartment)

survivors = [x for x in final_populations if len(x) > 0]
if len(survivors)>0:
    MINX = min(min(survivors, key = min))
    MAXX = max(max(survivors, key = max))
else:
    MINX = Ld
    MAXX = Lb
    

N_compartments = int((MAXX - MINX)/h_compartment)
summa = np.zeros(N_compartments)

for simulation in range(nsimulations):
    if final_populations[simulation]:
        xlocation, bars = plot_distribution(final_populations[simulation], h_compartment, MINX, MAXX)
        plt.plot(xlocation, bars, linewidth = 0.5, alpha = 0.4)
        summa = summa + bars
        number_of_survivals += 1
        
if number_of_survivals == 0:
    xlocation = Ld + (0.5 + np.arange(N_compartments))*h_compartment
Ih_y = integral_representation(PDEmodel_xmesh,PDEmodel_last_pop,xlocation,h_compartment)

plt.plot(xlocation,Ih_y, LineWidth= 3, color='green', label = 'pde model')
if number_of_survivals > 0:
    plt.plot(xlocation,summa/number_of_survivals, 'red', LineWidth= 2,label = 'mean simulations' )

plt.legend()
#plt.xlim([1.3*MINX, 1.3*MAXX])
plt.ylim([0, 100])
plt.title('time ' + str(round(Nt*Deltat,2)))
plt.title( ' extinctions: ' + str(np.sum(extinction)) + '/' + str(nsimulations) + '. R = ' + str(R) + ', q = ' + str(q) + ', D = ' + str(round(D,4))+ ', alpha = ' + str(alpha) ,fontsize = 16)
plt.xlabel('genotype space', fontsize = 16)
plt.ylabel('population distribution', fontsize = 16)
plt.savefig(filename + "-populationdistribution.png")

print(number_of_survivals, ' populations survived on ', nsimulations, ' independent simulated populations' )

begin_index += 1


import os
os.system("cp ./parameters.py ./"+filename+"-parameters.py")

