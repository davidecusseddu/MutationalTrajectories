import matplotlib
matplotlib.use('Agg')

# Code for PDE model
# below the code for thomas algorithm and  crout algorithm for solving the tridiagonal system and calculation of the integrals
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
now = datetime.now()

from parameters import *
from PDEmodelFunctions import *




# filename = 'Results/' + input('Give me the file name = ')

x_mesh_minus = np.arange(Ld,0.0,Delta_x)
x_mesh_plus = np.arange(0.0,Lb,Delta_x)
if x_mesh_minus[-1] == 0.0:
    x_mesh = np.concatenate((x_mesh_minus[0:-1],x_mesh_plus))
else:
    x_mesh = np.concatenate((x_mesh_minus,x_mesh_plus))

i_0 = np.where(x_mesh == 0)    

t_mesh = np.arange(0,Tfinal,Delta_t)
Nx = len(x_mesh)
Nt = len(t_mesh)

# Initial conditions for n
height_n = N0/Delta_x
n = np.zeros(len(x_mesh))
n[i_0] = height_n

solution_n = np.zeros([Nt,Nx])
solution_n[0,:] = n

# Generate matrix of the coefficients
upper_diagonal = np.diag(np.ones(Nx-1),1)
diagonal = np.diag(np.ones(Nx))
lower_diagonal = np.diag(np.ones(Nx-1),-1)

diff = D
advect = -q

cn = diff*Delta_t/(Delta_x**2)

An = - cn*upper_diagonal + (2*cn + 1)*diagonal - cn*lower_diagonal #implicit diffusion scheme
# more deleterious mutations
c_mu_n_centered = advect*Delta_t/(2*Delta_x)
An += c_mu_n_centered*upper_diagonal - c_mu_n_centered*lower_diagonal  # adds biased mutation

# LEFT BOUNDARY
# Dirichlet boundary condition - on the coefficient matrices
An[0,0] = 1.0
An[0,1:] = 0.0

# RIGHT BOUNDARY
# zero-flux - on the coefficient matrices
#An[Nx-1,:] = 0.0
An[Nx-1,Nx-1] = diff/Delta_x
An[Nx-1,Nx-2] = -diff/Delta_x
# if there are more deleterious mutations
An[Nx-1,Nx-1] += -advect

plt.figure(1,figsize=(8,6))
plt.plot(x_mesh, solution_n[0,:])#,'x-', lw=1)
my_xticks = ['$L_d$','0','$L_d$']
plt.xticks([x_mesh[0],0,x_mesh[-1]], my_xticks)
#plt.yticks([])
plt.tick_params(axis='both', which='major', labelsize=30)
plt.savefig(filename + "-ic.png")

# save parameters in text
x_mesh_positive = np.zeros(np.shape(x_mesh))
x_mesh_negative = np.zeros(np.shape(x_mesh))
for i in range(len(x_mesh)):
    if x_mesh[i]<0:
        # x_mesh_negative[i] = x_mesh[i]
        #x_mesh_negative[i] = np.arctan(x_mesh[i])
        x_mesh_negative[i] = np.log(x_mesh[i]-Ld) - np.log(-Ld)
    else:
        # x_mesh_positive[i] = x_mesh[i]
        #x_mesh_positive[i] = np.arctan(x_mesh[i])
        x_mesh_negative[i] = np.log(x_mesh[i]-Ld) - np.log(-Ld)


diag_x_negative = np.diag(x_mesh_negative)
diag_x_negative[0,0] = 0 # so only bc remains here
diag_x_negative[Nx-1,Nx-1] = 0    # so only bc remains here

diag_x_positive = np.diag(x_mesh_positive)
diag_x_positive[0,0] = 0 # so only bc remains here
diag_x_positive[Nx-1,Nx-1] = 0    # so only bc remains here

stop_simulation_index = Nt-1
# almost all implicit!
print('Final time = ', Nt*Delta_t)
    
n = solution_n[0,:]
integral_n = integral(solution_n[0,:], Delta_x)
population_integrals = np.zeros(Nt)
population_integrals[0] = integral_n
t = 0
stopping_flag = 0
while t < Nt - 1:
    
    t += 1
    
    print(t*Delta_t, end ='\r')
    
    previous_integral_n = integral_n

    ## CARRYING CAPACITY
    if growth_indicator == 'CC':
        integral_n = integral(n,Delta_x)
        integral_over_K = integral_n/K
        logistic = 1-integral_over_K
    ## PURELY EXPONENTIAL
    else:
        integral_n = integral(n,Delta_x)
        integral_over_K = 1
        logistic = 1
        
    population_integrals[t] = integral_n

    if logistic >= 0 and logistic <= 1:

        bn = bc_homo_L(bc_homo_R(n))    
        n_new = thomas(An - R*Delta_t*(integral_over_K*diag_x_negative + logistic*diag_x_positive), bn) 

        n = n_new

        solution_n[t,:] = n_new
    else:

        stop_simulation_index = t

        break
        
    if stopping_flag == 0 and (integral_n-previous_integral_n)/Delta_x < 1e-8 and integral_n > N0 and t < Nt/2:
        Nt = 2*t
        stopping_flag = 1
    
    
    
print('')

if stopping_flag == 1:
    population_integrals = population_integrals[:Nt]
    solution_n = solution_n[:Nt,:]
    t_mesh = t_mesh[:Nt]

f = open(filename + "-integraln.txt", "w")
for value in population_integrals[:-1]: 
    f.write(str(value) + '\n')
    f.write('\n')
f.close()

ff = open(filename + "-lastn.txt", "w")
for value in solution_n[-1,:]: 
    ff.write(str(value) + '\n')
    ff.write('\n')
ff.close()

fff = open(filename + "-xmesh.txt", "w")
for value in x_mesh: 
    fff.write(str(value) + '\n')
    fff.write('\n')
fff.close()

ffff = open(filename + "-tmesh.txt", "w")
for value in t_mesh: 
    ffff.write(str(value) + '\n')
    ffff.write('\n')
ffff.close()

            
print('filename = ' + filename)


plt.figure(figsize=(14,8))
plt.plot(t_mesh, np.ones(len(population_integrals)), 'k--', linewidth = 0.5)
plt.plot(t_mesh,population_integrals)
plt.ylim([0,10])
plt.savefig(filename + "-integral_zoomed.png")


plt.figure(figsize=(14,8))
plt.plot(t_mesh, np.ones(len(population_integrals)), 'k--', linewidth = 0.5)
plt.plot(t_mesh,population_integrals)
plt.savefig(filename + "-integral.png")