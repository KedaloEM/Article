import numpy as np
import random
import os.path

#constants
mN = 14*1.660539*10**(-27)
mRu = 101*1.660539*10**(-27)
w = 2774 # cm-1
c = 3*10**(10) # SOL
h = 6.626070*10**(-34) # Plank constant
k = 1.380648*10**(-23) # Boltzmann constant
d = 1.12 # N-N equilibrium distance
vel_scale = 10**(-5) # m/s -> A/fs
delta_x_max = 10**(10)*np.sqrt(4*h/(mN*w*c))/(2*np.pi) # maximum stretch during the vibrational period

# input of initial parameters:
print('Enter the initial height(A):')
#l = float(input())
l = 2.2
print('Enter the orientation angles:')
#alpha,theta = [float(x) for x in input().split()]
#theta = theta*np.pi/180
alpha = 0
alpha = alpha*np.pi/180
theta = 0
theta = theta*np.pi/180
print('Enter the initial phase (in pi units) and direction (+-1)')
#phase, dir = [float(x) for x in input().split()]
#phase = phase*np.pi
phase = 0
dir = 1
if (dir!=1) and (dir!=-1):
    raise ValueError ('Direction should be +-1!')
print('Enter the translational temperature(K):')
T = 700
#T = float(input()) # translational temperature

# Positions and velocities
r = d+delta_x_max*np.cos(phase) # N-N initial distance
vmod = np.sqrt(h*c*w/mN - ((c*w*2*np.pi*10**(-10)*(r-d))**2)/4)
print(mN*c**2*w**2*np.pi**2*10**(-20)*(r-d)**2)
sigma = np.sqrt(k*T/(2*mN))
sigma_Ru = np.sqrt(k*T/(mRu))
vtransz = np.abs(np.random.normal(0,sigma))
vtransx = np.random.normal(0,sigma)
vtransy = np.random.normal(0,sigma) # translational velocity
#vtransx = -364.99422#1020.27798#438.01235
#vtransy = 969.81114#597.49622#-187.450715 
#vtransz = 273.1566#794.847364#472.74629
x_0 = -6.9 # above the slab center
y_0 = 2.2
z_0 = 13.8+l
r1 = [x_0 + r*np.cos(theta)*np.sin(alpha)/2, y_0 + r*np.sin(alpha)*np.sin(theta)/2, z_0 + r*np.cos(alpha)/2]
r2 = [x_0 - r*np.cos(theta)*np.sin(alpha)/2, y_0 - r*np.sin(alpha)*np.sin(theta)/2, z_0 - r*np.cos(alpha)/2]
v1 = [vel_scale*(dir*vmod*np.cos(theta)*np.sin(alpha)+vtransx), vel_scale*(dir*vmod*np.sin(alpha)*np.sin(theta)+vtransy), vel_scale*(dir*vmod*np.cos(alpha)-vtransz)]
v2 = [vel_scale*((-dir)*vmod*np.cos(theta)*np.sin(alpha)+vtransx), vel_scale*((-dir)*vmod*np.sin(alpha)*np.sin(theta)+vtransy), vel_scale*(-dir*vmod*np.cos(alpha)-vtransz)]
sign1 = False
step = 0
Nat = 96+2 #number of atoms 
#writing into new POSCAR file
for i in range(21):
    path = '/home/kedaloem/VASP/POSCAR_%d' % i
    if not os.path.isfile(path):
        output = open('/home/kedaloem/VASP/POSCAR_%d'%i, 'w')
        break
input = open('/home/kedaloem/VASP/STEP','r')
for line in input:
    if sign1:
        step += 1
        if step <= 9:
            output.write(line[:len(line)-2]+" F F F" + "\n")
        if step<=Nat-2 and step > 9:
            output.write(line[:len(line)-2]+" T T T"+"\n")
        elif step==Nat-1:
            output.write("  "+"  ".join([str(q) for q in r1]) + " T T T" + "\n")
        elif step==Nat:
            output.write("  "+"  ".join([str(q) for q in r2]) + " T T T" + "\n")
    else:
        if line != 'Cartesian\n':
            output.write(line)
    if line == 'Cartesian\n':
        output.write('select\n')
        output.write(line)
        sign1 = True

output.write('Cartesian\n')
for i in range(Nat-2):
    Ru_vx = np.random.normal(0,sigma_Ru)
    Ru_vy = np.random.normal(0,sigma_Ru)
    Ru_vz = np.random.normal(0,sigma_Ru)
    Ru_v = vel_scale*np.array([Ru_vx,Ru_vy,Ru_vz])
   # output.write("  "+"0  0  0"+"\n") # for frozen surface
    output.write("  "+"  ".join([str(q) for q in Ru_v]) + "\n")
output.write("  "+"  ".join([str(q) for q in v1]) + "\n")
output.write("  "+"  ".join([str(q) for q in v2]) + "\n")
input.close()
output.close()

