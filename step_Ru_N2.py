from ase.build import fcc111, root_surface, hcp0001, hcp0001_root, cut, stack, root_surface_analysis,add_adsorbate,surface,bulk
from ase.calculators.vasp import Vasp
from ase.io import write
from ase.io import read
from ase.spacegroup import crystal
from ase.visualize import view
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.hexagonal import HexagonalClosedPacked
import numpy as np
from ase import Atoms
Rubulk = bulk('Ru', 'hcp', covera = 1.584, orthorhombic = True)
cell  = Rubulk.cell
xscale = 6
yscale = 2
zscale = 2
bulk1 = Rubulk*(xscale,yscale,zscale)
surf = surface(Rubulk,(1,1,5),12)
surf.center(vacuum=5, axis=2)
surf = surf*(3,3,1)
#cell[0,0] = 6.75
#cell[2,2] = 6.4152
for atom in bulk1:
    if atom.position[2] < 0.1:
        if atom.position[0] >= xscale*cell[0,0]/2:
                atom.position[2] = atom.position[2]+zscale*cell[2,2]
d = cell[2,2]
L = cell[0,0]*xscale
Lz = d*zscale
alpha = np.arctan(d/L)

for atom in bulk1:
    x1 = atom.position[0]*np.cos(alpha)+atom.position[2]*np.sin(alpha)
    z1 = -atom.position[0]*np.sin(alpha)+atom.position[2]*np.cos(alpha)
    atom.position[0] = x1
    atom.position[2] = z1


bulk1.cell[0,0] = L/np.cos(alpha)#+(Lz-L*np.tan(alpha))*np.sin(alpha) 
#bulk1.cell[2,2] = Lz/np.cos(alpha)+(L-Lz*np.tan(alpha))*np.sin(alpha)
bulk2 = bulk1*(1,1,1)




alphasurf = -np.arctan(37.3471/60.75)
del surf[[atom.index for atom in surf if atom.position[1] > 11.8]]
A = np.array([[np.cos(alphasurf),np.sin(alphasurf),0],[-np.sin(alphasurf), np.cos(alphasurf),0],[0,0,1]])
for atom in surf:
    x1 = atom.position[0]*np.cos(alphasurf)+atom.position[1]*np.sin(alphasurf)
    y1 = -atom.position[0]*np.sin(alphasurf)+atom.position[1]*np.cos(alphasurf)
    atom.position[0] = x1
    atom.position[1] = y1
surf.cell[0,:] = A.dot(surf.cell[0,:])
surf.cell[1,:] = A.dot(surf.cell[1,:]*(1/3))
del surf[[atom.index for atom in surf if atom.position[1] >5.65 ]]
surf.cell[0,:] = surf.cell[0,:]*(2/3)


# atom in surf: # CONVERTING TO DIFFERENT STEP TYPE
#    atom.position[0] = atom.position[0] - surf.cell[1,0]/2
#    if atom.position[1] < atom.position[0]*5.65613/9.20045:
#        atom.position[0] = atom.position[0] + surf.cell[1,0]


surf = surf*(1,1,1) # CHECKING BOUNDARY CONDITIONS

sx = -8.04622
sy = 1.48844
sz = 13.8363

sx1 = -10.7331
sy1 = 1.98911
sz1 = 17.2297

sx2 = -7.27966667 
sy2 = 1.53134
sz2 = 13.7023

sx3 = -9.58103
sy3 = 2.696125
sz3 = 17.2297

sx4 = -9.1977 
sy4 = 4.65019667
sz4 = 12.76533333

h = 2
rN = 1.12

surf.cell[2,:] = surf.cell[2,:]*1.2
surf1 = surf.copy()
surf2 = surf.copy()
surf3 = surf.copy()
surf4 = surf.copy()
surf6 = surf.copy()
surf7 = surf.copy()
mol_ontop = Atoms('2N', positions =[(sx,sy,sz+h),(sx,sy,sz+rN+h)])
mol1 = Atoms('2N', positions =[(sx,sy,sz+h+5),(sx,sy,sz+rN+h+5)])
atom1 = Atoms('N', positions = [(sx1,sy1,sz1+h)])
atom2 = Atoms('N', positions = [(sx2+0.18,sy2-0.4,sz2+0.9)])
atom3 = Atoms('N', positions = [(sx3,sy3,sz3+1.5)])
atom4 = Atoms('N', positions = [(sx4+0.2,sy4-0.2,sz4+1)])
surf1.extend(mol_ontop)
surf2.extend(mol1)
surf3.extend(atom1)
surf4.extend(atom2)
surf6.extend(atom3)
surf7.extend(atom2)
surf7.extend(atom4)
write('step_ontop',surf1,format = 'vasp')
write('step_before',surf2,format = 'vasp')
write('step_ads1',surf3,format = 'vasp')
write('step_ads2',surf4,format = 'vasp')
write('step_ads_f',surf7,format = 'vasp')
write('step_ads_3',surf6,format = 'vasp')
write('surf', surf,format = 'vasp')
write('bulk1',bulk2, format = 'vasp')
calc = Vasp(prec = 'Normal',kpts = (5,5,1), gamma = True, xc = "pbe", lreal = False)
bulk1.set_calculator(calc)
bulk1.get_potential_energy()
