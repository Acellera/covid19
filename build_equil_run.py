from htmd.ui import *
from os.path import join
from moleculekit.projections.metricdistance import MetricDistance
#Protein
prot = Molecule('noncovhits/Mpro-x0104.pdb')
prot.reps.add(sel='resname LIG', style='CPK')
prot.reps.add(sel='not protein', style='Licorice')
prot.reps.add(sel='protein', style='Lines')
prot.view(name='original')

prot.remove('resname LIG')
prot.filter('protein')
prot = proteinPrepare(prot, pH=7.0)


prot = autoSegment(prot, sel='protein')
prot.set('segid', 'W', sel='water')
prot.set('segid', 'CA', sel='resname CA')
prot.center()

#LIGAND
ligand = Molecule('/home/alejandro/covid19/MDrun/melatonin_cov/parameters/GAFF2/mol.mol2')
ligand.center()
ligand.set('segid','L')
ligand.set('resname','MOL')
ligand.rotateBy(uniformRandomRotation())

from moleculekit.util import maxDistance
D = maxDistance(prot, 'all')
D += 5
ligand.moveBy([D, 0, 0])

mol = Molecule(name='combo')
mol.append(prot)
mol.append(ligand)
mol.reps.add(sel='protein', style='NewCartoon', color='Secondary Structure')
mol.reps.add(sel='resname MOL', style='Licorice')
mol.view()


# We solvate with a larger box to fully solvate the ligand
DW = D + 3
smol = solvate(mol, minmax=[[-DW, -DW, -DW], [DW, DW, DW]])
smol.reps.add(sel='water', style='Lines')
smol.view()

topos_amber = amber.defaultTopo()
frcmods_amber = amber.defaultParam() + [join('melatonin_cov/parameters/GAFF2/', 'mol.frcmod')]
bmol_amber = amber.build(smol, topo=topos_amber, param=frcmods_amber, outdir='./build_amber')
bmol_amber.view(name='ready2run')

from htmd.protocols.equilibration_v2 import Equilibration
md = Equilibration()
md.runtime = 3
md.timeunits = 'ns'
md.temperature = 310
md.useconstantratio = False
md.write('./build_amber/', './equil')

local = LocalGPUQueue()
local.submit('./equil/')
local.wait()


from htmd.protocols.production_v6 import Production
md = Production()
md.runtime = 10
md.timeunits = 'ns'
md.temperature  = 310
md.acemd.bincoordinates = 'output.coor'
md.acemd.extendedsystem  = 'output.xsc'
md.write('equil','prod')

local = LocalGPUQueue()
local.submit('./prod/')
local.wait()




