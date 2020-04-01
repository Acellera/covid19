from moleculekit.molecule import Molecule
from moleculekit.projections.metricrmsd import MetricRmsd
from moleculekit.projections.metricdistance import MetricSelfDistance
from moleculekit.tools.sequencestructuralalignment import sequenceStructureAlignment
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from jobqueues.localqueue import LocalGPUQueue
from jobqueues.slurmqueue import SlurmQueue
from htmd.adaptive.adaptivegoal import AdaptiveGoal
from htmd.config import config
from htmd.ui import *

#for reference
fragalysis = Molecule('/home/alejandro/covid19/noncovhits/Mpro-x0104.pdb')


sims = simlist(glob('/shared/alejandro/corona/goal20/filtered/*/'), glob('/shared/alejandro/corona/goal20/filtered/filtered.pdb'))
metr = Metric(sims)
metr.set(MetricDistance('protein and noh' , 'resname MOL and noh', groupsel1='residue', groupsel2='residue', metric='contacts', threshold=5))
data = metr.project()
data.fstep = 0.1

data.plotTrajSizes()
data.dropTraj()

tica = TICA(data, 1, units='ns')
dataTica = tica.project(3)

dataBoot = dataTica.bootstrap(0.8)
dataBoot.cluster(MiniBatchKMeans(n_clusters=1000))
model = Model(dataBoot)
model.plotTimescales()

model.markovModel(2, 5, units='ns')
model.viewStates(ligand='resname MOL and noh')

abc=model.getStates()
prot_ref = abc[0].copy()
prot_ref.filter('protein')
frag_alig = sequenceStructureAlignment(fragalysis,prot_ref,molseg='', refseg='')[0]
frag_alig.view(sel='resname LIG', style='Licorice')


