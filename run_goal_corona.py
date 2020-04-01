from htmd.ui import *
from os.path import join
from moleculekit.projections.metricdistance import MetricDistance

#HID 42 and name CA #binding site
#MOL ligand

os.makedirs('./corona/goal20', exist_ok=True)
shutil.copytree('./prod/', './corona/goal20/generators')
os.chdir('./corona/goal20')

queue = SlurmQueue()
queue.jobname = 'covid19_20'
queue.partition = 'normalGPU'
queue.datadir = './data'
queue.envvars =  'ACEMD_HOME,HTMD_LICENSE_FILE,PATH'

#some kinases do not go into dfg-out
#https://www.sciencedirect.com/science/article/pii/S1074552113002147
def mygoalfunction(mol):
    distance_metric = MetricDistance('protein and resname HID and resid 42 and name CA', 'resname MOL and name C11')
    distance = distance_metric.project(mol)
    distance[distance < 20.0] = 1.0
    print('THE PROJECTION VALUES:',distance)
    return -distance  # or even 1/distance

adg = AdaptiveGoal()
adg.app = queue
adg.nmin = 10
adg.nmax = 20
adg.nepochs = 999
adg.generatorspath = './generators'
adg.updateperiod = 120  # execute every 2 minutes
adg.nosampledc = True
adg.goalfunction = mygoalfunction
adg.run()

