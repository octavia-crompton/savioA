# Save a dictionary into a pickle file.
import pickle
import os, sys
from commands import getoutput as cmd
import multiprocessing as mp

params = { 'ncol' : 10, 'nrow' : 50, 
        'dx' : 1., 'epsh' : .0002,
        'dt_sw'  : .05, 'tmax'  : 60*60., 
        'tr' : 20*60., 'rain' : 0.036/3600.,
        'ifricV': 1, 'ifricB': 1,
        'xnv' : 0.09, 'xni' : 0.02,
        'pveg' : .3, 'Dveg' : 8e-3,
        'itype1' : 0, 'itype3' : 1, 
        'itype2' : 1, 'itype4' : 1,
        'beta' : 1.0,       
        'slope' : 0.01, 'logscale' : 10,     
        'seed' : 5, 
        'append' : '',
        'iscouple': 0,
        'stop_tol'  : 1e-2, 'iscale' : 100,
        'htop'  : 0, 'infL'  : 2,
        'ponding' : 1, 'dz' : .4,    # [cm]
        'zmin' : 0. , 'zmax' : 20.,  # [cm]
        'hinit' : -100, 'depth'  : 0., # cm
        'ksatB' : 0.2/3600, 'ksatV' : 3./3600, 
        'zs' : 20, 
        'level' : 'B2_plane_randv_Sf',
        'weeds' : 0, 'runmodel' : 1,
          } 
          
folder = params['level']
print folder

cores = mp.cpu_count()
print 'There are %s cores on this machine, \n model runs will be performed on each core'%(str(cores))


cmd('mkdir {0}/'.format(folder))
cmd('mkdir {0}/input'.format(folder))
cmd('mkdir {0}/output'.format(folder))

fname = '{0}/input/params.p'.format(folder)
pickle.dump( params, open(fname, "wb" ) )

if os.getcwd().split('/')[1] == 'global':
  print 'savio run'
  print 'python build_model.py {0} 0'.format(folder)
  a = cmd('python build_model.py {0} 0'.format(folder))

else:
  print 'local run'
  print 'python build_model.py {0} 1'.format(folder)
  a = cmd('python build_model.py {0} 1'.format(folder))
  print a