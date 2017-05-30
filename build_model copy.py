import os
import sys
from os.path import dirname
import pickle
from commands import getoutput as cmd
#import pandas as pd
import shutil

parent_dir = dirname(dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'/modelfunctions'))

sys.path.append('modelfunctions/' )

mymodules = ['analytic_fxns',
             'output_dry', 
             'input_dry',              
             'misc_fxns', 
             'input_phi',
             'GWfunctions' ]
for mymod in mymodules:
    if mymod in sys.modules:  
        del sys.modules[mymod]
        

from GWfunctions import *
from analytic_fxns import *
from output_dry import *
from misc_fxns import *
from input_dry import *
from input_phi import *
from richards_functions import *

    
def main(argv):
    inputs = {}
    outputs = {}

    level =  sys.argv[1]  
    localrun =  sys.argv[2]      
    folder = level
    fname = '{0}/input/params.p'.format(folder)
    params = pickle.load( open( fname, "rb" ) )

    #fname = '{0}/input/rparams.p'.format(folder)
    #rparams = pickle.load( open( fname, "rb" ) )
        
    for key,val in params.items():
            exec(key + '=val')
                      
    if level[0] == 'R':
        case = 'richards'
        iscouple = 1
        ifixh  = 0
        slope = 0.00
        
    elif 'fixh' in level:
        case = 'fixed_h'
        icouple = 1
        ifixh  = 1      
        slope = 0.00        
    else:
        case = 'SVE'
        ifixh  = 0
        if level[2] == 'i':
            iscouple = 1
        else:
            iscouple = 0 
      

    topo = level.strip().split('_')[1]
    veg = level.strip().split('_')[2]
    if level[0] != 'A':  
        varcase =  level.strip().split('_')[3]
    else: 
        varcase = 'none' 
    
    print 'topography = {0}; veg = {1} '.format(topo, veg, varcase)
    print '    f(veg) = {0} '.format(varcase)    
    try:
        dim = int(level[1])
    except:
        dim = 1
        print 'Defaulting to 1D grid'.format(ncol)
    if dim == 2 and ncol == 2:
        ncol = 10          
        print '2D! Change ncol from 2 to 10'.format(ncol)
    

    #for key,val in rparams.items():
    #        exec(key + '=val')
    z =  np.arange(zmin, zmax+dz, dz)        
    nz = z.shape[0]            
    hinit = hinit*np.ones(nz)
    
    #Time parameters for python Richards solver
    dt  = dt_sw*iscale
    t  = np.arange(0, tmax + dt, dt)
    prate  = np.ones_like(t)*rain*100.
    prate[np.where(t>= tr)[0][0]:] = 0
    t  = t[1:]
    nt  = t.shape[0] 
    
    p = rain*3600
    
    nprt = int(60./dt_sw)
    print nprt
    if nprt < 1:
        nprt = 1
 
    Lx  = dx*ncol
    Ly  = dx*nrow
    
    make_folders(level)
    nt_sw = write_prate(folder = folder, tmax = tmax, dt_sw = dt_sw, tr= tr,
                        rain = rain)
    inputs['nt_sw'] = nt_sw
    xcc,ycc,zcc, xncc, isvegcc = build_coords(folder = folder, slope = slope, 
                                ncol = ncol, nrow = nrow, dx = dx, 
                               level = level, topo = topo, veg = veg,
                                seed = seed)    

    write_param(folder = folder,  dt_sw = dt_sw, 
                 tmax = tmax,  nprt = nprt, epsh = epsh,iscouple = iscouple,       
                 iscale= iscale, ifixh = ifixh, stop_tol = stop_tol, depth = depth,
                 beta = beta, ifricV = 1, ifricB = 1, 
                 xnv = xnv, xni = xni, pveg = pveg, Dveg = Dveg, tr = tr)                         
    inum, ipos, itype = write_dryin(folder = folder, ncol = ncol, nrow = nrow,
                                itype3=itype3, itype1 = itype1,
                                    itype2=itype2, itype4 = itype4) 
                                    
   ## richards input
    phi_veg, phi_bare  = make_phi( dz = dz, zmin = zmin, zmax = zmax, z = z,
              zs = 20, infL = infL, ksatB = ksatB, ksatV = ksatV)

    if level[0]  == 'B'and 'infl' not in varcase:
        phi_veg = phi_bare

    write_phi( phi_veg ,phi_bare , hinit,folder = folder, fname = 'vanG' , 
              infL = infL, dz = dz, nz = nz, zmax = zmax)
    
    DeltaPlus, DeltaMinus, MPlus, MMinus = Richards_matrices(nz, infL = infL)

    
    globals().update(locals())                      
    
    make_titlestr()     
    
    write_inc(folder = folder, nt = nt_sw, nz = nz,  ncol = ncol, nrow = nrow)                                      

    if runmodel == 1:
      
      shutil.copyfile('dry10.for', '{0}/dry.for'.format(folder))             

      if localrun == '1':
        print 'local run'
        a = cmd("gfortran -o {0}/sw  -framework accelerate {0}/dry.for".format(folder))
        print a
        a = cmd("cd {0} \n ./sw  \n cd ..".format(folder))
      else:
        print 'savio run'        
        a = cmd("module load gcc mkl \n gfortran  -lmkl_gf_lp64 -lmkl_core \
                -lmkl_sequential -lpthread -lm -ldl -o {0}/sw {0}/dry.for".format(folder))
        a = cmd("module load gcc mkl \n cd {0} \n ./sw  \n cd ..".format(folder))        
    
    read_output(folder)
    write_summary(folder, datastr)
    
    
    varsin = [ 'xcc', 'ycc', 'zcc', 'isvegcc', 'dx', 'iscouple',\
              'case', 'Lx', 'Ly', 'dim',\
               't', 'dt', 'prate', \
               'nt', 'nt_sw', \
              'titlestr', 'datastr',\
              'tr', 'slope', 'dt_sw', \
              # 'p', 'rain', 'epsh', 'beta',
              ]
    varsout =  [ 't_p', 'hydro_p', \
               'h', 'u', 'v', 'vmag', \
               'h_tr', 'u_tr', 'v_tr', \
               'zinflmap', 'eqlind', \
               ]
               
    weedvars =  [ 't_sw', 'cfl_sw', 'ftime_sw',
                'fluxin_for', 'fluxout_for', 'newmass_for', 
               'error_for2', 'error_for',               
               'dvol', 'flux', 'zinfl', 'rain_sw', 'total', \
               'flux1', 'flux2', 'flux3', 'flux4',
                ]
               
    allvars =  varsin + varsout # + weedvars
    
    outputs = dict([(name,globals()[name]) for name in allvars])
    
    outweeds = dict([(name,globals()[name]) for name in weedvars])
    
    
      
    if slope > 0 and iscouple == 1:

        phi = phi_bare
        for key in phi_veg.keys():
            phi[key] = phi_veg[key]*isvegcc.mean() + phi_bare[key]*(1 - isvegcc.mean())
        xn = xnv*isvegcc.mean() + xni*(1 - isvegcc.mean())

        alpha = np.sqrt(slope)/xn # [m^(1/3)/s]  sqrt(So)/n
        [Cinit,Kinit,thetainit] = vanGenuchten(hinit,phi) 
        Ksat=phi['ksat'][-1]/100.
        Ao = np.sqrt(2*(1-1/phi['n'][-1])**(4/3.)*Kinit[-1]*
             (phi['theta_S'][-1] -
              phi['theta_R'][-1])/phi['alpha'][-1])/100.
        try:
             [t_GW,q_GW,Perc_runoff] = Comparison_function2(rain,Ly,Ksat,Ao,
                                           tr,alpha)
             q_GW = q_GW/Ly*100.                                
             outputs['t_GW'] = t_GW
             outputs['q_GW'] = q_GW
        except:
          print 'GW failed. Check input params'
          
    elif slope > 0 and iscouple == 0:
      
        xn = isvegcc.mean()*xnv +(1- isvegcc.mean())*xni
        q0_p, tc = analytical_hydro(nrow = nrow, slope = slope, xn = xn, 
                            tr = tr, t = t_p, dx = dx, rain = rain)
    
        outputs['q0_p'] = q0_p
    
    
    fname = 'outputs/main/{1}.pkl'.format(folder, datastr)
    pickle.dump( outputs, open(fname, "wb" ) )
                           
    fname = 'outputs/weeds/{1}.pkl'.format(folder, datastr)
    pickle.dump( outweeds, open(fname, "wb" ) )
      
              
def make_folders(level):
    folder  =  level
    cmd('mkdir {0}/'.format(folder))
    cmd('mkdir {0}/input'.format(folder))
    cmd('mkdir {0}/output'.format(folder))
    cmd('mkdir {0}/summary'.format(folder))
    cmd('mkdir {0}/weeds'.format(folder))
        

    if os.path.isdir("outputs/main") == False:       
      cmd('mkdir outputs'.format(folder))            
      cmd('mkdir outputs/main'.format(folder))                
      cmd('mkdir outputs/weeds'.format(folder))                
      cmd('mkdir outputs/summaries'.format(folder))
      cmd('mkdir outputs/Csummaries'.format(folder))    

def read_output(folder):
  
    summarize_exit()
    t_p = read_time(folder = folder)
    t_sw, cfl_sw, ftime_sw = read_cfl(folder = folder)

    tdum, hydro_p = new_hydro(folder = folder, fname = 'hydro',Lx=Lx,Ly=Ly)  
    h,u,v = get_h(folder =  folder, ncol = ncol, nrow = nrow)
    vmag = np.sqrt(u**2+v**2)
    get_tp(folder)
    try:
        zc = get_zc(folder = folder)
    except:
        zc = zcc
    tdum, dvol, flux, zinfl, rain_sw, total, = get_dvol(folder = folder, 
                                                    Lx = Lx, Ly = Ly)

    flux1, flux2, flux3, flux4 = lateral_fluxes(folder)
 
    fluxin_map, fluxout_map, newmass_map, error_map = read_fluxmap(
                folder=folder, ncol=ncol, nrow=nrow, dt_sw=dt_sw)

    fluxin_for = fluxin_map.mean(1).mean(1)
    fluxout_for = fluxout_map.mean(1).mean(1)
    newmass_for = newmass_map.mean(1).mean(1)
    error_for = error_map.mean(1).mean(1)
    error_for2 = fluxin_for - newmass_for  + fluxout_for 
    
    zinflmap = read_zinflmap(folder, nrow, ncol, dt_sw, dx)        
    eqlind = int(np.where(t_p < tr)[0][-1])
    fluxin_11, fluxout_11, newmass_11, error_11 = read_soil_fluxes11(folder, dt)
    
    h_tr = h[eqlind]
    u_tr = u[eqlind]
    v_tr = v[eqlind]  
    zinflmap = zinflmap.sum(0)  
    globals().update(locals())      
    

def write_summary(folder, datastr):

    # H ending to represent height in cm
    volH = np.mean(h[-1])*100
    fluxH = -np.cumsum(flux)[-1]  # convert m^3 to m
    rainH = rain*tr*100   # units = cm

    f = open('{0}/summary{1}.txt'.format(subfolder(), append), 'w')
    g = open('{0}/Csummary.txt'.format(subfolder()), 'w')

    f.write('{0} \n '.format(level))

    f.write('\n ************ Results summary **************** \n')
    if slope > 0:
        f.write('\n ---  hydrograph summary  --- \n')
        f.write('  hydrograph peak  = {0:.3f}cm/hr\n'.format(hydro_p.max()*3600.))
        g.write('  hydro_peak  = {0}\n'.format(hydro_p.max()*3600.))

        time2peak =  t_p[hydro_p == hydro_p.max()][0]/60.
        f.write('  time to peak  = {0:.3f}min\n'.format(time2peak))
        g.write('  time2peak = {0}\n'.format(time2peak))

        time2peak90 = t_p[hydro_p >= hydro_p.max()*.9][0]/60.
        f.write('  time to 90% peak  = {0:.3f}min\n'.format(time2peak90))
        g.write('  time2peak90  = {0}\n'.format(time2peak90))

    if iscouple == 1:
        inflH = - np.cumsum(zinfl)[-1]  # units = cm
        f.write(' infiltration ratio  = {0:.2f}\n'.format( inflH/rainH))
        g.write(' infl_percent  = {0:.2f}\n'.format( inflH/rainH*100))

    if iscouple == 1 and isvegcc.std() > 0 :
        f.write('\n --- Veg/interspace infl partitioning --- \n')
        percent_veg = float(np.sum(isvegcc > 0.))/np.sum(isvegcc >=  0)
        f.write(' vegetation area  = {0:.2f}\n'.format( percent_veg))
        g.write(' veg_area = {0:.3f}\n'.format(percent_veg))

        f.write(' average infiltration depths \n')
        f.write('       over vegetation = {0:.3f}cm\n'.format(
                zinflmap[isvegcc > 0].mean()))


        f.write('       over bare ground = {0:.3f}cm\n'.format(
                zinflmap[isvegcc == 0].mean()))
        f.write(' infl over veg / total infl = {0:.3f}\n'.format(
                zinflmap[isvegcc >0].sum()/np.sum(zinflmap)))

        g.write(' veg_infl_percent = {0:.3f}\n'.format(
                zinflmap[isvegcc >0].sum()/np.sum(zinflmap)))

        infl_veg = zinflmap[isvegcc >0].sum()/np.sum(zinflmap)

        if tp_for:
            f.write('\n 1st ponding  = {0}s\n'.format( tp_for))
        #g.write('tp = {0}'.format(tp_for))
        if t_sat:
            f.write(' 1st soil saturation  = {0:.2}min\n'.format( t_sat/60.))
        #g.write('t_sat = {0}'.format(t_sat))
        if t_unsat:
            f.write(' soil unsaturation  = {0:.2}min\n'.format( t_unsat/60.))
        if t_final:

            f.write(' no more water  = {0:.2f}min\n'.format( t_final/60.))
            g.write(' t_final = {0}\n'.format(t_final))

    vmax = np.max(vmag)*100.
    hmax = np.max(h)*100.
    Remax = np.max(h*vmag/1e-6)
    f.write(' max v  = {0:.3f}cm/s\n'.format(vmax*100))
    g.write(' vmax  = {0}\n'.format(vmax))
    f.write(' max h  = {0:.3f}cm\n'.format(hmax*100))
    g.write(' hmax  = {0}\n'.format(hmax))
    f.write(' max Re  = {0:.3f}m2/s\n'.format(Remax))
    g.write(' Remax  = {0}\n'.format(Remax))


    #################### Output Summary  ####################
    f.write('\n --- Fortran run-time ---\n') # total infiltration from sw vs. richards time grids?
    f.write(' time = {0:.3}min \n'.format((ftime_sw)[-1]/60.))
    g.write(' runtime = {0:.3} \n'.format((ftime_sw)[-1]/60.))
    f.write('\n************ Error summary ******************\n')

    if iscouple == 1  :
        f.write('\n --- richards errors ---\n')
        #f.write('  richards solver error = {0:.2}cm\n'.format( np.cumsum(error_for)[-1]))
        richards_err =  np.cumsum((error_for2)[1:])[-1]  # in cm
        richards_infl =  np.cumsum((fluxin_for)[1:])[-1]  # in cm
        f.write(' richards error = {0:.3}cm\n'.format( richards_err) )
        g.write(' richards_err = {0:.3} \n'.format( richards_err) )
        f.write(' cumulative infl = {0:.3}cm \n'.format(richards_infl))
        g.write(' richards_infl = {0:.3} \n'.format(richards_infl))
        f.write(' error relative to surface flux = {0:.1f}%\n '.format(
            np.cumsum(error_for)[-1]/np.cumsum(fluxin_for)[-1]*100))


    f.write('\n --- mass balance summary --- \n')
    f.write(' rain depth = {0:.4}cm \n'.format(rainH))
    f.write(' final depth =  {0:.4}cm \n'.format(volH))
    f.write(' total runoff =  {0:.4}cm\n'.format(fluxH))
    if iscouple == 1:
        inflH = -np.cumsum(zinfl)[-1]  # units = cm
        f.write(' avg infl depth = {0:.2f}cm\n'.format(inflH))
        mass_bal = rainH - volH - fluxH - inflH
    else:
        mass_bal = rainH - volH - fluxH
    f.write(' mass balance = {0:.2}\n'.format(mass_bal))


    f.write('\n --- boundary fluxes ---\n')
    f.write(' flux2           = {0:.3}cm \n'.format(np.sum(abs(flux2))*100./Lx/Ly))
    f.write(' flux3 (top)     = {0:.3}cm \n'.format(np.sum(abs(flux3))*100./Lx/Ly))
    f.write(' flux4           = {0:.3}cm \n'.format(np.sum(abs(flux4))*100./Lx/Ly))
    f.write(' backwards flux1 = {0:.3}cm \n'.format(np.sum(flux1[flux1<0])*100./Lx/Ly))


    if iscouple == 1:
        f.write('\n --- iscale infl error ---\n') # total infiltration from sw vs. richards time grids?
        inflH_sw = np.cumsum(fluxin_for)[-1]
        f.write(' SVE - richards infl discrepency =  {0:.2} cm\n'.format(
                inflH_sw - inflH))
        f.write('      --      percent discrepency = {0:.2}%\n'.format(
                (inflH_sw - inflH)/inflH))

    if level[0] == 'A':
        f.write('\n ************ Compare to analytical ************ \n')

        if slope > 0 and iscouple == 1:
            f.write('\n --- GW error ---\n')
            dt_GW = np.diff(np.hstack((0, t_GW)))
            dt_p = t_p[1]- t_p[0]
            fluxH1= np.cumsum(flux1*100/Lx/Ly)[-1]
            GW_error = np.cumsum(q_GW*dt_GW)[-1] - fluxH1

            f.write(' total model runoff =  {0:.3}cm\n'.format(fluxH1))
            f.write(' total GW runoff = {0:.3}cm\n'.format(np.cumsum(q_GW*dt_GW)[-1]))
            f.write(' GW - model runoff   = {0:.3}cm\n'.format(GW_error))

            f.write('\n total model infl depth = {0:.3}cm\n'.format(inflH))
            f.write(' total GW infl depth = {0:.3}cm\n'.format(rainH - np.cumsum(q_GW*dt_GW)[-1]))
            f.write(' GW - model infl depth = {0:.3}cm\n'.format(rainH - np.cumsum(q_GW*dt_GW)[-1]-inflH))

    f.write('\n ************ Input parameters **************** \n')
    f.write('\n --- physical parameters: topography --- \n')
    f.write(' slope   = {0}% \n'.format(slope*100))
    f.write(' Ly      = {0}m \n'.format(Ly))
    f.write(' Lx      = {0}m \n'.format(Lx))

    f.write('\n --- physical parameters: rain --- \n')
    f.write(' P       =  {0}cm/hr \n'.format(rain*3.6e5))
    f.write(' tr      =  {0}min \n'.format(tr/60.))
    f.write(' tmax    =  {0}min \n'.format(tmax/60.))

    if iscouple == 1:
        f.write('\n --- physical parameters: soil --- \n')
        f.write(' Ksat (seal) =  {0}cm/hr \n'.format(phi_bare['ksat'][-1]*3600))
        f.write(' Ksat (veg)  =  {0}cm/hr \n'.format(phi_veg['ksat'][-1]*3600))

    f.write('\n ---  roughness parameters  --- \n')
    if np.mean(isvegcc) == 0:
        f.write(' ' + ' = '.join(fBstr(ifricB).split('='))+ '\n')
    elif np.mean(isvegcc) == 1:
        f.write(' ' + ' = '.join(fVstr(ifricV).split('='))+ '\n')
    else:
        f.write(' ' +' = '.join(fBstr(ifricB).split('='))+ '\n')
        f.write(' ' + ' = '.join(fVstr(ifricV).split('='))+ '\n')


    f.write('\n --- computational parameters  --- \n')
    f.write(' ncol    =  {0} \n'.format(ncol))
    f.write(' nrow    =  {0} \n'.format(nrow))

    f.write(' dt_sw   =  {0} \n'.format(dt_sw))
    f.write(' dx      =  {0} \n'.format(dx))
    f.write(' epsh    =  {0}mm \n'.format(epsh*1000.))
    f.write(' beta    =  {0}mm \n'.format(beta))

    if iscouple == 1:
        #f.write('\n --- Richards params  --- \n')
        f.write(' iscale   = {0} \n'.format(iscale))
        f.write(' stop_tol = {0} \n'.format(stop_tol))


    f.close()
    g.close()
    shutil.copyfile('{0}/summary{1}.txt'.format(subfolder(), append),
                    'outputs/summaries/{1}.txt'.format( folder, datastr))


    shutil.copyfile('{0}/Csummary.txt'.format(subfolder(), append),
                    'outputs/Csummaries/{1}.txt'.format( folder, datastr))

    globals().update(locals())

 
def get_tp(folder): 
    tp_for = None
    t_sat = None      
    t_unsat = None     
    t_final = None
    for line in open('{0}/output/summary.out'.format(folder), 'r'):
        a = line.strip().split(' ')
        a =  [b for b in a if b] 

        try:
            if a[0] == '1st':
                tp_for = float(a[-2])
            if a[0] == 'saturated':            
                t_sat = float(a[-1])
            if a[0] == 'unsaturated':
                t_unsat = float(a[-1])
            if a[0] == 'final':                            
                t_final = float(a[-1])                          
                
        except IndexError:
            continue
    globals().update(locals())
        
def make_titlestr():

    titlestr  = '{0}: tr={1:.1f},p={2},slope={3},'.format(
        level, tr/60., rain*3.6e5,slope*100)

    #if dim == 2:
    titlestr = titlestr + 'Ly={0},Lx={1},'.format(Ly, Lx)
    #elif dim == 1:
    #    titlestr = titlestr + 'Ly={0},Lx={1},'.format(Ly)            

    #if topo == 'log':
    #    titlestr = titlestr + 'logC={0},'.format(logscale)

    titlestr = titlestr + "\n dx={0},dt_sw={1},epsh={2},beta={3},".format(
                                                    dx,dt_sw,epsh*1000,beta)
    if level[0] == 'A' and 'bare' not in level:
        titlestr = titlestr + fVstr(ifricV) + ','
    elif level[0] =='A' and 'bare' in level:
        titlestr = titlestr + fBstr(ifricB) + ','
    elif level[0] == 'B' and '_infl_' in level:  
        titlestr = titlestr +  fVstr(ifricV) + ','
    else:
        titlestr = titlestr +  fVstr(ifricV)  + ',' + fBstr(ifricB) + ','
        
    if iscouple == 1:
        
        titlestr = titlestr +\
            '\n iscale={0},stop_tol={1},ksatB={2},ksatV={3},dz={4},zmax={5}'.\
            format(iscale, stop_tol, phi_bare['ksat'][-1]*3600.,
                                phi_veg['ksat'][-1]*3600., dz, zmax)
        

    if itype1 != 0:
        titlestr = titlestr + 'itype1={0},'.format(itype1)
    if itype3 != 1:
        titlestr = titlestr + 'itype3={0},'.format( itype3)        
    #datastr = ''.join(titlestr.split(': ')[1].split('\n '))
    datastr = ''.join(','.join(titlestr.split(': ')).split('\n '))
    globals().update(locals())


### Functions to write the frinction part of the titlestring
def fVstr(ifricV ):
    if ifricV == 1:
        return 'fV=Manning{0}'.format(xnv)
    
    elif ifricV == 2:
        return  'fV=Cd-array{0}'.format(pveg)        

    elif ifricV == 3:
        return  'fV=Poggi'
        
    elif ifricV == 10:
        return  'fV=Kirstetter'
    
    elif ifricV == 11:
        return  'fV=DW17'    

    elif ifricV == 12:
        return  'fV=Poisseuille' 
    
    else:   
        return  'fV=test'
    

def fBstr(ifricB):
    
    if ifricB == 1:
        return 'fB=Manning{0}'.format(xni)

    elif ifricB == 2:
        return  'fB=Cd-array{0}'.format(pveg)        

    elif ifricB == 3:
        return  'fB=Poggi'
              
    elif ifricB == 10:
        return  'fB=Kirstetter'
    
    elif ifricB == 11:
        return  'fB=DW_17'    

    elif ifricB == 12:
        return  'fB=Poisseuille' 
    else:
        return  'fV=test'
    
    
def ftype(fnum):
    if fnum == 1:
        return 'Manning'

    elif fnum == 2:
        return  'Cd-array'

    elif fnum == 3:
        return  'Poggi'
    
    elif fnum == 10:
        return  'Kirstetter' 
    
    elif fnum == 11:
        return  'DW17'    

    elif fnum == 12:
        return  'Poisseuille'
    else:
      return 'Cd-test'
      
def summarize_exit():
    early_exit = 0
    for line in open('{0}/output/summary.out'.format( folder, datastr)):
        if line[:10] == '    quit a':
            a = (line.split(' ')[-1])
            quitting_time =  float(a[:-2])
            early_exit = 1
            
    if  early_exit == 0:

        shutil.copyfile('{0}/output/summary.out'.format(folder),
                       '{0}/error_summary.txt'.format( subfolder(), append))

    else:
        print 'quit at {0}!!  failed run:   '.format(quitting_time/60.)
        print titlestr

        #shutil.copyfile('{0}/output/summary.out'.format(folder),
        #                '{0}/exit_summary.txt'.format( subfolder()))


def subfolder():
    subfolder = '{0}/single'.format(folder)
    return subfolder   
        
if __name__ == '__main__':
    main(sys.argv)