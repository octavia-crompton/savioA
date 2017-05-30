import numpy as np
import os, sys
# functions to read fortran output
def strip_line():
    line = f.next() 
    line = f.next()
    a = line.strip().split(" ")
    a = [np.float(b.replace('d', 'e')) for b in a if b]
    return a

def myfloat(b):
    
    try: 
        b = float(b)
    except ValueError:        
        b = [b for b in b.split('-') if b]
        b = float(b[0])*10**(-float(b[1]))
    return b

    
def read_param(folder = 'test') :   
 
    fname = '{0}/input/params.dat'.format(folder)    
    f = open(fname, 'r')
    a = strip_line()
    dt = a[1]
    tmax = a[2]
    a = strip_line()
    epsh = a[0]
    beta = a[1]
    a = strip_line()
    xk = a[0]
    nprt = a[1]
    a = strip_line()
    iscouple = a[0]
    ifixh = a[1]
    a = strip_line()
    iscale = a[0]
    a = strip_line()
    stop_tol = a[0]
    f.close()
    return dt, tmax, nprt, epsh, iscale, ifixh, stop_tol

def read_time(folder = 'test') :  
    """
    read time.out
    
    """
    tp = []
    itp = []  #  print step
    it = []   #  time step
    f =  open("{0}/output/time.out".format(folder), 'r')
    f.next()
    for line in f:
        a = (line.strip().split(" "))
        a = [b for b in a if b]
        tp.append(float(a[0]))
        itp.append(int(a[1]))
        it.append(int(a[2]))

    tp = np.array(tp)
    
    return tp

def read_cfl(folder = 'test'):
    f =  open("{0}/output/cfl.out".format(folder), 'r')
    cfl = []
    ta = []
    ftime = []
    for line in f:
        a = line.strip().split(" ")
        a = [np.float(b) for b in a if b]
        ta.append(a[0])
        cfl.append(a[1])
        ftime.append(float(a[2]))
    cfl =  np.array(cfl)
    ta = np.array(ta)
    ftime = np.array(ftime)
    
    return ta, cfl, ftime


def new_hydro(folder , fname = 'hydro', Lx=1, Ly=1):
    """ hydro units are f*ds = m^3/s, converted to cm/s"""
    tdum = []
    hydro = []

    for line in open("{0}/output/{1}.out".format(folder, fname), 'r'):
        a = (line.strip().split(" "))
        a = [float(b) for b in a if b]
        try:
            tdum.append(float(a[0]))
            hydro.append(float(a[1]))

        except IndexError:
            raise KeyboardInterrupt

    tdum = np.array(tdum)
    hydro = -np.array(hydro)/Lx/Ly*100  # now  cm/s   
    
    return tdum, hydro # now cm/s

def get_h(folder, fdir = 'output/h.out', ncol = 1, nrow = 1):
    h = []
    hdum =  np.zeros([ncol+2, nrow+2])
    v = []
    vdum =  np.zeros([ncol+2, nrow+2])
    u = []
    udum =  np.zeros([ncol+2, nrow+2])

    for line in open("{0}/{1}".format(folder, fdir), 'r'):
        a = (line.strip().split(" "))
        a = [float(b) for b in a if b]
        try:
            j = int(a[0])
            k = int(a[1])
            hdum[j, k] = a[2]
            udum[j, k] = a[3]
            vdum[j, k] = a[4]
        except IndexError:
            dumt = int(a[0])
            h.append(hdum.copy())    
            u.append(udum.copy())
            v.append(vdum.copy())

    h = np.array(h)
    h = h[:, 1:-1, 1:-1]
    u = np.array(u)
    u = u[:, 1:-1, 1:-1]
    v = np.array(v)
    v = v[:, 1:-1, 1:-1]
    return h,u,v

def get_zc(folder = 'test'):
    zcdum = np.zeros([ncol+2, nrow+2])
    for line in open("{0}/output/{1}.out".format(folder, 'zcc'), 'r'):
        a = (line.strip().split(" "))
        a = [str(b) for b in a if b]        
        j = int(a[0])
        k = int(a[1])    
        zcdum[j,k] = float(a[2])
    zc = zcdum[1:-1, 1:-1]
    
    return zc

        
def get_dvol(folder, Lx, Ly):
    # dvol, flux, infl, rain from fortran mass tracking
    # converted here from 
    ta = []
    dvol = []
    infl = [] 
    flux = []
    rain = []
    dts = []
    f = open('{0}/output/dvol.out'.format(folder), 'r'); 
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        ta.append(a[0])
        dvol.append(a[1])
        flux.append(a[2])
        infl.append(a[3])
        rain.append(a[4])
        dts.append(a[5])
    ta = np.array(ta)
    dvol = np.array(dvol)
    flux = np.array(flux)
    infl = np.array(infl)
    rain = np.array(rain)    
    dts = np.array(dts)    
    total = dvol  - flux - infl - rain 
    
    ti = 0; tf = len(ta)-1;  nprt = 1
    
    dvol = dvol/Lx/Ly*100
    flux = flux/Lx/Ly*100
    infl = infl/Lx/Ly*100
    rain = rain/Lx/Ly*100
    return ta, dvol, flux, infl, rain, total

def lateral_fluxes(folder):
    """ in units m^3"""
    flux1 = []
    flux2 = []
    flux3 = [] 
    flux4 = []
    f = open('{0}/output/fluxes1234.out'.format(folder), 'r'); 

    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        flux1.append(a[0])
        flux2.append(a[1])
        flux3.append(a[2])
        flux4.append(a[3])

    flux1 = -np.array(flux1)
    flux2 = -np.array(flux2)
    flux3 = -np.array(flux3)
    flux4 = -np.array(flux4)
    return flux1, flux2, flux3, flux4
       
def read_soil_fluxes11(folder, dt):
    """ fluxes from fortran subroutine.
        units = :
    """
    fluxin_11 = []
    fluxout_11 = []
    newmass_11 = []
    error_11 = []
    for line in open("{0}/output/soilfluxes.out".format(folder), 'r'):
        a = (line.strip().split(" "))
        a = [str(b) for b in a if b]        
        fluxin_11.append(float(a[1]))
        fluxout_11.append(float(a[2]))
        newmass_11.append(float(a[3]))
        error_11.append(float(a[4]))

    fluxin_11 = np.array(fluxin_11)*3.6e3/dt
    fluxout_11 = np.array(fluxout_11)*3.6e3/dt
    newmass_11 = np.array(newmass_11)*3.6e3/dt
    error_11 = np.array(error_11)*3.6e3/dt
    return fluxin_11, fluxout_11, newmass_11, error_11
   
def read_fluxmap(folder, ncol, nrow, dt_sw):
    """
    fluxmap : map of surface flux  from fortran richards solver.
        fortran units =  cm
        recorded only when Richards solver is called
        
    fluxmap : map of estimated error from fortran richards solver.        
    """
    fluxin_map = []
    fluxout_map = []
    newmass_map = []    
    error_map = []
    fluxin_dum  = np.ones([ncol, nrow])
    fluxout_dum  = np.ones([ncol, nrow])
    newmass_dum  = np.ones([ncol, nrow])
    error_dum  = np.ones([ncol, nrow])
    t0 = dt_sw
    for line in open("{0}/output/soilfluxmap.out".format(folder), 'r'):

        a = (line.strip().split(" "))
        a = [str(b) for b in a if b]   
        j = int(a[0])-1
        k = int(a[1])-1
        t1 =float(a[2])
        if t1<=t0:
            fluxin_dum[j,k] = float(a[3])
            fluxout_dum[j,k] = float(a[4])
            newmass_dum[j,k] = float(a[5])
            error_dum[j,k] = float(a[6])            
        elif t1>t0:
            t0 = t1       
            fluxin_map.append(fluxin_dum.copy())
            fluxout_map.append(fluxout_dum.copy())
            newmass_map.append(newmass_dum.copy())
            error_map.append(error_dum.copy())
            
            fluxin_dum[j,k] = float(a[3])
            fluxout_dum[j,k] = float(a[4])
            newmass_dum[j,k] = float(a[5])
            error_dum[j,k] = float(a[6])     

    fluxin_map.append(fluxin_dum.copy())        
    fluxout_map.append(fluxout_dum.copy())        
    newmass_map.append(newmass_dum.copy())        
    error_map.append(error_dum.copy())        
    
    fluxin_map = np.array(fluxin_map)
    fluxout_map = np.array(fluxout_map)    
    newmass_map = np.array(newmass_map)    
    error_map = np.array(error_map)
    

    return fluxin_map,fluxout_map, newmass_map, error_map
     
def read_zinflmap(folder, nrow, ncol, dt_sw, dx):
    """
     zinflmat : map of infiltration from fortran mass tracking. 
         fortran units = m^3
         converted here to cm/hr
         recorded every time step  (code redundancy to fluxmap) 
    """
    zinflmap = []
    zinfldum  = np.ones([ncol, nrow])
    t0 = dt_sw

    for line in open("{0}/output/zinflmap.out".format(folder), 'r'):
        a = (line.strip().split(" "))
        a = [str(b) for b in a if b]   
        j = int(a[0])-1
        k = int(a[1])-1
        t1 =float(a[2])
        zinfldum[j,k] = float(a[3])
        if t1>t0:
            t0 = t1
            zinflmap.append(zinfldum.copy())
    zinflmap.append(zinfldum.copy())        
    zinflmap = -np.array(zinflmap)/dx**2*100   # convert  units m^3 to cm
    return zinflmap
  