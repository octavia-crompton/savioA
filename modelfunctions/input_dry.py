import numpy as np
import scipy as sp
import os, sys

def write_prate(folder = 'test', tmax = 100, dt_sw = 1., tr = 20*60, rain = 0.036/3600.):
    t = np.arange(0., tmax + dt_sw, dt_sw)
    nt_sw = t.shape[0]-1
    prate = np.ones_like(t)*rain
    prate[np.where(t>=tr)[0][0]:] = 0

    f = open('{0}/input/prate.dat'.format(folder), 'w')

    f.write('{0:<13} \n'.format(nt_sw))
    for n in range( nt_sw +1):
        f.write('{0} {1} \n'.format(t[n], prate[n])) 

    f.close()
    return nt_sw

def build_coords(folder = 'test', ncol = 2, nrow = 2, dx = 1, 
                 slope = 0., seed = 0, level = '', topo = 'plane', 
                veg = 'none'):
    np.random.seed(seed)
    npt = (ncol+1)*(nrow+1)  # number of points
    ne = nrow*ncol  # number of edges
    nbcell = 2*ncol + 2*nrow - 4  # number of boundary cells

    xdum = np.arange(0, (ncol+1)*dx - 1e-10, dx )
    ydum = np.arange(0, (nrow+1)*dx - 1e-10, dx )
    ydum, xdum = np.meshgrid(ydum, xdum)

    zymax = slope*(np.max(ydum) - np.min(ydum))
    zxmax = slope*(np.max(xdum) - np.min(xdum))*2


    if topo == 'gnoise': 
        rowNoise = sp.randn(nrow+1)/20.
        rowGauss =  gaussian_filter(rowNoise, 5, order=1, \
                                output=None, cval=0.0, truncate=20.0)
        zrow = rowGauss + np.linspace(0,zymax, nrow+1)

        colNoise = sp.randn(ncol+1)/20.
        colGauss =  gaussian_filter(colNoise, 2, order=1, \
                                output=None, cval=0.0, truncate=20.0)
        zcol = colGauss + np.abs(np.linspace(-zxmax, zxmax, ncol+1))

        zrow2 = np.tile(zrow, [ncol+1]).reshape([ncol+1, nrow+1])
        zcol2 = np.tile(zcol, [nrow+1]).reshape([nrow+1, ncol+1])

        zdum = zrow2 + 0*zcol2.T
        
    elif topo == 'plane':
        zrow =  np.linspace(0,zymax, nrow+1)
        zrow2 = np.tile(zrow, [ncol+1]).reshape([ncol+1, nrow+1])
        zdum = zrow2 
        
    elif topo == 'sine':
        zrow =  np.linspace(0,1, nrow+1)
        zrow = zymax*np.sin(zrow*2*np.pi/4)
        zrow2 = np.tile(zrow, [ncol+1]).reshape([ncol+1, nrow+1])
        zdum = zrow2 
        
    elif topo == 'log':
        
        zrow =  np.linspace(0,zymax, nrow+1)
        zcol = np.abs(np.linspace(-zxmax, zxmax, ncol+1))
        y0 = np.log(np.abs(np.linspace(-zxmax, zxmax, ncol+1))*logscale+1)
        y0 = y0*zxmax/np.max(y0)/2
        y1 = np.ones_like(y0)*np.max(zrow)
        zdum = np.ones_like(xdum)
        for i in range(len(y0)):
            zdum[i, :] = (np.linspace(y0[i], y1[i], nrow+1))
    elif topo == 'flat':
        zdum = np.zeros_like(xdum)
        
    xndum = np.ones_like(zdum)
    isveg_dum = np.zeros_like(zdum)    
    
    # Infer vegetation field from level
    if level[0] == 'A':
        if 'bare' not in level:
            isveg_dum = np.ones_like(zdum)    
        print 'level = ', level, '!'
    if level[0] == 'R' and  veg == 'veg' :
        isveg_dum = np.ones_like(zdum)    
    if veg == 'stripe1HL':
        print 'level = ', level, '!'  
        isveg_dum[:, :nrow/2]  = 1   
    elif veg == 'stripe2HL':
        print 'level = ', level, '!'  
        isveg_dum[:, :nrow/4]  = 1           
        isveg_dum[:, nrow/2:nrow*3/4]  = 1                   
    elif veg =='randv':
        isveg_dum = sp.rand(ncol+1, nrow+1) > .9
    
    elif veg == 'upslope':
        isveg_dum = (zdum > np.mean(zdum))
    xcc, ycc, zcc, xncc, isvegcc = write_coords(folder = folder, 
                                    ncol = ncol, nrow= nrow, dx = dx, 
                                     xdum = xdum, ydum = ydum, xndum = xndum,
                                     zdum = zdum, isveg_dum = isveg_dum)
    isvegcc[isvegcc > 0] = 1
    return   xcc, ycc, zcc, xncc, isvegcc
    
def write_coords(folder, ncol, nrow, dx, xdum, ydum, zdum, xndum, isveg_dum ):        
    
    npt = (ncol+1)*(nrow+1)  # number of points
    ne = nrow*ncol  # number of edges
    nbcell = 2*ncol + 2*nrow - 4  # number of boundary cells

    
    x = np.zeros(npt + 1)
    y = np.zeros(npt + 1)
    z = np.zeros(npt + 1)
    xn = np.zeros(npt + 1)
    isveg = np.zeros(npt + 1)

    x[1:] = xdum.ravel()
    y[1:] = ydum.ravel()
    z[1:] = zdum.ravel()
    xn[1:] = xndum.ravel()
    isveg[1:] = isveg_dum.ravel()

    # (ncol+1) by (nrow+1)  -  node numbers
    nodes = np.arange(1, npt+1, dtype = int).reshape([ncol+1, nrow+1])

    nop = np.zeros([ncol+1, nrow+1, 4], dtype = int)
    for j in range(ncol):
        for k in range(nrow):
            nop[j+1, k+1] =  nodes[j,k], nodes[j+1, k], nodes[j+1,k+1], nodes[j,k+1]
           

    fname = '{0}/input/coords.dat'.format(folder)
    f = open(fname, 'w')
    f.write('{0:<13}   {1:<13}\n'.format(npt, ne))

    # write x, y, z
    for n in range(1, npt+1):
        f.write('{0:<13.6f} {1:<13.6f} {2:<13.6f} {3:<13.6e} {4:<13.6e}  \n'.format(
                    x[n],y[n],z[n],xn[n],isveg[n])) 

    # write node numbers  
    for j in range(1, ncol+1):
        for k in range(1, nrow+1):
            n1 = nop[j, k, 0] 
            n2 = nop[j, k, 1]       
            n3 = nop[j, k, 2]        
            n4 = nop[j, k, 3] 
            f.write('{0:<10} {1:<10}  {2:<10} {3:<10}\n'.format(n1, n2, n3, n4)) 
    f.close()  

    # get cell center values:
    xcc  = np.zeros([ncol+2, nrow+2])    
    ycc  = np.zeros([ncol+2, nrow+2])
    zcc  = np.zeros([ncol+2, nrow+2])
    xncc  = np.zeros([ncol+2, nrow+2])
    isvegcc  = np.zeros([ncol+2, nrow+2])

    for j in range(1, ncol+1):
        for k in range(1, nrow+1):
            n1 = nop[j, k, 0] 
            n2 = nop[j, k, 1]       
            n3 = nop[j, k, 2]        
            n4 = nop[j, k, 3]  
            xcc[j,k] = 0.25*(x[n1] + x[n2] + x[n3] + x[n4])  
            ycc[j,k] = 0.25*(y[n1] + y[n2] + y[n3] + y[n4])
            zcc[j,k] = 0.25*(z[n1] + z[n2] + z[n3] + z[n4])   
            xncc[j,k] = 0.25*(xn[n1] + xn[n2] + xn[n3] + xn[n4])   
            isvegcc[j,k] = 0.25*(isveg[n1] + isveg[n2] + isveg[n3] + isveg[n4])   
                
    return xcc[1:-1, 1:-1], ycc[1:-1, 1:-1], zcc[1:-1, 1:-1], xncc[1:-1, 1:-1], \
            isvegcc[1:-1, 1:-1]
            
def write_param(folder = 'test',  dt_sw = 1, 
                 tmax = 100,  nprt = 1, epsh = 0.0025, iscouple = 0,
                 iscale= 1, ifixh = 0, stop_tol = 1e-5, depth = 0.,
                 beta = 1.0, ifricV = 1, ifricB = 1, 
                 xnv = 0.1, xni = 0.03, pveg = .2, Dveg = 8e-3, tr = 20):
    
    fname = '{0}/input/params.dat'.format(folder)    
    f = open(fname, 'w')
    f.write('gravity     dt           \n')
    f.write('9.806d0     {0}          \n'.format(dt_sw))
    f.write('tmax       tr    \n')
    f.write('{0}       {1}     \n'.format( tmax,tr))    
    f.write(' epsh      beta     \n')  
    f.write('{0}        {1}     \n'.format(epsh, beta))
    f.write('xk          nprt      \n')
    f.write('3.9217d-4   {0}       \n'.format(int(nprt)))   
    f.write(' iscouple    ifixh \n')
    f.write(' {0}        {1} \n'.format( iscouple, ifixh))   
    f.write(' iscale   \n')
    f.write(' {0}       \n'.format( iscale))  
    f.write(' stop_tol   \n')
    f.write(' {0}       \n'.format( stop_tol))      
    f.write('h0      u0    v0   \n ')
    f.write('{0}    0.0    0.0  \n '.format(depth))   
    f.write('ifricV     ifricB   \n ')
    f.write('{0}        {1}      \n '.format(ifricV, ifricB))       
    f.write('xnv     xni   \n')
    f.write(' {0}    {1}\n'.format( xnv, xni))  
    f.write('pveg     D   \n')
    f.write('{0}    {1}\n'.format( pveg, Dveg))      
    f.close()
    
def write_dryin(folder = 'test', ncol = 2, nrow = 2, itype3 = 1, itype1 = 0,
               itype2 = 1,itype4 = 1,):
    inum = np.zeros([ncol+1, nrow+1], dtype = int)
    inum[1:, 1] = 1
    inum[1:, -1]= 1
    inum[1, 1:] = 1
    inum[-1, 1:] = 1
    inum[1, 1] = 2
    inum[1, -1] = 2
    inum[-1, -1] = 2
    inum[-1, 1] = 2
    
    ipos = np.zeros( [ncol+1, nrow+1, 2], dtype = int)
    # bottom boundary
    ipos[2:-1, 1,0] = 1
    ipos[1, 1,1] = 1
    ipos[-1, 1,1] = 1

    # right boundary
    ipos[-1, 1:-1, 0] = 2
    ipos[-1, -1,1] = 2

    # left boundary
    ipos[1, 1:, 0] = 4

    # top boundary
    ipos[2:, -1,0] = 3
    ipos[1, -1,1] = 3
    
    itype = np.zeros([ncol+1, nrow+1, 2], dtype = int)
    # bottom boundary
    itype[2:-1, 1,0] = itype1
    itype[1, 1,1] = itype1
    itype[-1, 1,1] = itype1

    # right boundary
    itype[-1, 1:-1, 0] = itype2
    itype[-1, -1,1] = itype2

    # left boundary
    itype[1, 1:,0] = itype4

    # top boundary
    itype[2:, -1,0] = itype3
    itype[1, -1,1] = itype3
    
    npt = (ncol+1)*(nrow+1)  # number of points
    ne = nrow*ncol  # number of edges
    nbcell = 2*ncol + 2*nrow - 4  # number of boundary cells
    
    fname = '{0}/input/boundary.dat'.format(folder)    
    f = open(fname, 'w')
    f.write('number of boundary cell \n') 
    f.write('  {0} \n'.format(nbcell))
    f.write(' j    k          inum    itype             ipos \n')
    # f.write(' j \t k \tinum    itype \t\t ipos')
    j = 1
    for k in range(1, nrow+1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    for j in range(2, ncol+1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    for k in range(nrow-1,0,-1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))
            
    for j in range(ncol-1,1,-1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1], 
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    kbeg = np.ones(ncol+1, dtype = int)
    kend = np.ones(ncol+1, dtype = int)*nrow
   
    f.write('ncol\n')
    f.write("{0}\n".format(ncol))
    f.write('nrow\n')
    f.write("{0}\n".format(nrow))    
    f.write('j     kbeg          kend \n')
    for j in range(1, ncol+1):
        f.write( '{0:>5}  {1:>5} {2:>13}   \n'.format(
                    j, kbeg[j],kend[k] ))

        
    f.close()
    boundary_fix(folder = folder, nrow = nrow, ncol = ncol, )
    return inum, ipos, itype

def boundary_fix(folder = 'test', nrow = 0, ncol = 0, Ly = 0, Lx = 0):  # add dx to get Ly and Lx
    fname = '{0}/input/boundary.dat'.format(folder)    
    f = open(fname, 'a')
    fixj = np.array([1,2])
    fixk = np.array([nrow, nrow])
    fixh = np.array([0,0]) 
    fixu = np.array([0,0])
    fixv = np.array([-1., -1.])*1e-5  # q = .36cm/hr
    
    #fixv = np.array([-1., -1.])/3.6e3/Ly  # q = 1cm/hr    
    #q = 1.0e-3/Lx/Ly
    #fixv = np.array([-q, -q])  # q = 1cm/hr    
    ndir = len(fixj)
    f.write('number of fixed bc cells, ndir    \n ')
    f.write('{0}  \n '.format(ndir))   
    f.write('j     k    fix h    fix u    fix v	\n')
    for i in range(ndir): 
        f.write('{0}     {1}     {2}     {3}     {4}     \n'.format(
                fixj[i],fixk[i],fixh[i],fixu[i],fixv[i]))
    f.close()
    
def write_inc(folder= 'test', nt = 100, nz = 101, ncol = 2, nrow = 2):
    #  Avoid allocating to much space 
    with open('dry1.inc', 'r') as input_file, open('{0}/dry.inc'.format(folder), 'w') as output_file:
        for line in input_file:
            if line[6:15] == 'parameter':
                # print 'old line:', line
                a = (line.strip().split(" ")[2].split(","))
                nn0 = int(a[0].split('=')[-1])
                nt0 = int(a[1].split('=')[-1])
                nx0 = int(a[2].split('=')[-1])
                ny0 = int(a[3].split('=')[-1])
                nz0 = int(a[4].split('=')[-1])

                newline = '      parameter ( nn={0},ntp={1},nx={2},ny={3},nz={4} )'.format(
                     (ncol+1)*(nrow+1)+1, nt+1, ncol+4, nrow+4, nz)

                output_file.write(newline)
            else:
                output_file.write(line)

    