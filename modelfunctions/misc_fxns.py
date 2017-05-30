import numpy as np
import pickle

def mybin(myarray, nbins = 10):
    len0 = len(myarray)
    cut = len0 - np.mod(len0, nbins)
    if cut < len0:
        myarray = myarray[:cut]

    if len0 -cut>1:
        print 'removing {0} rows from end'.format(len0-cut) 
        
    try:
        return myarray.reshape( myarray.shape[0]/nbins, 
                           nbins, myarray.shape[1]).mean(1)
    except:
        return myarray.reshape(myarray.shape[0]/nbins, nbins).mean(1)

def get_scalefactor(zcc, hdum):
    scalefactor = int(np.round(np.ptp(zcc)/np.ptp(hdum)))/3
    if scalefactor < 3:
        scalefactor = 1
    elif scalefactor > 3 and scalefactor < 8:
        scalefactor = 5
    elif scalefactor > 8 and scalefactor < 15:
        scalefactor  = 10
    elif scalefactor < 100:
        scalefactor = round(scalefactor, -1)
    else:
        scalefactor = round(scalefactor, -2)        
    return scalefactor        
    
### Functions to write the frinction part of the titlestring
