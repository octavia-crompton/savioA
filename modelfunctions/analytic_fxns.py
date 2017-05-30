def analytical_hydro(t , slope , nrow , tr ,  xn, dx, rain):
    """  
    analytical solution based on the kinematic wave theory, 
    provided by Stephenson and Meadows [1986] """
    import scipy.optimize
    import numpy as np
    alphaf  = (slope**.5)/xn 
    L = (nrow)*dx
    m = 5/3.
    if slope <= 0:
        print 'zero slope'           
        return np.zeros_like(t), 0
    try:
        tc = (L/alphaf/(rain)**(m-1))**(1./m)
        #print 'tc = {0:.1f} min, tr = {1:.1f} min'.format(tc/60., tr/60.)
    except: 
        return np.zeros_like(t), 0
    q = np.zeros_like(t)
    q[t<=tc] = alphaf*(rain*t[t<=tc])**m
    q[t>tc] = alphaf*(rain*tc)**m

    trdum = t[t>tr]- tr
    xx = np.zeros_like(trdum)
    for i in range(len(trdum)):
        def F(x):
            return rain*L - rain*m*alphaf**(1./m)*x**(1-1/m)*(t[t>tr]- tr)[i] - x
        x= scipy.optimize.broyden1(F, np.zeros_like(t[t>tr][1]), f_tol=1e-14)
        xx[i] = x

    q[t>tr] = xx
    q = q/L*3.6e5  # converting to cm/hr
    return q, tc
