import numpy as np

def vanGenuchten(h,phi) :
    alpha   = phi['alpha']
    theta_S = phi['theta_S']
    theta_R = phi['theta_R']
    n       = phi['n']
    m       = phi['m']
    Ksat    = phi['ksat'] 
    # Compute the volumetric moisture content
    theta = (theta_S - theta_R)/(1 + (alpha*abs(h))**n)**m + theta_R
    # Compute the effective saturation
    Se = ((theta - theta_R)/(theta_S - theta_R))
    # Compute the hydraulic conductivity
    K = Ksat*Se**(1./2)*(1 - (1 - Se**(1./m))**m)**2
    # Compute the specific moisture storage
    C =  -alpha*n*np.sign(h)*(1./n - 1)*(alpha*abs(h))**(n - 1)*(theta_R - 
         theta_S)*((alpha*abs(h))**n + 1)**(1/n - 2)
    try:
        for i in range(len(h)):
            if h[i] > 0:
                K[i] = Ksat[i]
                C[i] = 0.
                theta[i] = theta_S[i]
    except TypeError:
        if h > 0:
            K = Ksat[i]
            C = 0.
            theta = theta_S[i]
    return [C,K,theta]

def Richards_matrices(nz, infL = 2):
    # Define matrices that we'll need in solution
    DeltaPlus  = np.diag(-np.ones(nz)) + np.diag(np.ones(nz-1), 1)
    DeltaPlus[0,:] = 0
    DeltaPlus[nz-1,:] = 0

    DeltaMinus = np.diag(np.ones(nz)) + np.diag(-np.ones(nz-1),-1);
    DeltaMinus[0,:] = 0
    DeltaMinus[nz-1,:] = 0

    MPlus = np.diag(np.ones(nz))+np.diag(np.ones(nz-1),1)
    MPlus[0,0] = 2
    MPlus[0,1:nz-1] = 0
    MPlus[nz-1,nz-1] = 2
    MPlus[nz-1,:nz-1] = 0

    MMinus = np.diag(np.ones(nz)) + np.diag(np.ones(nz-1),-1)
    MMinus[0,0] = 2
    MMinus[0,1:nz-1] = 0 
    MMinus[nz-1,nz-1] = 2
    MMinus[nz-1,:nz-1] = 0 

    if infL == 2:
        MPlus[-2, -2] = 2
        MPlus[-2, -1] = 0
        
    return DeltaPlus, DeltaMinus, MPlus, MMinus

def pot_infl(hnp1m, thetan, phi, htop = 0, infL = 2, stop_tol = 1e-5, inv = 0):
    
    stop_flag = 0
    niter = 0

    while(stop_flag == 0):
        [cnp1m,knp1m,thetanp1m] = vanGenuchten(hnp1m,phi)
        Cdiag = np.diag(cnp1m) 
        kbarplus = (1/2.)*MPlus.dot(knp1m)
        Kbarplus = np.diag(kbarplus)
        kbarminus = (1/2.)*MMinus.dot(knp1m)
        Kbarminus = np.diag(kbarminus)
        A = (1./dt)*Cdiag - 1./((dz)**2)*(Kbarplus.dot(DeltaPlus) - 
                                          Kbarminus.dot(DeltaMinus)) 
        #  Compute the residual of MPFD (RHS)
        R_MPFD = (1./(dz**2))*(Kbarplus.dot(DeltaPlus).dot(hnp1m) - \
                               Kbarminus.dot(DeltaMinus).dot(hnp1m)) + \
             (1./dz)*(kbarplus - kbarminus) - (1./dt)*(thetanp1m - thetan); 
        # Compute deltam for iteration level m+1

        if inv == 0:
            Ainv = np.linalg.pinv(A)
        
        elif inv == 1:
            Ainv = np.linalg.inv(A)
            
        deltam = Ainv.dot(R_MPFD)

        # Increment iteration counter and display number of iterations
        niter = niter + 1;
        if niter > 100.:
            stop_tol = stop_tol*10.
            niter = 0 
        if (max(abs(deltam[1:(nz-1)]))<stop_tol):

            stop_flag = 1
            hnp1mp1 = hnp1m + deltam 
            hnp1mp1[0] = hnp1mp1[1]  # free drainage BC
            hnp1mp1[-1] = htop
            [cnp1m,knp1m,thetanp1m] = vanGenuchten(hnp1mp1,phi);
            knp1mp1 = knp1m
            cnp1mp1 = cnp1m            
            hnp1m = hnp1mp1
            kt = (knp1mp1[-infL] + knp1mp1[-2])/2.
            PI = kt*((hnp1mp1[-1] - hnp1mp1[-2])/dz + 1.)
        else:
            hnp1mp1 = hnp1m + deltam
            hnp1m = hnp1mp1
            hnp1m[0] = hnp1m[1] 
            hnp1m[-1] = htop
    return PI

def timestep(hnp1m, thetan, phi, setflux = 1, flux = 0, htop = 0., infL = 2,
              stop_tol = 1e-5, inv = 0):
    """"""
    stop_flag = 0
    niter = 0
    while(stop_flag == 0):
        [cnp1m,knp1m,thetanp1m] = vanGenuchten(hnp1m,phi)
        Cdiag = np.diag(cnp1m) 
        kbarplus = (1/2.)*MPlus.dot(knp1m)
        Kbarplus = np.diag(kbarplus)
        kbarminus = (1/2.)*MMinus.dot(knp1m)
        Kbarminus = np.diag(kbarminus)
        A = (1./dt)*Cdiag - 1./((dz)**2)*(Kbarplus.dot(DeltaPlus) - 
                                          Kbarminus.dot(DeltaMinus)) 
        R_MPFD = (1./(dz**2))*(Kbarplus.dot(DeltaPlus).dot(hnp1m) - \
                               Kbarminus.dot(DeltaMinus).dot(hnp1m)) + \
                 (1./dz)*(kbarplus - kbarminus) - (1./dt)*(thetanp1m - thetan) 

        if inv == 0:
            Ainv = np.linalg.pinv(A)
            
        elif inv == 1:
            Ainv = np.linalg.inv(A)
            
            
        deltam = Ainv.dot(R_MPFD)
        niter = niter + 1
        if niter > 100.:
            stop_tol = stop_tol*10.
            niter = 0        
            # print 'niter > 100, stop_tol increased to ', stop_tol
        if (max(abs(deltam[1:-1]))<stop_tol):
            stop_flag = 1
            hnp1mp1 = hnp1m + deltam 
            hnp1mp1[0] = hnp1mp1[1]  # free drainage BC
            if setflux == 1:                
                kt = (knp1m[-infL] + knp1m[-2])/2.
                hnp1mp1[-1] =  hnp1mp1[-2] - dz - flux*dz/kt
            else:
                hnp1mp1[-1] = htop
    
            [cnp1m, knp1m,thetanp1m] = vanGenuchten(hnp1mp1,phi);
            knp1mp1 = knp1m
            cnp1mp1 = cnp1m  
            hnp1m = hnp1mp1
 
        else:      
            hnp1mp1 = hnp1m + deltam
            hnp1m = hnp1mp1
            hnp1m[0] = hnp1m[1] 
            if setflux == 1:            
                kt = (knp1m[-infL] + knp1m[-2])/2.
                hnp1m[-1] = hnp1m[-2] - dz - flux*dz/kt
            else:
                hnp1m[-1] = htop
                     
    errornp1 =  np.sum((A.dot(deltam)))*dt*dz
    return hnp1mp1,cnp1mp1,knp1mp1, thetanp1m, errornp1, niter, stop_tol
    
    
def Richards(hinit, phi, nz = 101, nt= 100, depth = 0, infL = 2,
             stop_tol = 1e-5, htop = 0, ifixh = 0):

    # Define a storage container to store the pressure heads and soil moistures
    [Cinit,Kinit,thetainit] = vanGenuchten(hinit,phi) 
    H = np.zeros([nz,nt+1])
    H[:,0] = hinit
    THETA = np.zeros([nz,nt+1])
    THETA[:,0] = thetainit
    K = np.zeros([nz,nt+1])
    K[:,0] = Kinit

    depths = np.zeros([nt])
    iterations = np.zeros([nt])
    ktop = np.zeros([nt])
    kbot = np.zeros([nt])
    error = np.zeros([nt])
    iterations = np.zeros([nt])
    # Define the container for an iteration counter
    tp = 0
    start_time = time.time()

    # Initialize the Picard iteration solver
    for i in range(0, nt): 

        hnp1m =  H[:,i]  
        thetan = THETA[:,i]   
        if ifixh == 1 :
            hnp1mp1,cnp1mp1,knp1mp1, thetanp1mp1, errornp1,niter,tol = \
                        timestep(hnp1m,thetan,phi,setflux=0,
                                 htop=htop,stop_tol=stop_tol, infL = infL)
    
        elif (depth <= 0) and (prate[i] > 0):
            PI = pot_infl(hnp1m, thetan, phi,htop = htop, infL = infL, 
                          stop_tol = stop_tol)
            
            if  PI < prate[i]:
                if tp == 0:
                    print 'ponding at time t=', t[i]
                    tp = i
                flux = - PI
                hnp1mp1,cnp1mp1,knp1mp1,thetanp1mp1,errornp1,niter,tol  = \
                        timestep(hnp1m, thetan,phi,setflux = 0, 
                        htop = htop, stop_tol=stop_tol, infL = infL)
                if ponding == 1:
                    depth += prate[i]*dt + flux*dt
            else:
                flux = - prate[i]
                hnp1mp1,cnp1mp1,knp1mp1,thetanp1mp1,errornp1,niter,tol  = \
                   timestep(hnp1m,thetan,phi,setflux = 1, flux = flux, 
                            stop_tol=stop_tol, infL = infL)

        elif (depth <= 0) and (prate[i] == 0):   
            hnp1mp1,cnp1mp1,knp1mp1,thetanp1mp1,errornp1,niter,tol = \
               timestep(hnp1m,thetan,phi,setflux = 1, flux = 0., 
                        stop_tol=stop_tol, infL = infL)                

        elif (depth > 0):  
            hnp1mp1,cnp1mp1,knp1mp1,thetanp1mp1,errornp1,niter,tol  = \
                            timestep(hnp1m,thetan, phi,setflux=0, htop = depth, 
                                     stop_tol=stop_tol, infL = infL)
            kt = (knp1mp1[-infL] + knp1mp1[-2])/2.       
            if ponding == 1:
                depth +=  prate[i]*dt - kt*((hnp1mp1[-1] - hnp1mp1[-2])/dz + 1.)*dt 

        if depth <0.:
            depth = 0.

        THETA[:,i+1] = thetanp1mp1
        H[:,i+1] = hnp1mp1
        K[:,i+1] = knp1mp1

        kbot[i] =  (K[1, i+1] + K[0, i+1])/2.  
        ktop[i] =  (K[-infL, i+1] + K[-2, i+1])/2.  
        error[i] = errornp1
        iterations[i] = niter
        depths[i] = depth

        if np.mod(i,1) == 10:  
            if ifixh == 1:
                print 't = {0:.2f}; '.format( t[i])

            elif prate[i] > 0:
                print 't = {0:.2f}; depth = {1:.2f}cm; PI = {2:.2f}cm/hr; rain = {3:.2f}cm/hr '.format(
                     t[i], depth, PI*3600,prate[i]*3600)
            else:
                print 't = {0:.2f}; depth = {1:.2f}cm;  '.format(
                     t[i], depth)

    t_elapse = (time.time() - start_time)
    
    newmass =  (THETA[:, 1:] - THETA[:, :-1]).sum(0)*dz   #  change in mass      

    H = H[:, 1:]
    THETA = THETA[:, 1:]
    K = K[:, 1:]
    fluxin = ktop*((H[-1]-H[-2])/dz + 1.)*dt   # top flux (cm)
    fluxout = - kbot*((H[1] - H[0])/dz + 1. )*dt # bottom flux (cm)

    return H, THETA, K, fluxin, fluxout, newmass, error, depths, t_elapse