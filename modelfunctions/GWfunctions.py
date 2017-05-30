import numpy as np
# GW functions
def findh(teq, param):
    # param=[tp,tca,Ksat,I,Ao,Kr_flat]
    tp=param[0]
    tca=param[1]
    Ksat=param[2]
    I=param[3]
    Ao = param[4]
    Kr = param[5]
    res=((I - Ksat)*(teq-tp)-Ao*(np.sqrt(teq-tp+tca)-np.sqrt(tca)))
    return res

def findt(teq,param):
    tp=param[0]
    tca=param[1]
    Ksat=param[2]
    I=param[3]
    Ao = param[4]
    Kr = param[5]
    x= param[6]
    xo = param[7]
    dt =  (teq-tp)/100.
    t = np.arange(tp, teq+dt/10., dt)
    h = ((I - Ksat)*(t-tp)-Ao*(np.sqrt(t-tp+tca)-np.sqrt(tca)))
    h[h<0]=0
    hpower = h**(2./3)
    hint = np.sum((hpower[:-1]+hpower[1:])/2.*(t[1:]-t[:-1]))+\
            hpower[0]/2.*(t[1]-t[0]) + hpower[-1]/2.*(t[-1]-t[-2])
    res=5/3.*Kr*hint-x+xo
    return res


def findtinit(teq,param):
    tpd=param[0]
    tca=param[1]
    Ksat=param[2]
    I=param[3]
    Ao = param[4]
    Kr = param[5]
    L = param[6]
    tp = param[7]  
    dt =  (tpd-teq)/100.
    t = np.arange(teq, tpd+dt/10., dt)
    h = ((I - Ksat)*(t-teq)-Ao*(np.sqrt(t-tp+tca)-np.sqrt(teq-tp+tca)))    
    h[h<0]=0
    hpower = h**(2./3)
    hint = np.sum((hpower[:-1]+hpower[1:])/2.*(t[1:]-t[:-1]))+\
            hpower[0]/2.*(t[1]-t[0]) + hpower[-1]/2.*(t[-1]-t[-2])
               
    res=5/3.*Kr*hint-L
    return res

def findto_out(teq,param):
    to=param[0]
    tca=param[1]
    Ksat=param[2]
    I=param[3]
    Ao = param[4]
    Kr = param[5]
    L= param[6]
    tp=param[7]
    dt = (teq-to)/100.
    t = np.arange(to, teq+dt/10., dt )

    h = (I-Ksat)*(t-to)-Ao*np.sqrt(t-tp+tca)+Ao*np.sqrt(to-tp+tca)
    h[h<0]=0.
    hpower = h**(2./3)
    hint = np.sum((hpower[:-1]+hpower[1:])/2.*(t[1:]-t[:-1]))+\
            hpower[0]/2.*(t[1]-t[0]) + hpower[-1]/2.*(t[-1]-t[-2])

    res=5./3*Kr*hint-L
    return res

def  findt_star(teq,param):
    h_star=param[0]
    Ksat=param[1]
    tpd=param[2]
    Ao=param[3]
    tp=param[4]
    tca=param[5]
    res= h_star-Ksat*(teq-tpd)-Ao*(teq-tp+tca)**(1./2)+Ao*(tpd-tp+tca)**(1./2)
    return res

def findt_cross(teq,param):
    tpd=param[0]
    tca=param[1]
    Ksat=param[2]
    x_star=param[3]
    Ao = param[4]
    Kr = param[5]
    L= param[6]
    tp=param[7]
    h_star=param[8]
    dt = (teq-tpd)/100
    t = np.arange(tpd,teq+dt/10., dt)
    h = h_star-Ksat*(t-tpd)-Ao*np.sqrt(t-tp+tca)+Ao*np.sqrt(tpd-tp+tca)
    h[h<0]=0
    hpower = h**(2./3)
    hint = np.sum((hpower[:-1] +hpower[1:] )/2.*(t[1:]-t[:-1]))+\
                  hpower[0]/2.*(t[1] - t[0])+hpower[-1]/2.*(t[-1]-t[-2])
    res = 5./3*Kr*hint-L+x_star
    return res

def Comparison_function2(I,L,Ksat,Ao,tpd,alpha):
    #   function called as Comparison_function2(I,L,Ksat,Ao,tf,alpha) 
    #   tpd : storm duration (s), tf
    #tpd = tr
    g = 9.81  # m/s^2 Gravity
    dx =  L/500.
    x = np.arange(0, L+dx, dx)
    m = len(x)
    if I > Ksat:
        tp_flat = Ao**2/(I-Ksat)**2
    else:
        tp_flat = tpd*10.

    tca_flat = 1/Ksat**2*(- Ao/2+(Ao**2/4+ I*tp_flat*Ksat)**(1/2.))**2

    if tpd>tp_flat:
        tpd_flat = tpd-tp_flat  #   s      Storm Duration from time of ponding onwards
        # Kr_flat = So**(1/2.)*1./n;  # Manning computation of the resistance
        Kr_flat = alpha           # alpha is the hydraulic resistance 

        [tout_flat,qout_flat] = Steady_flat_w_infil2(I,g,L,tpd,tp_flat,tca_flat,Ksat,Ao,Kr_flat);

    elif tpd<tp_flat:       #  Time of ponding occurs after storm ends!
        qout_flat=np.zeros(100)
        tout_flat = np.arange(1,101)
    ## COMPUTE PARTITIONING
    # Total Rainfall 
    PPT = I*tpd*L     # m^3/m
    # Area under hydrograph
    Q_flat= sum((qout_flat[:-1]+qout_flat[1:])/2.*(tout_flat[1:]-tout_flat[:-1]))+\
            (qout_flat[0]/2.*(tout_flat[1]-tout_flat[0])+ \
            qout_flat[-1]/2.*(tout_flat[-1]-tout_flat[-2]))

    # Percent runoff
    Perc_runoff_flat = Q_flat/PPT
    
    output = np.vstack((tout_flat, qout_flat)).T
    output = output[output[:, 0].argsort()]
    tout = output[:,0]
    qout = output[:,1]

    return [tout,qout,Perc_runoff_flat] 

def Steady_flat_w_infil2(I,g,L,tpd,tp,tca,Ksat,Ao,Kr_flat):
    #  Steady_flat_w_infil2(I,g,L,tpd,tp_flat,tca_flat,Ksat,Ao,Kr_flat);
    # tp = tp_flat
    # tca = tca_flat
    
    x = L
    param=[tp,tca,Ksat,I,Ao,Kr_flat]

    from scipy.optimize import fsolve
    tmin = fsolve(findh, tp, param)[0]
    xo = 0
    param=[tp,tca,Ksat,I,Ao,Kr_flat,x,xo]

    Teq=fsolve(findt,tmin+.1,param)  

    if tpd>Teq:
        T = np.arange(np.ceil(tp), np.round(Teq),(Teq-tp)/100.)
        h = np.ones_like(T)
        q = np.ones_like(T)
        for j in range(len(T)):
            t=T[j]
            h[j]=((I - Ksat)*(t-tp)-Ao*(np.sqrt(t-tp+tca)-np.sqrt(tca)))
            if h[j] < 0:
                h[j] = 0    
            q[j]=Kr_flat*h[j]**(5./3)

        dim = max(np.round(tp)-1, 0)
        q = np.hstack((np.zeros(np.int(dim)), q))
        t1= np.hstack((np.arange(1,np.int(np.round(tp))), T))
        # Domain 2:  Characteristics originating from x=0 for t>Teq
        # Find the last time at which a characteristic originating from x=0 at
        # t>Teq can cross the x=L line before the rainfall ends

        param=[tpd,tca,Ksat,I,Ao,Kr_flat,L,tp]
        tinit = fsolve(findtinit,tp,param)

        To = np.arange( np.round(tp)+1, np.round(tinit),(tinit-tp)/100. )
        to_out = np.ones_like(To)
        h_out = np.ones_like(To)
        q_out = np.ones_like(To)
        for j in range(len(To)):
            to=To[j]
            # Find the value of t at which the characteristic originating at to
            # reaches x = L
            param=[to,tca,Ksat,I,Ao,Kr_flat,L,tp]
            to_out[j] = fsolve(findto_out,tpd,param)
            h_out[j] = (I-Ksat)*(to_out[j]-to)-Ao*(np.sqrt(to_out[j]-tp+tca)-np.sqrt(to-tp+tca))
            q_out[j] = Kr_flat*h_out[j]**(5./3)

        # Domain 3: Characteristics originating from x = 0 which persist into the
        # falling limb of the hydrograph
        dt = (tpd-tinit)/100
        To2 = np.arange(tinit+(tpd-tinit)/100, tpd-(tpd-tinit)/100+ dt/10., dt)
        x_star = np.ones_like(To2)
        h_star = np.ones_like(To2)
        t_star = np.ones_like(To2)
        x_zw = np.ones_like(To2)
        t_cross = np.zeros_like(To2)
        h_cross = np.zeros_like(To2)    
        q_cross = np.zeros_like(To2)        
        for j in range(len(To2)):
            to = To2[j]
            # Find the value of x (x_star) at which the characteristic reaches t = tpd
            dt = (tpd-to)/100.
            t = np.arange(to, tpd + dt/10., dt)

            h = (I-Ksat)*(t-to)-Ao*np.sqrt(t-tp+tca)+Ao*np.sqrt(to-tp+tca)
            hpower = h**(2./3)
            hint = np.sum((hpower[:-1]+hpower[1:])/2.*(t[1:]-t[:-1]))+ \
                    hpower[0]/2.*(t[1]-t[0]) + hpower[-1]/2.*(t[-1]-t[-2])
            x_star[j] = 5./3*Kr_flat*hint
            # Find the value of h (h_star) at this point, (x_star,tpd) on the
            # characteristic originating at to
            h_star[j] = (I-Ksat)*(tpd-to)-Ao*np.sqrt(tpd-tp+tca)+Ao*np.sqrt(to-tp+tca)
            # Now find the value of time where h-h* = 0        
            param = [h_star[j],Ksat,tpd,Ao,tp,tca]
            t_star[j] = fsolve(findt_star,10*tpd,param)
            dt = (t_star[j]-tpd)/100.
            t2 = np.arange(tpd, t_star[j]+ dt/10., dt)
            h2 = h_star[j]-Ksat*(t2-tpd)-Ao*np.sqrt(t2-tp+tca)+Ao*np.sqrt(tpd-tp+tca);
            h2[-1]=0
            hpower2 = h2**(2./3)
            hint2 = np.sum((hpower2[:-1]+hpower2[1:])/2.*(t2[1:]-t2[:-1]))+\
                hpower2[0]/2.*(t2[1]-t2[0])+hpower2[-1]/2.*(t2[-1]-t2[-2])
            x_zw[j] = 5./3*Kr_flat*hint2+x_star[j]
            # If the characteristic crosses L, let's find the time of that crossing
            if x_zw[j]>L:        
                param=[tpd,tca,Ksat,x_star[j],Ao,Kr_flat,L,tp,h_star[j]]
                t_cross[j] = fsolve(findt_cross,(tpd*1.00000000001 + t_star[j])/2,param)

                # Now determine the depth of the characteristic crossing L
                h_cross[j]  = h_star[j]-Ksat*(t_cross[j]-tpd)-\
                            Ao*np.sqrt(t_cross[j]-tp+tca)+Ao*np.sqrt(tpd-tp+tca)

                # Now determine the flow associated with this depth
                q_cross[j] = Kr_flat*h_cross[j]**(5./3)
        t_cross = t_cross[t_cross>0]
        q_cross = q_cross[q_cross>0]
        tout = np.hstack((t1, to_out,t_cross))            
        qout = np.hstack((q, q_out,q_cross))
        
    elif tpd<Teq:
        # Up to t = tpd, characteristics originating on t-tp = 0 and x = xo can
        # cross L
        dt = (tpd-tp)/100.
        T = np.arange(round(tp), round(tpd)+dt/10., dt) 
        h = np.ones_like(T)
        q = np.ones_like(T)
        for j in range(len(T)):
            t=T[j]
            h[j]=((I - Ksat)*(t-tp)-Ao*(np.sqrt(t-tp+tca)-np.sqrt(tca)))
            if (h[j])<0:
                h[j]=0.
            q[j]=Kr_flat*h[j]**(5./3)
        q = np.hstack((np.zeros(int(np.round(tp)-1)), q))
        t1 = np.hstack(( np.arange(1,np.round(tp)), T))

        # Find xo for the last catchment that can cross x=L at time < tpd
        dt = (tpd-tp)/100.
        t = np.arange(tp, tpd+dt/10., dt)
        h = ((I - Ksat)*(t-tp)-Ao*(np.sqrt(t-tp+tca)-np.sqrt(tca)))
        h[h<0]=0
        hpower = h**(2./3)
        hint = np.sum((hpower[:-1] + hpower[1:])/2.*(t[1:] - t[:-1]))+\
                hpower[0]/2.*(t[1]-t[0])+hpower[-1]/2.*(t[-1]-t[-2])
        xo = L - 5./3*Kr_flat*hint
        # for xo located between o and xo, characteristics evolve to a point at
        # time tpd, and then decay.
        XO = np.arange(0, xo + xo/1000., xo/100.)  
        x_star = np.zeros_like(XO[:-1])
        h_star = np.zeros_like(XO[:-1])
        t_star = np.zeros_like(XO[:-1])
        x_zw = np.zeros_like(XO[:-1])
        t_cross = np.zeros_like(XO[:-1])
        h_cross = np.zeros_like(XO[:-1])
        q_cross = np.zeros_like(XO[:-1])
        for j in range(len(XO)-1):
            xo = XO[j]
            dt = (tpd-tp)/100.
            t = np.arange(tp, tpd + dt/10., dt)

            h = (I-Ksat)*(t-tp)-Ao*np.sqrt(t-tp+tca)+Ao*np.sqrt(tca)
            h[h<0] = 0.
            hpower = h**(2./3)
            hint = np.sum((hpower[:-1]+hpower[1:])/2.*(t[1:]-t[:-1]))+\
                     + hpower[0]/2.*(t[1]-t[0]) + hpower[-1]/2.*(t[-1]-t[-2])
            x_star[j] = xo + 5./3*Kr_flat*hint

            # So at tpd, the characteristic is at xstar.
            # Find the value of h (h_star) at this point, (x_star,tpd) on the
            # characteristic originating at to
            h_star[j] = (I-Ksat)*(tpd-tp)-Ao*np.sqrt(tpd-tp+tca)+Ao*np.sqrt(tca)

            if h_star[j]>0:
                # Now find the value of time where h-h* = 0
                param=[h_star[j],Ksat,tpd,Ao,tp,tca]
                t_star[j] = fsolve(findt_star,0.95*tpd,param)
                dt = (t_star[j]-tpd)/100.
                t2 = np.arange(tpd, t_star[j]+dt/10., dt)
                h2 = h_star[j]-Ksat*(t2-tpd)-Ao*np.sqrt(t2-tp+tca)+\
                        Ao*np.sqrt(tpd-tp+tca)
                h2[-1]=0
                hpower2 = h2**(2./3)
                hint2 = sum((hpower2[:-1]+hpower2[1:])/2.*(t2[1:]-t2[:-1]))+\
                        hpower2[0]/2.*(t2[1]-t2[0])+hpower2[-1]/2.*(t2[-1]-t2[-2])
                x_zw[j] = 5./3*Kr_flat*hint2+x_star[j]   

                # If the characteristic crosses L, let's find the time of that crossing
                if x_zw[j] > L*1.01: 
                    param=[tpd,tca,Ksat,x_star[j],Ao,Kr_flat,L,tp,h_star[j]]
                    t_cross[j] = fsolve(findt_cross,tpd*1.00000000001,param)

                    # Now determine the depth of the characteristic crossing L
                    h_cross[j] = h_star[j]-Ksat*(t_cross[j]-tpd)-Ao*\
                                np.sqrt(t_cross[j]-tp+tca)+Ao*np.sqrt(tpd-tp+tca)

                    # Now determine the flow associated with this depth
                    q_cross[j] = Kr_flat*h_cross[j]**(5./3)
                else: 
                    h_cross[j]=0
                    q_cross[j]=0
                    x_zw[j]=0


        if len(np.where(x_zw>L*1.01)[0])>=1:
            tout = np.hstack((t1, t_cross))
            qout = np.hstack((q, q_cross))
        else:
            tout=t1 
            qout=q
    return tout,qout
                

