************************************************************************

      parameter ( nn=10000,ntp=20000,nx=410,ny=410,nz=101 )

************************************************************************

      implicit real*8(a-h,o-z)
      
      common/m/      grav,amax,epsh

      common/reals/  precip

      common/reals/  flux, tp, start, finish

      common/reals/  tmax,dt,t,tc,tr,beta, 
     &               xk, dz, dtinfl,
     &               r8fluxin, r8fluxout, r8error, r8newmass

      common/reals/  xni, xnv, pveg, Dveg, dxdum 

      common/array/  tt(0:ntp), prate(0:ntp)

      common/array/  xc(0:nx,0:ny), yc(0:nx,0:ny), zc(0:nx,0:ny), 
     &               dzz(0:nx,0:ny,2), xnc(0:nx,0:ny), 
     &               vegc(0:nx,0:ny)

      common/array/  f(0:nx,0:ny,3,2), q(0:nx,0:ny,3), fdum(3), qs(3) 
    
      common/array/  dh(0:nx,0:ny,2), du(0:nx,0:ny,2), dv(0:nx,0:ny,2)
      
      common/array/  itype(nx,ny,4), ipos(nx,ny,4), inum(nx,ny),
     &               nop(nx,ny,4), kbeg(nx), kend(nx)


      common/array/  sn(0:nx,0:ny,2), cn(0:nx,0:ny,2), 
     &               ds(0:nx,0:ny,2), dxi(0:nx,0:ny,2),
     &               deta(0:nx,0:ny,2), area(0:nx,0:ny), 
     &               zinflmap(nz, ny)

      common/srctrm/ fix(nx,ny,3),
     &               sx(0:nx,0:ny), sy(0:nx,0:ny), 
     &               t0(0:nx,0:ny)

      common/ints/   ncol, np, ne, nbcell, ilim,
     &               nprt, it, itp, nt, iprt, ntinfl,
     &               nrow, iflagUp, iflagDown, isV1, isB1,
     &               isat_flagV,  isat_flagB

      common/ints/   inv,  isetflux, niter, iscouple, iscale,
     &               ifixh, infL,isveg,ifric,ifricB, ifricV

      common/reals/  htop, hbottom, stop_tol, stop_tol0

      common/array/  h(0:nx,0:ny), u(0:nx,0:ny), v(0:nx,0:ny),
     &               hn(nn), un(nn), vn(nn),
     &               hp(0:nx,0:ny), up(0:nx,0:ny), vp(0:nx,0:ny),
     &               eta(0:nx,0:ny),  etap(0:nx,0:ny)
    
      common/array/  alpha(nz,2), theta_S(nz,2), theta_R(nz,2),
     &               r8lambda(nz,2), r8n(nz,2), r8m(nz,2), 
     &               r8Ksat(nz,2)

      common/array/  alpha_b(nz), theta_S_b(nz), theta_R_b(nz),
     &               r8lambda_b(nz), r8n_b(nz), r8m_b(nz), 
     &                r8Ksat_b(nz)

      common/array/  DeltaPlus(nz,nz), DeltaMinus(nz,nz),
     &               APlus(nz,nz), AMinus(nz,nz)
     
      common/array/  r8H(0:nx,0:ny,nz), r8THETA(0:nx,0:ny,nz), 
     &               r8K(0:nx,0:ny,nz),  oldTHETA(0:nx,0:ny,nz)
                
      common/array/  z(nz),  iterations(ntp)   
      