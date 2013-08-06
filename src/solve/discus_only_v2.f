	subroutine approximate(cbound,u,d,numbox,numtim,xdelta,
     $                           tdelta,k1,k2,d1,dk2,cdown)

       IMPLICIT NONE
*...........................................................*
* no transient d1 0 k1 to 0 and k2 to 1
* d1 has to do with decay (set to 0)
* k1 
* k2 
* So what this version does is return dk1 as 0 if it passed or dk2 = 1 if it fails
*
*     DISCUS - Domain of Influence Search 
*            for Convective Unconditional Stability
*
*     Copyright - J.Russell Manson
*
*     DISCUS is a numerical method for solving the
*     transport equations which are used to model 
*     various water quality problems.
*
*     This code uses the DISCUS methodology in a one-dimensional 
*     framework making it suitable for modeling the fate and 
*     transport of contaminants in rivers, estuaries and long 
*     narrow lakes.
*
*     Last modified August 19th, 1998.
*...........................................................*

*======================================
*     Declare variables 
*======================================

      real*8 u, d, k1, k2, d1, dk2, rmu2
      real*8 rk1, rk2, rinit, rmass, rpc
      real*8 c(2500), cp(2500)
      real*8 cd(2500), cdp(2500)
      real*8 cnew(2500), dt(2500)

      real*8 cbound(1500), divb(1500), cdown(1500)
      real*8 ctemp(1500,12), ctemp2(1500,12)
      real*8 cout(1500)

      real*8 canal(2500)
      real*8 cinit(2500)
      real*8 div(2500), diva(2500), divcd(2500), divcda(2500)
      real*8 x(2500), deltax(2500), dtbydx(2500)
      real*8 veloc(2500)
      real*8 area(2501)
      real*8 flow(2500)
      real*8 disper(2501)
	real*8 expdif(2501) 
      real*8 decay(2500)

      real*8 dzex(2500), dzarea(2500)

      real*8 mass1, mass2, marea1, marea2, tsarea1, tsarea2
	real*8 mass1l, mass2l, masscd1, masscd2
      real*8 area_old, dzarea_old, cstar, cdstar
      real*8 ratup, ratdn, divd, divu, xdn, xup
      real*8 tcompu, tcompd
      real*8 disp

      real*8 load(2500,12)

      real*8 cg1(2500)
      real*8 cg10(2500)
      real*8 xdelta, tdelta
      real*8 vel, theta, neta, velc
      real*8 mass,mu,var,csa, decr

      real*8 peclet(2501), exch(2501)
      real*8 rpecmx, rtsum, rpcmx, rt, rtav, uo, xdash
      real*8 d0, dimax, cmax2, dxavg, emax, simtim
      real*8 exp1, exp2, cn5max, areaav, cr, crmax

      real*8 cmass

      real*8 PHI
      real*8 kuwah, minmod
      real*8 lin
      real*8 n3init, n3finl

      real*8 dc(4) 
      real*8 a2(2500),b2(2500),c2(2500),d2(2500)
      real*8 dco(2500) 
      real*8 rp
      real*8 coeff1, coeff2, coeff3, coeff0

      real*8 chkmss(9)

      real*8 p(5), dp(4)
      real*8 xp(5)

      real*8 r1, r2, r3
      real*8 s1, s2, s3
	real*8 dfrac, dzfrac, dx, q, ti, tf, time
 
      real*8 table2

      real*8 diff1, diff4

      integer adv, al, it, i, l, numcow, kmax
      integer limiter, itemp, imax, ib
      integer numtim
      INTEGER SHAPE
	integer igeom, iscen1, iscen2, ibnd, ilin
	integer numbox
      integer ns1, ns2, ns3, ns4, ns5, ns6, ns7
      integer ns8, ns9, numpts
      integer jcntup, jcntdn, j, k, nodes

      CHARACTER*50 str
      character*3 ch1

c      set dk2 to pass
       dk2 = 0.0;
*======================================
*     Problem description and setup
*======================================

c      call mexWarnMsgTxt(u,d,k1,k2)

	dfrac = 1.0
	dzfrac = 1.0

      igeom  = 1
	iscen1 = 1
	iscen2 = 1
	ibnd   = 2
	ilin   = 3
	
	if(igeom.eq.1)then

*=========================================
*     Constant coefficients
*=========================================

	x(1) = 0.0
      do i=1,numbox

	x(i+1) = x(i) + xdelta
      deltax(i) = xdelta
    
      end do 
 
c      q = 0.0621
      q = u
      csa = q / (u+0.000000001)

      disp = dble(d)
      velc = dble(u)

	do i=1,numbox

      veloc(i) = velc
      disper(i) = disp
      
	flow(i) = q
	area(i) = flow(i) / veloc(i)

      k1 = k1 + 0.000000001
      k2 = k2 + 0.000000001

      dzex(i) = k1 
	dzarea(i) = dzex(i) * area(i) / k2

c      print*,veloc(i)

      end do

      nodes = numbox + 1

 	endif
       
       
*======================================
*     Other parameters
*======================================

      theta = 0.5
 
      dx = deltax(1)
c	csa = area(1)
c	disp = disper(1) 

      simtim = tdelta * dble(numtim)
	ti = 0.0


      ns1=2
	ns2=4
	ns3=6
	ns4=8
	ns5=10
	ns6=12
      ns7=14
	ns8=16
	ns9=18
      
	var=154400.0
	mu=200.0
	
      decr = d1
	rk1 = k1
	rk2 = k2


      ti = 2.0
      tf = ti + simtim


*======================================
*     Boundary conditions
*======================================


      if(iscen1.eq.1)then

c	do i=1,numtim

c      cbound(i) = 0.0

c      end do

	elseif(iscen1.eq.2)then

      call bound(cbound,numtim,tdelta,1.0D0,100.0D0,10.0D0)

c     do ib = 1, numtim 
c     cbound(ib) = 0.0D0 * exp( -1.0D0 * ( dble(ib) - 20.0D0 ) ** 2 )
c     end do

      elseif(iscen1.eq.3)then

      open(15,file='Rhine.trace')

      do i=1,297

      read(15,*) (ctemp(i,j),j=1,12)

      end do

      close(15)

      do i=1,296

      cbound(i) = 0.5 * ( ctemp(i,1) + ctemp(i+1,1) ) 

      cout(i) = 0.5 * ( ctemp(i,7) + ctemp(i+1,7) ) 

      end do

      
      do i=297,numtim
	cbound(i) = 0.0
	cout(i) = 0.0
      end do

      elseif(iscen1.eq.4)then

      open(15,file='Uvas.trace')

	read(15,*) numpts

c	print*, numpts

      do i=1,numpts

      read(15,*) (ctemp(i,j),j=1,5)

c	print*, i, ctemp(i,1)

      a2(i) = ctemp(i,1)
	b2(i) = ctemp(i,2)

      end do

      close(15)

      do i=1,numtim+1

c      print*,i

      ctemp2(i,1) = dble(i-1) * (tdelta/3600.0) + ctemp(1,1)

	ctemp2(i,2) = table2(numpts,a2,b2,ctemp2(i,1))

c	print*,ctemp2(i,1),ctemp2(i,2)
	
      end do

c      pause
      
      do i=1,(numtim-1)
      cbound(i) = 0.5 * ( ctemp2(i,2) + ctemp2(i+1,2) ) 
      cout(i) = 0.5 * ( ctemp2(i,5) + ctemp2(i+1,5) ) 
      end do

C      do i=154,numtim
C      	cbound(i) = 3.7
C	cout(i) = 3.7
C      end do

	elseif(iscen1.eq.5)then

C      open(15,file='Hart.trace')
C      do i=1,47
C      read(15,*) (ctemp(i,j),j=1,3)
C      end do
C      close(15)
C      do i=1,46
C      cbound(i) = 0.5 * ( ctemp(i,2) + ctemp(i+1,2) ) 
C      cout(i) = 0.5 * ( ctemp(i,3) + ctemp(i+1,3) ) 
C      end do
C      do i=47,numtim
C      cbound(i) = 0.0
C      cout(i) = 0.0
C      end do

      open(15,file='Hart.trace')
	open(25,file='Hart.trace.new')

	read(15,*) numpts

c	print*, numpts

      do i=1,numpts
c	print*,i

      read(15,*) (ctemp(i,j),j=1,3)

c	print*, i, ctemp(i,1)

      a2(i) = ctemp(i,1)
	b2(i) = ctemp(i,2)
	c2(i) = ctemp(i,3)

      end do

      close(15)

      do i=1,numtim+1

c      print*,i

      ctemp2(i,1) = dble(i-1) * tdelta + ctemp(1,1)

	ctemp2(i,2) = table2(numpts,a2,b2,ctemp2(i,1))

	ctemp2(i,3) = table2(numpts,a2,c2,ctemp2(i,1))

c	print*,ctemp2(i,1),ctemp2(i,2)
c	write(25,*) ctemp2(i,1),ctemp2(i,2),ctemp2(i,3)
	
      end do
      
      do i=1,(numtim-1)
      cbound(i) = 0.5 * ( ctemp2(i,2) + ctemp2(i+1,2) ) 
      cout(i) = 0.5 * ( ctemp2(i,3) + ctemp2(i+1,3) ) 
      end do

	elseif(iscen1.eq.6)then

      do i=1,numtim
	 if(i.eq.3)then
	 cbound(i) = 2.894
	 else
	 cbound(i) = 0.0
	 endif
      end do

      elseif(iscen1.eq.7)then

      open(15,file='cow.trace')

      read(15,*) numcow

      do i=1,numcow

      read(15,*) (ctemp(i,j),j=1,4)

      end do

      close(15)

      do i=1,(numcow-1)

      cbound(i) = 0.5 * ( ctemp(i,2) + ctemp(i+1,2) ) 

      cout(i) = 0.5 * ( ctemp(i,4) + ctemp(i+1,4) ) 

      end do
	
      do i=numcow,numtim
	cbound(i) = 0.0
	cout(i) = 0.0
      end do

	endif

      

      divb(1) = 0.0

      do 50 ib=1,numtim

      divb(ib+1) = divb(ib) + tdelta * cbound(ib) * area(1)

c      print*, ib, cbound(ib), divb(ib)

   50 continue

c      pause


      open(unit=10,file='discus.output')

      open(unit=20,file='mass.wsp')

      open(11,file='discus.boundary')

*======================================
*     Initial conditions
*======================================


      rmu2 = mu + simtim * velc

      mass = 250000.0

      time = ti

      tf = time + dble(numtim) * tdelta

*==================================
*     Initial Conditions
*==================================


      if(ibnd.eq.0)then

      if(iscen2.eq.1)then

*======================================
*     Specify initial variance
*======================================

      call initpw(c,numbox,var,mu,mass,csa,dx)
      call initpw(cinit,numbox,var,mu,mass,csa,dx)
      call initpw(canal,numbox,var,rmu2,mass,csa,dx)

	elseif(iscen2.eq.2)then

*======================================
*     Specify initial time-dispersion
*======================================

      call analpw(c,numbox,disp,mu,mass,csa,dx,ti)
      call analpw(cinit,numbox,disp,mu,mass,csa,dx,ti)
      call analpw(canal,numbox,disp,rmu2,mass,csa,dx,tf)
  
  	elseif(iscen2.eq.3)then

*======================================
*     Specify initial condition for Zoppou
*      - pure advection
*======================================
    
      call zopp1t(c,numbox,0.2D0,0.04D0,dx,0.0D0)
      call zopp1t(cinit,numbox,0.2D0,0.04D0,dx,0.0D0)
      call zopp1t(canal,numbox,0.2D0,0.04D0,dx,tf)

	elseif(iscen2.eq.4)then

*======================================
*     Specify initial condition for Zoppou
*      - advection-diffusion
*======================================

      call zoppad(c,numbox,uo,d0,dx,ti)
      call zoppad(cinit,numbox,uo,d0,dx,ti)
      call zoppad(canal,numbox,uo,d0,dx,tf)

      elseif(iscen2.eq.5)then
 
*======================================
*     Specify initial condition as box
*======================================

	do i=1,numbox
	c(i) = 0.0
	cd(i) = 0.0D0 
	end do      

*      imax = int((50.0/deltax(1))+0.5)

      imax = 1
      do i=1,imax
C	c(i) = 100000000.0D0 / ( 80.0D0 * deltax(i) * 1000.0D0 * tdelta)
	c(i) = 1.0D0 

      end do

      endif

	else

	do i=1,numbox
	c(i) = 0.0D0
	cd(i) = 0.0D0 
	end do

      endif

      crmax = -1000.0
      dimax = -1000.0
      rpecmx = -1000.0

      do 100 i=1,numbox

      decay(i) = decr

	dtbydx(i) = tdelta / deltax(i)

	peclet(i) = area(i) * disper(i) * tdelta / (deltax(i)*deltax(i))

       diff1 = disper(i) * tdelta / deltax(i)**2 
       diff4 = area(i) * disper(i) * tdelta / deltax(i)**2 
	    cr = veloc(i) * tdelta / deltax(i)
	   rpc = cr / diff1

c      print*,i, veloc(i), area(i), disper(i), peclet(i)

      cp(i) = 0.0
      cnew(i) = 0.0

      if(cr.gt.crmax)then
      crmax = cr
      endif

      if(diff1.gt.dimax)then
      dimax = diff1
      endif

	if(diff4.gt.rpecmx)then
      rpcmx = diff4
      endif

  100 continue


c      print*, crmax, dimax, rpcmx
	
c      pause

*======================================
*     Initialise certain variables
*======================================

      rtsum = 0.0D0 

	do i=1,9
	chkmss(i) = 0.
	end do

c      write(11,'(11f10.3)') 0.00, cbound(1), c(ns1), c(ns2) 
c     $                      ,c(ns3), c(ns4), c(ns5), c(ns6),
c     $                      c(ns7), c(ns8), c(ns9) 

      cdown(1) = c(numbox)

*======================================
*     Start time stepping
*======================================

      cn5max = -100.

      time = ti

      ch1 = char(numtim)

       write(20,*) 'Got here'
       write(20,*) u,d,k1,k2,numbox,numtim,xdelta,tdelta

c      call mexWarnMsgTxt(ch1)

      do 250 it=1,numtim

c      print*, it, numbox

      time = time + dble(it) * tdelta

*=========================================
*     Construct discrete integral variable
*=========================================

      div(1) = 0.0
      diva(1) = 0.0
      divcd(1) = 0.0
      divcda(1) = 0.0
	
      do 900 i=2,nodes

      div(i) = div(i-1) +  c(i-1) * area(i-1) * deltax(i-1)
      diva(i) = diva(i-1) +  area(i-1) * deltax(i-1)
      divcda(i) = divcda(i-1) +  dzarea(i-1) * deltax(i-1)
      divcd(i) = divcd(i-1) +  cd(i-1) * dzarea(i-1) * deltax(i-1)
      
  900 continue

      emax = -100.

	expdif(1) = disper(1) * 0.0

	exch(1) = 0.0

      do 901 i=2,nodes-1

      areaav = 0.5 * ( area(i-1) + area(i) )

	dxavg = 0.5 * ( deltax(i-1) + deltax(i) )

	exch(i) = areaav * disper(i) / dxavg

	expdif(i) =   (1.0-theta) 
     $            * tdelta * exch(i)
     $            * ( c(i) - c(i-1) )
     
      if(expdif(i).gt.emax)then
	emax = expdif(i)
	endif	 
  
  901 continue

C      print*, emax

      expdif(nodes) = disper(nodes) * 0.0
	exch(nodes) = 0.0

       write(20,*) 'Got here'


*======================================
*     Set dirichlet boundary conditions
*======================================

      tcompu = tdelta
      jcntup = 1
      xup = x(2)

  910 if(tcompu.gt.(deltax(1)/veloc(1)))then

      tcompu = tcompu - deltax(1) / veloc(1) 

c     jcntup = jcntup - 1

c     xup = xup - xdelta

c     On temporal domain

      ratup = tcompu / tdelta

         if(it.eq.1.or.it.eq.numtim)then

         divu = lin(divb(it),divb(it+1),ratup) 

         cp(1) = ( divb(it+1) - divu ) / ( tdelta - tcompu )
		 
	   dt(1) = tdelta - 0.5 * ( tcompu )
	   dtbydx(1) = dt(1) / deltax(1)

         else

         p(1) = divb(it-1)
         p(2) = divb(it)
         p(3) = divb(it+1)
         p(4) = divb(it+2)

         xp(1) = -tdelta
         xp(2) = 0.0
         xp(3) = tdelta
         xp(4) = 2.0*tdelta

         call poly2(p,xp,ratup*tdelta,divu)

C	   divu = lin(p(2),p(3),ratup)

         cp(1) = ( divb(it+1) - divu ) / ( tdelta - tcompu ) 

	   dt(1) = tdelta - 0.5 * ( tcompu )
	   dtbydx(1) = dt(1) / deltax(1)

         endif

      else

c     On spatial domain

         xdash = veloc(1) * tcompu

         xup = xup - xdash

         ratup = 1.0D0 - xdash / deltax(1)

         divu = lin(div(jcntup),div(jcntup+1),ratup) 

         cp(1) = ( (divu-div(1)) 
     $              + (divb(it+1)-divb(it)) * (xdash/tdelta) ) 
     $            / deltax(1)

         dt(1) = tdelta
      
	endif


      

*======================================
*     Compute first estimate of future 
*     values of dependent variable 
*======================================


      do 400 j=2,numbox

       write(20,*) 'Got here'

      tcompu = tdelta
      jcntup = j
      xup = x(jcntup+1)

  980 if(tcompu.gt.(deltax(jcntup)/veloc(jcntup)).and.jcntup.gt.0)then

c         print*,jcntup,deltax(jcntup),veloc(jcntup)
c     $    ,deltax(jcntup)/veloc(jcntup),tcompu
      
      tcompu = tcompu - deltax(jcntup) / veloc(jcntup) 

      jcntup = jcntup - 1

c         print*,jcntup,deltax(jcntup),veloc(jcntup)
c     $    ,deltax(jcntup)/veloc(jcntup),tcompu

      if(jcntup.eq.0)then
	  xup = 0.0
	  else
      xup = xup - deltax(jcntup)
      endif

      goto 980

      else

         if(jcntup.gt.0)then
         
	   xdash = veloc(jcntup) * tcompu

         xup = xup - xdash
	   
	   endif

      endif

      if(jcntup.gt.0)then
       ratup = 1.0D0 - xdash / deltax(jcntup)
      endif
 
     
      if(jcntup.eq.1.or.jcntup.eq.numbox)then

      mass1 = lin(div(jcntup),div(jcntup+1),ratup) 

      marea1 = lin(diva(jcntup),diva(jcntup+1),ratup) 

      tsarea1 = lin(divcda(jcntup),divcda(jcntup+1),ratup) 

      masscd1 = lin(divcd(jcntup),divcd(jcntup+1),ratup) 

      elseif(jcntup.ge.2)then

      p(1) = div(jcntup-1)
      p(2) = div(jcntup)
      p(3) = div(jcntup+1)
      p(4) = div(jcntup+2)

      xp(1) = -deltax(jcntup-1)
      xp(2) = 0.0
      xp(3) = deltax(jcntup)
      xp(4) = deltax(jcntup) + deltax(jcntup+1)

      dp(1) = expdif(jcntup-1)
      dp(2) = expdif(jcntup)
      dp(3) = expdif(jcntup+1)
      dp(4) = expdif(jcntup+2)


         if(ilin.eq.3)then

         limiter = 1
         call polylim2(p,xp,ratup*deltax(jcntup),mass1,limiter)

	   elseif(ilin.eq.1)then

         mass1 = lin(div(jcntup),div(jcntup+1),ratup)
 
         else

         call poly2(p,xp,ratup*deltax(jcntup),mass1)

         endif

C 	call poly2e(p,xp,coeff0,coeff1,coeff2,coeff3)

C 	mass1 = poly2c(ratup*deltax(jcntup),coeff0,coeff1,coeff2,coeff3)

C 	call poly2(dp,xp,ratup*deltax(jcntup),exp1)

 	exp1 = lin(dp(2),dp(3),ratup)

      marea1 = lin(diva(jcntup),diva(jcntup+1),ratup) 

      tsarea1 = lin(divcda(jcntup),divcda(jcntup+1),ratup) 

      masscd1 = lin(divcd(jcntup),divcd(jcntup+1),ratup)

      else

c     Temporal domain

         ratup = tcompu / tdelta

cc         if(it.eq.1.or.it.eq.numtim)then

         divu = lin(divb(it),divb(it+1),ratup) 

         mass1 = ( divu - divb(it) ) 

cc         else

cc         p(1) = divb(it-1)
cc         p(2) = divb(it)
cc         p(3) = divb(it+1)
cc         p(4) = divb(it+2)

cc         xp(1) = -tdelta
cc         xp(2) = 0.0
cc         xp(3) = tdelta
cc         xp(4) = 2.0*tdelta

cc         call poly2(p,xp,ratup*tdelta,divu)

C		   divu = lin(p(2),p(3),ratup)


         mass1 = ( divu - divb(it) ) 
  
cc         endif

      endif


      tcompd = tdelta
      jcntdn = j-1
      xdn = x(jcntdn+1)

  990 if(tcompd.gt.(deltax(jcntdn)/veloc(jcntdn)).and.jcntdn.gt.0)then

      tcompd = tcompd - deltax(jcntdn) / veloc(jcntdn) 

      jcntdn = jcntdn - 1

      if(jcntdn.eq.0)then
	xdn = 0.0
	else
      xdn = xdn - deltax(jcntdn)
      endif

      goto 990

      else

         if(jcntdn.gt.0)then

         xdash = veloc(jcntdn) * tcompd

         xdn = xdn - xdash
	   
	   endif

      endif

      if(jcntdn.gt.0)then
      ratdn = 1.0D0 - xdash / deltax(jcntdn)
      endif

     
      if(jcntdn.eq.1.or.jcntdn.ge.(numbox-1))then

      mass2 = lin(div(jcntdn),div(jcntdn+1),ratdn) 

      marea2 = lin(diva(jcntdn),diva(jcntdn+1),ratdn) 

      tsarea2 = lin(divcda(jcntdn),divcda(jcntdn+1),ratdn) 

      masscd2 = lin(divcd(jcntdn),divcd(jcntdn+1),ratdn)


      elseif(jcntdn.ge.2)then

      p(1) = div(jcntdn-1)
      p(2) = div(jcntdn)
      p(3) = div(jcntdn+1)
      p(4) = div(jcntdn+2)

	xp(1) = -deltax(jcntdn-1)
      xp(2) = 0.0
      xp(3) = deltax(jcntdn)
      xp(4) = deltax(jcntdn) + deltax(jcntdn+1)

      dp(1) = expdif(jcntdn-1)
      dp(2) = expdif(jcntdn)
      dp(3) = expdif(jcntdn+1)
      dp(4) = expdif(jcntdn+2)

       
	  if(ilin.eq.3)then

	         limiter = 1
         call polylim2(p,xp,ratdn*deltax(jcntdn),mass2,limiter)

	  elseif(ilin.eq.2)then

	  p(1) = div(jcntdn-1)
        p(2) = div(jcntdn)
        p(3) = div(jcntdn+1)
        p(4) = div(jcntdn+2)
	  p(5) = div(jcntdn+3)

	  xp(1) = -deltax(jcntdn-1)
        xp(2) = 0.0
        xp(3) = deltax(jcntdn)
        xp(4) = deltax(jcntdn) + deltax(jcntdn+1)
	  xp(5) = deltax(jcntdn) + deltax(jcntdn+1) + deltax(jcntdn+2)

        call poly5(p,xp,ratdn*deltax(jcntdn),mass2)

        call poly5(p,xp,deltax(jcntdn)+ratup*deltax(jcntdn+1),mass1)

	  elseif(ilin.eq.1)then

        mass2 = lin(div(jcntdn),div(jcntdn+1),ratdn) 
    
        else

        call poly2(p,xp,ratdn*deltax(jcntdn),mass2)

        endif


C 	call poly2e(p,xp,coeff0,coeff1,coeff2,coeff3)
C 	mass2 = poly2c(ratup*deltax(jcntdn),coeff0,coeff1,coeff2,coeff3)

C     call poly2(dp,xp,ratup*deltax(jcntdn),exp2)

	exp2 = lin(dp(2),dp(3),ratdn)

      marea2 = lin(diva(jcntdn),diva(jcntdn+1),ratdn) 

      tsarea2 = lin(divcda(jcntdn),divcda(jcntdn+1),ratdn) 

      masscd2 = lin(divcd(jcntdn),divcd(jcntdn+1),ratdn)

      else

c     Temporal domain

         ratup = tcompd / tdelta

cc         if(it.eq.1.or.it.eq.numtim)then

         divd = lin(divb(it),divb(it+1),ratup) 

         mass2 = ( divd - divb(it) )

cc         else 

cc         p(1) = divb(it-1)
cc         p(2) = divb(it)
cc         p(3) = divb(it+1)
cc         p(4) = divb(it+2)

cc         xp(1) = -tdelta
cc         xp(2) = 0.0
cc         xp(3) = tdelta
cc         xp(4) = 2.0*tdelta

cc         call poly2(p,xp,ratup*tdelta,divd)

C		   divd = lin(p(2),p(3),ratup)


cc         mass2 = ( divd - divb(it) ) 
  
cc         endif

      endif


c===========================
c      Accounting for mass
c===========================


       if(jcntup.ge.1.and.jcntdn.ge.1)then

	 dt(j) = tdelta

	 dtbydx(j) = dt(j) / deltax(j)

       elseif(jcntup.eq.1.and.jcntdn.eq.0)then

	 dt(j) = tdelta

	 dtbydx(j) = dt(j) / deltax(j)

	 elseif(jcntup.eq.0.and.jcntdn.eq.0)then

	 dt(j) = tdelta - 0.5 * ( tcompu + tcompd )

	 dtbydx(j) = dt(j) / deltax(j)

       endif


       if(jcntup.ge.1.and.jcntdn.ge.1)then

       cp(j) = ( mass1 - mass2 ) / deltax(j) 

 	 cp(j) = cp(j) + (1.0-theta) * ( exp1 - exp2 ) / ( xup - xdn )

c
c     Note that decay should be evaluated at foot of characteristic
c

       area_old = (marea1 - marea2) / ( xup - xdn )
       dzarea_old = (tsarea1 - tsarea2) / ( xup - xdn )
       cstar = ( ( mass1 - mass2 ) / ( xup - xdn ) ) / area_old
       cdstar = ( ( masscd1 - masscd2 ) / ( xup - xdn ) ) / dzarea_old


	 cp(j) = cp(j) - (1.0-theta) * decay(j) * tdelta * 
     $                  ( mass1 - mass2 ) / ( xup - xdn )
     $  + (1.0-theta) * dzex(j) * area_old * tdelta * ( cdstar - cstar ) 

       elseif(jcntup.eq.1.and.jcntdn.eq.0)then

       cp(j) = ((mass1-div(1))+(mass2*veloc(j)/deltax(j))) / deltax(j)

       elseif(jcntup.eq.0.and.jcntdn.eq.0)then

       cp(j) = ( mass2 - mass1 ) / ( tcompd - tcompu )

       endif

c       print*,mass1,mass2,cp(j),area(j)

c        print*,cp(j)

  400 continue


c.....DIFFUSION PART


        if(theta.gt.0.49)then

c             call mexWarnMsgTxt('Still going diff')

c     	       call diff3(cp,exch,dtbydx,area,theta,dt,
c     $         	 decay,dzex,dzarea,numbox,cd,cdp,cnew)

     	       call diff5(cp,exch,dtbydx,area,theta,dt,
     $         	 decay,dzex,dzarea,numbox,cd,cdp,cnew)


c	       do itemp=1,numbox
c 	       cnew(itemp) = cp(itemp) 
c	       end do

  	else

	       do itemp=1,numbox
 	       cnew(itemp) = cp(itemp) / area(itemp)
	       end do
      
  	endif



      IF(IT.EQ.numtim)THEN

      do 450 l=1,numbox

      write(10,449)0.5*(x(l)+x(l+1)),cinit(l),canal(l),cnew(l)

  449 format(5f15.8)

  450 continue

      ENDIF

c      print*, numbox

      rmass = cmass(cnew,area,deltax,numbox)

c      print*, numbox

      if(it.eq.1)then
      rinit = rmass
      endif

      cmax2 = -10000.
      do 500 k=1,numbox

C     c(k) = cp(k)

      c(k) = cnew(k)
	cd(k) = cdp(k)

      if(cnew(k).gt.cmax2)then
      cmax2=cnew(k)
	kmax = k
      endif

  500 continue

c      print*,kmax,cmax2

c       print*,numbox

      write(20,'(i4,2f15.3,i5,f15.3)') it,rmass,100.0*(rmass/rinit),
     $	kmax,cmax2

c      cg1(it) = cnew(n1)/area(n1)
c      cg10(it) = cnew(n2)/area(n2)

****  call comp(canal,c,numbox,rt)

      rtsum = rtsum + rt

cc    write(11,'(5f15.6)') dble(it-1), cbound(it), cp(1), cp(5), cp(150) 

CCC   write(11,'(5f15.6)') dble(it-1), cbound(it), cout(it), cnew(num2)

      write(11,'(11f10.3)') dble(it-1), cbound(it), cnew(ns1), cnew(ns2) 
     $                      ,cnew(ns3), cnew(ns4), cnew(ns5), cnew(ns6),
     $                      cnew(ns7), cnew(ns8), cnew(ns9) 

      if(cnew(ns5).gt.cn5max)then
      cn5max = cnew(ns5)
	endif

	chkmss(1) = chkmss(1) + ( cnew(ns1) * flow(ns1) )
      chkmss(2) = chkmss(2) + ( cnew(ns2) * flow(ns2) )
      chkmss(3) = chkmss(3) + ( cnew(ns3) * flow(ns3) )
      chkmss(4) = chkmss(4) + ( cnew(ns4) * flow(ns4) )
      chkmss(5) = chkmss(5) + ( cnew(ns5) * flow(ns5) )
      chkmss(6) = chkmss(6) + ( cnew(ns6) * flow(ns6) )
      chkmss(7) = chkmss(7) + ( cnew(ns7) * flow(ns7) )
      chkmss(8) = chkmss(8) + ( cnew(ns8) * flow(ns8) )
      chkmss(9) = chkmss(9) + ( cnew(ns9) * flow(ns9) )

C      write(11,'(11f10.3)') dble(it-1), cbound(it), cnew(ns2), cnew(ns3)

      if(cmax2.gt.1000000.)then
	print*,'Computations terminated due to instability'
	dk2 = 1.0	
c	stop
	return
	endif

c      print*,numbox

      cdown(it+1) = cnew(numbox)

c      numbox = 10

c      print*,it,numbox,cdown(it+1)

  250 continue

      
c	print*, 'Answer', cn5max


      rtav = rtsum / dble(numtim)

c	write(6,*) (chkmss(i),i=1,9)

	write(20,'(a5)') '====='

	write(20,'(9f8.2)') (chkmss(i),i=1,9)

      write(20,*) 'Got to end'

	close(11)
	close(20)

**    rtav = rt

      end

      real*8 function lin(udn,uup,rat)
      real*8 udn,uup,rat
      lin = rat * uup + (1.0-rat)*udn
      return
      end

      subroutine poly(p,dpbydx,x,xc,pcomp)
c
c
c
      real*8 p(2)
      real*8 dpbydx(2)
      real*8 x(2)
      real*8 xc
      real*8 pcomp
      real*8 ab(4,5)
      integer ip(2,2)
c
      irow = 0
c
      do 100 i=1,2
c
      do 50 k=1,2
c
      if(k.eq.1)then
      irow = irow + 1
      ab(irow,1) = 1.0
      ab(irow,2) = x(i)
      ab(irow,3) = x(i)*x(i)
      ab(irow,4) = x(i)*x(i)*x(i)
      ab(irow,5) = p(i)
      elseif(k.eq.2)then
      irow = irow + 1
      ab(irow,1) = 0.0
      ab(irow,2) = 1.0
      ab(irow,3) = 2.0*x(i)
      ab(irow,4) = 3.0*x(i)*x(i)
      ab(irow,5) = dpbydx(i)
      endif
   50 continue
  100 continue

c     call matrix solver

      call elim(ab,4,5,4)

       pcomp =
     $   ab(1,5) * 1.0D0 +
     $   ab(2,5) * xc +
     $   ab(3,5) * xc*xc +
     $   ab(4,5) * xc*xc*xc 

      return
      end

      subroutine poly2(p,x,xc,pcomp)

c

c

c

      real*8 p(4)

      real*8 x(4)

      real*8 xc

      real*8 pcomp

      real*8 ab(4,5)

c

      irow = 0

c

      do 100 i=1,4

c

c

      irow = irow + 1

      ab(irow,1) = 1.0

      ab(irow,2) = x(i)

      ab(irow,3) = x(i)*x(i)

      ab(irow,4) = x(i)*x(i)*x(i)

      ab(irow,5) = p(i)



  100 continue



c     call matrix solver



      call elim(ab,4,5,4)



      pcomp =

     $   ab(1,5) * 1.0D0 +

     $   ab(2,5) * xc +

     $   ab(3,5) * xc*xc +

     $   ab(4,5) * xc*xc*xc 

       

      return

      end


	subroutine poly6(p,x,xc,pcomp)

c
c
c

      real*8 p(6)

      real*8 x(6)

      real*8 xc

      real*8 pcomp

      real*8 ab(6,7)

c

      irow = 0

c

      do 100 i=1,6

c
c

      irow = irow + 1

      ab(irow,1) = 1.0

      ab(irow,2) = x(i)

      ab(irow,3) = x(i) * ab(irow,2)

      ab(irow,4) = x(i) * ab(irow,3) 

	ab(irow,5) = x(i) * ab(irow,4) 

	ab(irow,6) = x(i) * ab(irow,5)


      ab(irow,7) = p(i)



  100 continue



c     call matrix solver



      call elim(ab,6,7,6)



      pcomp =

     $   ab(1,7) * 1.0D0 +

     $   ab(2,7) * xc +

     $   ab(3,7) * xc*xc +

     $   ab(4,7) * xc*xc*xc +
     
     $   ab(5,7) * xc*xc*xc*xc + 

     $   ab(6,7) * xc*xc*xc*xc*xc

       

      return

      end

      subroutine poly5(p,x,xc,pcomp)

c
c
c

      real*8 p(5)

      real*8 x(5)

      real*8 xc

      real*8 pcomp

      real*8 ab(5,6)

c

      irow = 0

c

      do 100 i=1,5

c
c

      irow = irow + 1

      ab(irow,1) = 1.0

      ab(irow,2) = x(i)

      ab(irow,3) = x(i) * x(i)

      ab(irow,4) = x(i) * x(i) * x(i)

	ab(irow,5) = x(i) * x(i) * x(i) * x(i) 

      ab(irow,6) = p(i)



  100 continue



c     call matrix solver



      call elim(ab,5,6,5)



      pcomp =

     $   ab(1,6) * 1.0D0 +

     $   ab(2,6) * xc +

     $   ab(3,6) * xc*xc +

     $   ab(4,6) * xc*xc*xc +
     
     $   ab(5,6) * xc*xc*xc*xc 

       

      return

      end

      subroutine poly2e(p,x,c0,c1,c2,c3)

c

c

c

      real*8 p(4)

      real*8 x(4)

      real*8 ab(4,5)

	real*8 c0, c1, c2, c3

c

      irow = 0

c

      do 100 i=1,4

c

c

      irow = irow + 1

      ab(irow,1) = 1.0

      ab(irow,2) = x(i)

      ab(irow,3) = x(i)*x(i)

      ab(irow,4) = x(i)*x(i)*x(i)

      ab(irow,5) = p(i)



  100 continue



c     call matrix solver



      call elim(ab,4,5,4)



      c0 = ab(1,5)
	c1 = ab(2,5)
	c2 = ab(3,5)
	c3 = ab(4,5)

      return

      end

      real*8 function poly2c(xc,c0,c1,c2,c3)

c

c

c

      real*8 xc, c0, c1, c2, c3

c
      poly2c =

     $   c0 * 1.0D0 +

     $   c1 * xc +

     $   c2 * xc*xc +

     $   c3 * xc*xc*xc 

       

      return

      end



      subroutine elim(ab,n,np,ndim)

c

c-------------------------------------------------------------------

c

c     subroutine elim:

c                     This subroutine solves a set of linear 

c     equations and gives an LU decomposition of the coefficient

c     matrix. The gauss elimimation method is used, partial pivoting

c     is used aswell. Multiple right hand sides are permitted, they

c     should be supplied as columns that augment the coefficient

c     matrix.

c

c-------------------------------------------------------------------

c

c     Parameters are:

c

c     AB    - coefficient matrix augmented with R.H.S. vectors

c     N     - number of equations

c     NP    - total number of columns in the augmented matrix

c     NDIM  - first dimension of matrix AB in the calling program

c

c     The solution vector(s) are returned in the augmentation columns

c     of AB.

c

c-------------------------------------------------------------------

c

      real*8 ab(ndim,np)

      integer n,np,ndim

      real*8 save,ratio,value

      integer nm1,ipvt,ip1,j,nvbl,l,kcol,jcol,jrow

c

c-------------------------------------------------------------------

c

c     Begin the reduction

c

      nm1 = n - 1

      do 35 i = 1,nm1

c

c     Find the row number of the pivot row. We will then interchange

c     rows to put the pivot element on the diagonal. 

c

      ipvt = i

      ip1 = i + 1

      do 10 j = ip1,n

         if(abs(ab(ipvt,i)).lt.abs(ab(j,i))) ipvt = j

   10 continue

c

c-------------------------------------------------------------------

c

c     Check for a near singular matrix

c

      if(abs(ab(ipvt,i)).lt.1.0E-21)then

      print 100

      return

      endif

c

c     Now, interchange except if the pivot element is already on

c     diagonal - don't need to.

c

c

      if(ipvt.ne.i)then

         do 20 jcol=1,np

         save = ab(i,jcol)

         ab(i,jcol) = ab(ipvt,jcol)

         ab(ipvt,jcol) = save

   20 continue

      endif

c

c-------------------------------------------------------------------

c

c     Now reduce all elements below the diagonal in the i-th row

c     Check first to see if a zero already present. If so, can skip

c     reduction on that row.     

c

      do 32 jrow=ip1,n

        if(ab(jrow,i).eq.0) goto 32

           ratio = ab(jrow,i) / ab(i,i)

           ab(jrow,i) = ratio

           do 30 kcol = ip1,np

           ab(jrow,kcol) = ab(jrow,kcol) - ratio * ab(i,kcol)

   30      continue

   32    continue

   35 continue

c

c-------------------------------------------------------------------

c

c     We still need to check ab(n,n) for size

c

      if(abs(ab(n,n)).le.1.0E-21) then

      print 100

      return

      endif

c

c-------------------------------------------------------------------

c

c     Now we back substitute

c

      np1 = n + 1

     

      do 50 kcol=np1,np

        ab(n,kcol) = ab(n,kcol) / ab(n,n)

        do 45 j=2,n

        nvbl = np1 - j

        l = nvbl + 1

        value = ab(nvbl,kcol)

        do 40 k=l,n

           value = value - ab(nvbl,k) * ab(k,kcol)

   40   continue

        ab(nvbl,kcol) = value / ab(nvbl,nvbl)

   45   continue


   50 continue

      return

c

  100 format(/' SOLUTION NOT FEASIBLE. A NEAR ZERO PIVOT ',

     $       'WAS ENCOUNTERED.')

      end


      subroutine initpw(c,numbox,var,mu,mass,csa,dx)

      integer numbox 

      real*8 c(numbox)

      real*8 mu, var, mass, csa

      real*8 movera, u, sigma, sigsqr

      real*8 pi, n3, dx, peak, x, xdn, xup, cdn, cup 



      movera = mass / csa

      sigma = sqrt(var)

      sigsqr = var

      pi = 4.0D0 * atan(1.0)

      n3 = 4.0*sigma/dx

      peak = movera / sqrt(pi*var)



      do 100 i=1,numbox

      

      xdn = dble(i-1)*dx - mu

      xup = dble(i)*dx - mu



      cdn = peak * exp( -1.0*xdn*xdn/sigsqr ) 

      cup = peak * exp( -1.0*xup*xup/sigsqr ) 



      c(i) = 0.5 * ( cdn + cup ) 



  100 continue



      return

      end


      subroutine analpw(c,numbox,disp,mu,mass,csa,dx,time)

      integer numbox 

      real*8 c(numbox)

      real*8 mu, var, mass, csa

      real*8 movera, u, sigma, sigsqr

      real*8 pi, n3, dx, time, peak, x, xdn, xup, cdn, cup 

      real*8 disp



      movera = mass / csa

      sigsqr = 4.0D0 * disp * time

      pi = 4.0D0 * atan(1.0)

      n3 = 4.0*sigma/dx

      peak = movera / sqrt(4.0*pi*disp*time)



      do 100 i=1,numbox

      
      xdn = dble(i-1)*dx - mu

      xup = dble(i)*dx - mu



      cdn = peak * exp( -1.0*xdn*xdn/sigsqr ) 

      cup = peak * exp( -1.0*xup*xup/sigsqr ) 



      c(i) = 0.5 * ( cdn + cup ) 



  100 continue



      return

      end

      subroutine zopp1i(c,numbox,mu,var,dx)

      integer numbox 
      real*8 c(numbox)
      real*8 mu, var, mass, csa
      real*8 movera, u, sigma, sigsqr
      real*8 pi, n3, dx, time, peak, x 
      real*8 disp

      co = 100.

      sigsqr = 2.0D0 * var

      c(1) = 0.0
      do 100 i=2,numbox

      
      xdn = dble(i-1)*dx 
      xup = dble(i)*dx 

      xdn2 = log( xdn / mu )
      xup2 = log( xup / mu )

      cdn = co * exp( -1.0*xdn2*xdn2/sigsqr ) 
      cup = co * exp( -1.0*xup2*xup2/sigsqr ) 

      c(i) = 0.5 * ( cdn + cup ) 

  100 continue

      return

      end

      subroutine zopp1t(c,numbox,mu,var,dx,time)

      integer numbox 
      real*8 c(numbox)
      real*8 mu, var, mass, csa
      real*8 movera, u, sigma, sigsqr
      real*8 pi, n3, dx, time, peak, x 
      real*8 disp

      co = 20.

      sigsqr = 2.0D0 * var

      c(1) = 0.0
      do 100 i=2,numbox

      
      xdn = dble(i-1)*dx 
      xup = dble(i)*dx 

      xdn2 = log( xdn / mu ) - 0.1 * time
      xup2 = log( xup / mu ) - 0.1 * time

      cdn = (co/xdn) * exp( -1.0*xdn2*xdn2/sigsqr ) 
      cup = (co/xup) * exp( -1.0*xup2*xup2/sigsqr ) 

      c(i) = 0.5 * ( cdn + cup ) 

  100 continue

      return

      end

      subroutine zoppad(c,numbox,uo,d0,dx,time)

      integer numbox 
      real*8 c(numbox)
      real*8 mu, var, mass, csa
      real*8 movera, uo, d0, sigma, sigsqr
      real*8 pi, n3, dx, time, peak, x 
      real*8 disp

      co = 0.

      sigsqr = 4. *  d0 * time

      peak = co / sqrt( sigsqr ) 

      c(1) = 0.0
      do 100 i=2,numbox

      
      xdn = dble(i-1)*dx 
      xup = dble(i)*dx 

      xdn2 = log( xdn ) - ( uo - d0 ) * time
      xup2 = log( xup ) - ( uo - d0 ) * time

      cdn = peak * exp( -1.0*xdn2*xdn2/sigsqr ) 
     $           * exp( -1.0D0 * uo * time )
      cup = peak * exp( -1.0*xup2*xup2/sigsqr ) 
     $           * exp( -1.0D0 * uo * time )

      c(i) = 0.5 * ( cdn + cup ) 

  100 continue

      return

      end


      subroutine initsq(c,nodes,hei,i1,i2)

      integer nodes 

      real*8 c(nodes)

      real*8 hei

      integer i1, i2



      do 100 i=1,nodes

      

      if(i.ge.i1.and.i.le.i2)then



      c(i) = hei 



      else



      c(i) = 0.0



      endif



  100 continue



      return

      end


      subroutine trid(sub,diag,sup,b,n)

      integer n, i

      real*8 b(n), diag(n), sub(n), sup(n) 

c    

c     The tridiagonal linear system

c

c     sub(i) * x(i-1) + diag(i) * x(i) + sup(i) * x(i+1) = b(i)  i=1,n

c

c     with (sub(1) taken to be zero and sup(n) taken to be zero)

c     is solved by factorization and substitution.  The factorization

c     is returned in sub, diag, sup and the solution is returned in b.



      if(n.le.1)then

        b(1) = b(1) / diag(1)

        return

      endif 



      do 100 i=2,n

 

         sub(i) = sub(i) / diag(i-1) 

         diag(i) = diag(i) - sub(i)*sup(i-1) 

  100    b(i) = b(i) - sub(i)*b(i-1)

       

      b(n) = b(n) / diag(n)



      do 200 i=n-1,1,-1



  200    b(i) = (b(i) - sup(i)*b(i+1)) / diag(i)

      

      return

      

      end



      real*8 function hercub(x1,x2,y1,y2,dy1,dy2,x)



      real*8 x1,x2

      real*8 y1,y2

      real*8 dy1,dy2

      real*8 x, sigma, xdel, h01, h02, h11, h12



      sigma = (1.0/abs(x2-x1)) * (2.0*x-(x1+x2))

      xdel  = abs(x2-x1)



      h01 = 0.25*(sigma**3-3.0*sigma+2.0)

      h02 = -0.25*(sigma**3-3.0*sigma-2.0) 



      h11 = 0.125*xdel*(sigma**3-sigma**2-sigma+1.0)

      h12 = 0.125*xdel*(sigma**3+sigma**2-sigma-1.0)



      hercub = y1 * h01 + y2 * h02

     $       + dy1 * h11 + dy2 * h12



      return

      end


      double precision function ultim(for,cour,vel,cup,ccn,cdn)

      real*8 for,fornew,cup,ccn,cdn,cour
      real*8 rphif, rphic, y1rc, y2rc, y3rc

      rphif = (for-cup)/((cdn-cup)+0.00000001)

      rphic = (ccn-cup)/((cdn-cup)+0.00000001) 

      
      y1rc = rphic
      y2rc = rphic/cour
      y3rc = 1.0

      if(rphic.lt.0.0)then

         rphif = y1rc

      elseif(rphic.lt.cour.and.rphic.gt.0.0)then

         if(rphif.gt.y1rc.and.rphif.lt.y2rc)then

         else

           if(rphif.lt.y1rc)then

           rphif = y1rc

           elseif(rphif.gt.y2rc)then

           rphif = y2rc

           endif

         endif

 

      elseif(rphic.gt.cour.and.rphic.lt.1.0)then

 

         if(rphif.gt.y1rc.and.rphif.lt.y3rc)then

         else

           if(rphif.lt.y1rc)then

           rphif = y1rc

           elseif(rphif.gt.y3rc)then

           rphif = y3rc

           endif

         endif

 

      elseif(rphic.gt.1.0)then

 

         rphif = y1rc

 

      endif

 

      ultim =  rphif * ( cdn - cup ) + cup 

 

      return

      end




      subroutine diff(c,peclet,theta,nodes,cp)



      real*8 c(nodes), cp(nodes)
      real*8 a(2500),b(2500),c1(2500),d(2500),w(2500)
      real*8 e(2500), f(2500)



      integer rank 



      real*8 peclet



c     -------------

c     SET UP MATRIX

c     -------------



      rank  = 0

 

      do 900 l=1, nodes



c     Test for the type of link and insert

c     appropriate equation coefficients

c     ------------------------------------



         if (l.eq.1) then

            rank = rank + 1

            a(rank) =  0.0

            b(rank) = 1.0

            c1(rank) =  0.0

            d(rank) = 0.0

         elseif (l.eq.nodes) then

            rank = rank + 1

            a(rank) =  0.0

            b(rank) = 1.0

            c1(rank) =  0.0

            d(rank) = c(nodes)

         else

            rank = rank + 1

c...........node -1

            a(rank) = theta * peclet 

c...........node 1

            b(rank) = 1.0D0 + 2.0D0 * theta * peclet 

c...........node +1

            c1(rank) = theta * peclet 

c

C            d(rank) = c(l) + peclet * (1.0-theta) * 

C     $                ( c(l+1) - 2.0*c(l) + c(l-1) ) 

	      d(rank) = c(l) 



         endif 

         

  900 continue



c     ----------------------

c     MATRIX SOLVING SECTION

c     ----------------------



c     Preparation of variables for matrix routines.

c     ---------------------------------------------



         nrank = rank



c     Calls to matrix routines

c     ------------------------



c           L-U Decompostion

c           ----------------



      call gtri(a,b,c1,d,rank,e,f,w,1,0.0D0,0.0D0,1,c(nodes),0.0D0)



c     ---------------------

c     END OF MATRIX SECTION

c     ---------------------



c       Assign values from solution vector to 

c       cp values.

c       -------------------------------------



           do 1000 i = 1,nodes



           cp(i) = w(i)



 1000 continue



      end



      subroutine gtri(a,b,c,d,md,e,f,w,l1,a1,q1,lm,am,qm)

c     general tri-diagonal solver

      real*8 a(md), b(md), c(md), d(md), e(md), f(md), w(md)
      real*8 den, a1, am, q1, qm

      if(l1.eq.1) e(1) = 0.0

      if(l1.eq.1) f(1) = a1

c     above overwritten if lm.ne.2

      if(l1.eq.2) e(1) = 1.0

      if(l1.eq.2) f(1) = -a1

      if(l1.eq.3) e(1) = a1/(a1-1.0)

      if(l1.eq.3) f(1) = q1/(1.0-a1)

      mm=md-1

      do 1 m=2,mm

      den=b(m)-c(m)*e(m-1)

      e(m)=a(m)/den

    1 f(m)=(d(m)+c(m)*f(m-1))/den

      if(lm.eq.1) w(md)=am

      if(lm.eq.2) w(md)=(f(mm)+am)/(1.0-e(mm))

      if(lm.eq.3) w(md)=(f(mm)+qm/am)/((1.0+am)/am-e(mm))

      do 2 mk=1,mm

      m=md-mk

    2 w(m)=e(m)*w(m+1)+f(m)

      return

      end  

        

      real*8 function table2(numpts,xarray,yarray,x)



c     This function returns the value of ordinate 

c     which corresponds to abscisa x.

c     It will only find one value for the ordinate 

c     hence should only be used for single valued 

c     functions of x.

c     Calls to the real*8 function interp.



c     Declare variables

c     -----------------



c     numpts    :the number of elements in the 

c                interpolation table.

c     xarray    :array of independent variable.

c     yarray    :array of dependent variable.

c     x         :the independent variable for 

c                which the corresponding ordinate

c                is sought.

c     table     :the corresponding ordinate.

c     i         :a counter



      integer numpts

      real*8 xarray(numpts)

      real*8 yarray(numpts)

      real*8 x

      real*8 interp

      integer i    



c     If only one entry in interpolation table

c     then trivial solution.

c     ----------------------------------------

 

      if(numpts.eq.1) then

         table2 = yarray(1)

      endif



c     If x is not in the range x(1) to x(numpts)

c     then print error message

c     ------------------------------------------

 

      if(x.lt.xarray(1).or.x.gt.xarray(numpts)) then

         print*, 'Error in interpolation routine - TABLE.FOR'
         print*, 'Independent variable sent to routine is   '
         print*, 'outwith the range in the interpolation table'

         stop

      endif 



c     Loop through table to find which two 

c     points x lies between.

c     ------------------------------------



      do 10 i = 1, ( numpts - 1 )

         if(x.ge.xarray(i).and.x.lt.xarray(i+1)) then



c     Interpolate for dependent variable

c     corresponding to x

c     ----------------------------------



      table2=interp(xarray(i),yarray(i),xarray(i+1),yarray(i+1),x)



        elseif(x.eq.xarray(i+1)) then


c     Interpolate for dependent variable
c     corresponding to x
c     ----------------------------------



      table2=interp(xarray(i),yarray(i),xarray(i+1),yarray(i+1),x)


         endif



   10 continue



      return

      end



      real*8 function interp(x1,y1,x2,y2,x)



c     This function returns the y value corresponding

c     to x ( which lies between x1 and x2 ) by using

c     linear interpolation.



c     Declare variables

c     -----------------



c     x1       :abscisa of the first point.

c     y1       :ordinate of the first point.

c     x2       :abscisa of the second point.

c     y2       :ordinate of the second point.

c     x        :abscisa for which the ordinate 

c               is to be calculated.

c     deltax   :the distance between x1 and x2.

c     t        :the ratio of (x-x1) to deltax.

c     interp   :the ordinate to be calculated.  



      real*8 x1, x2

      real*8 y1, y2

      real*8 x

      real*8 deltax, t



c     Calculate deltax and t

c     ----------------------



      deltax = x2 - x1



      t = (x-x1) / deltax



c     Calculate the ordinate corresponding

c     to x.

c     ------------------------------------



      interp = y1 + t * ( y2 - y1 )



c     Return to calling routine

c     -------------------------



      return

      end



      real*8 function table3(numpts,xarray,yarray,x)



c     This function returns the value of ordinate 

c     which corresponds to abscisa x.

c     It will only find one value for the ordinate 

c     hence should only be used for single valued 

c     functions of x.

c     Calls to the real*8 function interp.



c     Declare variables

c     -----------------



c     numpts    :the number of elements in the 

c                interpolation table.

c     xarray    :array of independent variable.

c     yarray    :array of dependent variable.

c     x         :the independent variable for 

c                which the corresponding ordinate

c                is sought.

c     table     :the corresponding ordinate.

c     i         :a counter



      integer numpts

      real*8 xarray(numpts)

      real*8 yarray(numpts)

      real*8 x

      real*8 interp

      integer i    



c     If only one entry in interpolation table

c     then trivial solution.

c     ----------------------------------------

 

      if(numpts.eq.1) then

         table2 = yarray(1)

      endif



c     If x is not in the range x(1) to x(numpts)

c     then print error message

c     ------------------------------------------

 

      if(x.lt.xarray(1).or.x.gt.xarray(numpts)) then

         print*, 'Error in interpolation routine - TABLE.FOR'

         print*, 'Independent variable sent to routine is   '

         print*, 'outwith the range in the interpolation table'

         stop

      endif 



c     Loop through table to find which two 

c     points x lies between.

c     ------------------------------------



      do 10 i = 1, ( numpts - 1 )

         if(x.ge.xarray(i).and.x.lt.xarray(i+1)) then



c     Interpolate for dependent variable

c     corresponding to x

c     ----------------------------------



      table3 = ( yarray(i+1) - yarray(i) ) / 

     $         ( xarray(i+1) - xarray(i) )



         endif



   10 continue



      return

      end



      subroutine comp(canal,c,nodes,rt)

      integer nodes

      real*8 canal(nodes)

      real*8 c(nodes)

      real*8 rt, sum1, sum2



      sum1 = 0.

      sum2 = 0.



      do 100 i=1,nodes

      

      sum1 = sum1 + ( canal(i) - c(i) ) ** 2

      sum2 = sum2 + c(i) 



  100 continue



      avg = sum2 / dble(nodes)



      sum3 = 0.



      do 200 i=1,nodes

      

      sum3 = sum3 + ( avg - c(i) ) ** 2



  200 continue



      rt = 1.0D0 - sum1 / sum3



      return

      end



      subroutine varcmp(c,nodes,dx,sd)

      integer nodes

      real*8 c(nodes)

      real*8 sd

      real*8 sum0, sum1, sum2, x

      sum0 = 0.

      sum1 = 0.

      sum2 = 0.



      do 100 i=1,nodes

      

      x=dble(i-1)*dx



      sum0 = sum0 + c(i) 

      sum1 = sum1 + c(i) * x

      sum2 = sum2 + c(i) * x ** 2



  100 continue



      avg = sum1 / sum0 



      sd = sqrt( (sum2/sum0) - avg ** 2 )



      return

      end



      real*8 function cmass(c,area,dx,nodes)

      integer nodes
      real*8 c(nodes), area(nodes)
      real*8 dx(nodes)

      cmass = 0.0

      do 100 i=1,nodes

      cmass = cmass + c(i) * area(i) * dx(i)

  100 continue



      end



      subroutine diff2(c,diffus,area,theta,nodes,cp)

c     Remember c is mass per unit length
c     but cp is concentration

      real*8 c(nodes), cp(nodes)

      real*8 a(nodes),b(nodes),c1(nodes),d(nodes),w(nodes)

      real*8 e(nodes), f(nodes)

      integer rank 

      real*8 diffus(nodes+1)

      real*8 area(nodes+1)


c     -------------

c     SET UP MATRIX

c     -------------

      rank  = 0

      do 900 l=1, nodes

c     Test for the type of link and insert

c     appropriate equation coefficients

c     ------------------------------------

         if (l.eq.1) then

            rank = rank + 1

c...........node -1

            a(rank) = 0.0

c...........node 1

            b(rank) = area(l) + theta * diffus(l+1)

c...........node +1

            c1(rank) = -1.0D0 * theta * diffus(l+1)

c

C            d(rank) = c(l) + diffus(l+1) * (1.0-theta) * 
C     $                ( c(l+1) - c(l) ) 

             d(rank) = c(l) 

         elseif(l.ge.2.and.l.le.(nodes-1))then

            rank = rank + 1

c...........node -1

            a(rank) = -1.0D0 * theta * diffus(l)

c...........node 1

            b(rank) = area(l) + theta * diffus(l) + theta * diffus(l+1)

c...........node +1

            c1(rank) = -1.0D0 * theta * diffus(l+1)

c

C            d(rank) = c(l) + diffus(l) * (1.0-theta) * 
C     $                ( c(l-1) - c(l) ) 
C     $                + diffus(l+1) * (1.0-theta) * 
C     $                ( c(l+1) - c(l) ) 

             d(rank) = c(l)  


         elseif (l.eq.nodes) then

            rank = rank + 1

c...........node -1

            a(rank) = -1.0D0 * theta * diffus(l)

c...........node 1

            b(rank) = area(l) + theta * diffus(l) 

c...........node +1

            c1(rank) = 0.0

c

C            d(rank) = c(l) + diffus(l) * (1.0-theta) * 
C     $                ( c(l-1) - c(l) ) 

             d(rank) = c(l) 


         endif 

  900 continue



c     ----------------------

c     MATRIX SOLVING SECTION

c     ----------------------



c     Preparation of variables for matrix routines.

c     ---------------------------------------------



         nrank = rank



c     Calls to matrix routines

c     ------------------------

c           L-U Decompostion

c           ----------------



*     call gtri(a,b,c1,d,rank,e,f,w,1,0.0,0.0,1,c(nodes),0.0)


      call trid(a,b,c1,d,rank)



c     ---------------------

c     END OF MATRIX SECTION

c     ---------------------

c       Assign values from solution vector to 

c       cp values.

c       -------------------------------------



           do 1000 i = 1,nodes

           cp(i) = d(i)

 1000 continue

      end


      subroutine bound(c,numtim,dt,peak,mutime,var)

      integer numtim 

      real*8 c(numtim)

      real*8 mutime, var, peak, dt, tun, tdn, cdn, cup

      do 100 i=1,numtim

      
      tdn = dble(i-1)*dt - mutime

      tup = dble(i)*dt - mutime


      cdn = peak * exp( -1.0*tdn*tdn/var ) 

      cup = peak * exp( -1.0*tup*tup/var ) 



      c(i) = 0.5 * ( cdn + cup ) 



  100 continue



      return

      end

      subroutine diff3(c,exch,dtbydx,area,theta,dt,
     $                	decay,dzex,dzarea,nodes,cd,cdp,cp)

c     Remember c is mass per unit length
c     but cp is concentration

      integer rank, nodes, l  

      real*8 c(nodes), cp(nodes), cd(nodes), cdp(nodes)
      real*8 a(nodes),b(nodes),c1(nodes),d(nodes),w(nodes)
      real*8 e(nodes), f(nodes), dt(nodes)
      real*8 exch(nodes+1)
      real*8 dtbydx(nodes) 
	real*8 area(nodes), decay(nodes)
	real*8 dzex(nodes), dzarea(nodes)
      real*8 x,y,z, theta

c     rk1 is alpha
c     rk2 is areas

c     -------------
c     SET UP MATRIX
c     -------------


c     Test for the type of link and insert
c     appropriate equation coefficients
c     ------------------------------------

            rank = 1
               l = 1

	      x = dzarea(l) + dt(l) * dzex(l) * area(l)
     $	          + dt(l) * 0.0D0 * dzarea(l) 

	      y = cd(l) * dzarea(l) / x

	      z = dt(l) * dzex(l) * area(l) / x


c...........node -1

            a(rank) = 0.0

c...........node 1

            b(rank) = area(l) * ( 1.0D0 + dt(l) * dzex(l)
     $		   - dt(l) * dzex(l) * z  + theta * decay(l) * dt(l) )
     $	  + theta * dtbydx(l) * exch(l+1)

c...........node +1

            c1(rank) = -1.0D0 * theta * dtbydx(l) * exch(l+1)

            d(rank) = c(l) + dzex(l) * area(l) * y * dt(l)




         do l=2, nodes-1

            rank = rank + 1

	      x = dzarea(l) + dt(l) * dzex(l) * area(l)
     $	          + dt(l) * 0.0D0 * dzarea(l) 

	      y = cd(l) * dzarea(l) / x

	      z = dt(l) * dzex(l) * area(l) / x

c...........node -1

            a(rank) = -1.0D0 * theta * dtbydx(l) * exch(l)

c...........node 1

            b(rank) = area(l) * ( 1.0D0 + dzex(l) * dt(l) 
     $                		  - dzex(l) * z * dt(l)
     $                		  + theta * decay(l) * dt(l) )
     $                		  + theta * dtbydx(l) * exch(l) 
     $                          + theta * dtbydx(l) * exch(l+1)

c...........node +1

            c1(rank) = -1.0D0 * theta * dtbydx(l) * exch(l+1)

c

            d(rank) = c(l) + dzex(l) * area(l) * y * dt(l)
		   
		end do    


            l = nodes

            rank = rank + 1

	      x = dzarea(l) + dt(l) * dzex(l) * area(l)
     $	          + dt(l) * 0.0D0 * dzarea(l) 

	      y = cd(l) * dzarea(l) / x

	      z = dt(l) * dzex(l) * area(l) / x

c...........node -1

            a(rank) = -1.0D0 * theta * dtbydx(l) * exch(l)

c...........node 1

            b(rank) = area(l)* ( 1.0D0 + dzex(l) * dt(l) 
     $                		  - dzex(l) * z * dt(l)  
     $            		     + theta * decay(l) * dt(l) ) 
     $	                     + theta * dtbydx(l) * exch(l) 

c...........node +1

            c1(rank) = 0.0

c

            d(rank) = c(l) + dzex(l) * area(l) * y * dt(l)





c     ----------------------

c     MATRIX SOLVING SECTION

c     ----------------------



c     Preparation of variables for matrix routines.

c     ---------------------------------------------



         nrank = rank



c     Calls to matrix routines

c     ------------------------

c           L-U Decompostion

c           ----------------



*     call gtri(a,b,c1,d,rank,e,f,w,1,0.0,0.0,1,c(nodes),0.0)


      call trid(a,b,c1,d,rank)



c     ---------------------

c     END OF MATRIX SECTION

c     ---------------------

c       Assign values from solution vector to 

c       cp values.

c       -------------------------------------



           do i = 1,nodes

            cp(i) = d(i)

            x = dzarea(i) + dt(i) * dzex(i) * area(i)
     $	          + dt(i) * 0.0D0 * dzarea(i) 

	      y = cd(i) * dzarea(i) / x

	      z = dt(i) * dzex(i) * area(i) / x


            cdp(i) = y + z * cp(i)

           end do

      end



	subroutine polylim2(p,x,xc,pcomp,limiter)

c     deltax is constant and stored at x(3) 

c     cubic interpolation with universal flux limiter

c     only work at case 1 and case 2

      real*8 p(4)

      real*8 x(4)

      real*8 xc

      real*8 pcomp

      real*8 ab(4,5)

	integer irow,i

      integer  limiter,flag
	real*8 xposs1,xposs2
	real*8 slope1,slope2,slope3
	real*8 xneed,pneed1,pneed2,temp,xneed1
c	rp
	integer steps
c	rp
	steps=0

c     limiter other than 1 means normal option

	if (limiter.eq.1) then

c     x(3)=tdelta or x(3)=xdelta

c     get the position of A
cc    determine the case 1 or 2  

      slope1=(p(2)-p(1))/x(3)
      slope2=(p(3)-P(2))/x(3)
      slope3=(p(4)-p(3))/x(3)
      
      xposs1=p(2)+xc*slope2
	xposs2=0.0
      xneed = x(3)/2.0
	xneed1= xneed
      pneed1=p(2)
      pneed2=p(3)

	flag=0

	toler = 0.001
c     case 1      

      if (slope1.lt.slope2.and.slope2.lt.slope3) then

	   flag=2
         
c   50   if (abs(pneed1-pneed2).ge.toler) then
   50   if (abs(pneed1-pneed2).ge.toler.and.steps.lt.100) then
c	rp
	  steps=steps+1
          pneed1= p(2)+xneed*slope1
          pneed2= p(3)-(x(3)-xneed)*slope3
           if (pneed1.gt.pneed2) then 
	      xneed1=xneed1/2.0
            xneed=xneed+xneed1
           else
	      xneed1=xneed1/2.0
            xneed=xneed-xneed1 
           endif
          goto 50
	    endif
c	rp
c	print*,steps
	   if ( xc.ge.xneed ) then
         xposs2 = p(3)-(x(3)-xc)*slope3
         else
         xposs2 = p(2)+xc*slope1
         endif


c	   if ( xc.ge.xneed) then
c         xposs2 = p(3)-(x(3)-xc)*slope3
c	     if (slope1.lt.2.0*slope3)then
c	      xposs2=0.5*(xposs2+xoss1)
c	     endif
c         else
c         xposs2 = p(2)+xc*slope1
c		 if (slope1.lt.2.0*slope3)then
c	      xposs2=0.5*(xposs2+xoss1)
c	     endif
c         endif

c        xposs2=xposs1

c     case 2

      elseif (slope1.gt.slope2.and.slope2.gt.slope3) then  
	
         flag=2 

c	rp
	steps = 0
c   55   if (abs(pneed1-pneed2).ge.toler) then
   55   if (abs(pneed1-pneed2).ge.toler.and.steps.lt.100) then

c         rp
          steps = steps + 1          
          pneed1= p(2)+xneed*slope1
          pneed2= p(3)-(x(3)-xneed)*slope3
           if (pneed1.lt.pneed2) then 
	      xneed1=xneed1/2.0
            xneed=xneed+xneed1
           else 
	      xneed1=xneed1/2.0
            xneed=xneed-xneed1
           endif
            goto 55
         endif
c         rp
c	print*,steps
	   if ( xc.ge.xneed ) then
         xposs2 = p(3)-(x(3)-xc)*slope3
         else
         xposs2 = p(2)+xc*slope1
         endif


c     case 3 and 4 is a warning, can be deleted and use linear interpolation

c      elseif (slope1.eq.slope2.or.slope2.eq.slope3
c     $.or.slope2.eq.0.0) then 
	
c        xposs2=xposs1
      
      else
c      print*, 'adjustment maybe need'
       xposs2=xposs1
      endif

c     sort these possible of two values

         if (xposs1.lt.xposs2) then
	   temp= xposs1
	   xposs1=xposs2
	   xposs2=temp
         endif

      

	endif


      irow = 0

c

      do 100 i=1,4

c

c

      irow = irow + 1

      ab(irow,1) = 1.0

      ab(irow,2) = x(i)

      ab(irow,3) = x(i)*x(i)

      ab(irow,4) = x(i)*x(i)*x(i)

      ab(irow,5) = p(i)



  100 continue



c     call matrix solver



      call elim(ab,4,5,4)



      pcomp =

     $   ab(1,5) * 1.0D0 +

     $   ab(2,5) * xc +

     $   ab(3,5) * xc*xc +

     $   ab(4,5) * xc*xc*xc 

c     determine the value
      
c	if (flag.eq.2) then
	if (limiter.eq.1) then
         if (pcomp.gt.xposs1) then
	    pcomp=xposs1
        elseif (pcomp.lt.xposs2) then
	    pcomp=xposs2	
         endif     
      endif

      return

      end


      subroutine diff5(c,exch,dtbydx,area,theta,dt,
     $                	decay,dzex,dzarea,nodes,cd,cdp,cp)

c     Remember c is mass per unit length
c     but cp is concentration

      integer rank, nodes, l  

      real*8 c(nodes), cp(nodes), cd(nodes), cdp(nodes)
      real*8 a(nodes),b(nodes),c1(nodes),d(nodes),w(nodes)
      real*8 e(nodes), f(nodes), dt(nodes)
      real*8 exch(nodes+1)
      real*8 dtbydx(nodes) 
	real*8 area(nodes), decay(nodes)
	real*8 dzex(nodes), dzarea(nodes)
      real*8 x,y,z, theta
      real*8 rk1, rk2, alpha, beta

c     rk1 is alpha
c     rk2 is areas


c     -------------
c     SET UP MATRIX
c     -------------


c     Test for the type of link and insert
c     appropriate equation coefficients
c     ------------------------------------

            rank = 1
               l = 1

            rk1 = dzex(l)
            rk2 = dzex(l) * ( area(l) / dzarea(l) )

            alpha = theta*dt(l)*rk2 / (1.0D0 + theta*dt(l)*rk2)
            beta = (cd(l)+(1.0D0-theta)*dt(l)*rk2*(c(l)/area(l)-cd(l))) 
     $                     / (1.0D0 + theta*dt(l)*rk2)
         
            

c	      x = dzarea(l) + dt(l) * dzex(l) * area(l)
c     $	          + dt(l) * 0.0D0 * dzarea(l) 

c	      y = cd(l) * dzarea(l) / x

c	      z = dt(l) * dzex(l) * area(l) / x


c...........node -1

            a(rank) = 0.0

c...........node 1

            b(rank) = area(l) * 
     $       ( 1.0D0 - dt(l) * theta * dzex(l) * (alpha -1.0D0)
     $		   + theta * decay(l) * dt(l) )
     $	  + theta * dtbydx(l) * exch(l+1)

c...........node +1

            c1(rank) = -1.0D0 * theta * dtbydx(l) * exch(l+1)

            d(rank) = c(l) + beta * dzex(l) * area(l) * theta * dt(l)




         do l=2, nodes-1

            rank = rank + 1

c	      x = dzarea(l) + dt(l) * dzex(l) * area(l)
c     $	          + dt(l) * 0.0D0 * dzarea(l) 

c	      y = cd(l) * dzarea(l) / x

c	      z = dt(l) * dzex(l) * area(l) / x

            rk1 = dzex(l)
            rk2 = dzex(l) * ( area(l) / dzarea(l) )

            alpha = theta*dt(l)*rk2 / (1.0D0 + theta*dt(l)*rk2)
            beta = (cd(l)+(1.0D0-theta)*dt(l)*rk2*(c(l)/area(l)-cd(l))) 
     $                     / (1.0D0 + theta*dt(l)*rk2)


c...........node -1

            a(rank) = -1.0D0 * theta * dtbydx(l) * exch(l)

c...........node 1

            b(rank) = area(l) * 
     $       ( 1.0D0 - dt(l) * theta * dzex(l) * (alpha -1.0D0)
     $		   + theta * decay(l) * dt(l) )
     $                		  + theta * dtbydx(l) * exch(l) 
     $                          + theta * dtbydx(l) * exch(l+1)

c...........node +1

            c1(rank) = -1.0D0 * theta * dtbydx(l) * exch(l+1)

c

            d(rank) = c(l) + beta * dzex(l) * area(l) * theta * dt(l)
		   
		end do    


            l = nodes

            rank = rank + 1

c	      x = dzarea(l) + dt(l) * dzex(l) * area(l)
c     $	          + dt(l) * 0.0D0 * dzarea(l) 

c	      y = cd(l) * dzarea(l) / x

c	      z = dt(l) * dzex(l) * area(l) / x

            rk1 = dzex(l)
            rk2 = dzex(l) * ( area(l) / dzarea(l) )

            alpha = theta*dt(l)*rk2 / (1.0D0 + theta*dt(l)*rk2)
            beta = (cd(l)+(1.0D0-theta)*dt(l)*rk2*(c(l)/area(l)-cd(l))) 
     $                     / (1.0D0 + theta*dt(l)*rk2)


c...........node -1

            a(rank) = -1.0D0 * theta * dtbydx(l) * exch(l)

c...........node 1

            b(rank) = area(l) *    
     $       ( 1.0D0 - dt(l) * theta * dzex(l) * (alpha -1.0D0)
     $		   + theta * decay(l) * dt(l) ) 
     $	                     + theta * dtbydx(l) * exch(l) 

c...........node +1

            c1(rank) = 0.0

c

            d(rank) = c(l) + beta * dzex(l) * area(l) * theta * dt(l)





c     ----------------------

c     MATRIX SOLVING SECTION

c     ----------------------



c     Preparation of variables for matrix routines.

c     ---------------------------------------------



         nrank = rank



c     Calls to matrix routines

c     ------------------------

c           L-U Decompostion

c           ----------------



*     call gtri(a,b,c1,d,rank,e,f,w,1,0.0,0.0,1,c(nodes),0.0)


      call trid(a,b,c1,d,rank)



c     ---------------------

c     END OF MATRIX SECTION

c     ---------------------

c       Assign values from solution vector to 

c       cp values.

c       -------------------------------------



           do i = 1,nodes

            cp(i) = d(i)

c            x = dzarea(i) + dt(i) * dzex(i) * area(i)
c     $	          + dt(i) * 0.0D0 * dzarea(i) 

c	      y = cd(i) * dzarea(i) / x

c	      z = dt(i) * dzex(i) * area(i) / x


            rk1 = dzex(i)
            rk2 = dzex(i) * ( area(i) / dzarea(i) )

            alpha = theta*dt(i)*rk2 / (1.0D0 + theta*dt(i)*rk2)
            beta = (cd(i)+(1.0D0-theta)*dt(i)*rk2*(c(i)/area(i)-cd(i))) 
     $                     / (1.0D0 + theta*dt(i)*rk2)

            cdp(i) = beta + alpha * cp(i)

           end do

      end

