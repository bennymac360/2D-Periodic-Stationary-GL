      implicit real*8 (a-h,o-z)
      dimension
     1  x(256),y(256),ft1(256),ft2(256),
     2  qf(256,256),qqf(256,256),pp(256,256),hz(256,256),ppp(256,256),
     3  vx(256,256),vy(256,256),curr(256,256),vvx(256,256),vvy(256,256),
     4  sx(256,256),sy (256,256),ex(256,256),ey(256,256),currr(256,256),
     5  free_t(0:500),rnd(256,256),rnd2(256,256),sxx(256,256),
     6  aac(100000),bbc(100000),ccc(100000),ddc(100000),eec(100000),
     7  fc(100000),f(100000),cen(100000),cc(100000),amag(0:500),
     8  dpx(256),dpy(256),iss(256,256),ww(256,256),syy(256,256),
     9  dummy(256,256),dum2(256,256),dum3(256,256),php(256,256),
     :  xvp(100),yvp(100)

      complex*16 qf,aac,bbc,ccc,ddc,f,fc,aim,eec,bim,phase,qfn,qqf
      real*8 rndm,free_new,free_old
      character*7 text1,text2
      character*10 text3
              
********RANDOM GENERATOR********
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
      integer*2 seed(3),iii
      COMMON /aleat/seed
      seed(1) = 11
      seed(2) = 54
      seed(3) = -27
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
********************************

      idd=0			        !!! job counter
      iddm=200	            !!! maximum number of jobs 
      nloopfirst=8
      nloopafter=8 

      field=0.00d0		    !!! Applied field in units of Hc2(T)
      dfield=0.005d0		!!! Step for increasing magnetic field
      
      size_x=32.d0		    !!! Size of the simulation region in x-direction in units of xi(T)
      size_y=32.d0		    !!! Size of the simulation region in y-direction in units of xi(T)
       
      wstr=20.0d0			!!! Width of the stipe in units of xi(T)
      thick=1.5d0 		    !!! Thickness of the film in units of xi(T)
      ak1=4.0d0			    !!!	GL parameter kappa
	  
      
      iter=250000
      nit=500
      icc=6532
      free_new=1000.d0
      er=1.d-9
      err=1.d-1
      error=3.d-5           !!! Convergence error for first job
      error2=3.d-5          !!! Convergence error for subsequent jobs
      demp=2.0d0
      dempv=0.05d0
      nn=128				!!! number of grid points, has to be a power of 2 for fft
      iself=0				!!! to include the second equation (1) or not (0)
      iread=0				!!! to use the input (continue) file (iread=1) or not (iread=0); possibility to start from a state with Nv vortices, if iread=2
      icont=1               !!! write out continue files at every icont job
      defected=0            !!! Creates a defected stripe, 1=defected, 0=perfect
      
      
      
      pi=4*atan(1.d0)
	  
      hmin_x=size_x/nn
      hmin_y=size_y/nn
        
      akx=1/hmin_x**2
      aky=1/hmin_y**2
      
      do n=1,nn           !!! Here the code froms from 1-nn, so we must subtract the half length.
      x(n)=hmin_x*n-0.5d0*size_x
      y(n)=hmin_y*n-0.5d0*size_y
      end do
	
      addx=0.00d0           !!! This adds an applied electric field to the vector potential, acting similar to a current but without time evol.
      daddx=-0.005d0
      addy=0.d0
       
       
!!!!!! Now we create the stripe !!!!!!
**************************************

      do n=1,nn        !This creates defects in the stripe.
      do k=1,nn
      if(abs(y(k)).lt.0.5d0*wstr)then
	
      iss(k,n)=1
      
      ww(k,n)=thick
      
      
      if(defected.eq.1)then
      call random_number(dummy(k,n))
      call random_number(dum2(k,n))
      call random_number(dum3(k,n))
      ww(k,n)=((1.d0 -0.9d0)*dummy(k,n) + 0.9d0) *thick	!why randomize the thickness so much? up to 1.8 times the predefined thickness 'thick'? Makes little sense to me.
       
      if(iss(k,n).eq. 1 . and . iss(k+1,n).eq.0)ww(k,n)=((1.d0 -0.9d0)*
     :dum2(k,n) + 0.9d0)*ww(k,n)
      if(iss(k,n).eq. 1 . and . iss(k+1,n).eq.0)ww(k-1,n)=((1.d0-0.9d0)*
     :dum2(k,n)+0.9d0)*ww(k-1,n)
      if(iss(k,n).eq. 1 . and . iss(k+1,n).eq.0)ww(k-2,n)=((1.d0-0.9d0)*
     :dum2(k,n)+0.9d0)*ww(k-2,n)						!if k=1 then k-2=-1 so array bounds of ww are exceeded (defined as ww(256,256))
       
      if(iss(k,n). eq. 1 . and . iss(k-1,n).eq.0)ww(k,n)=((1.d0 -0.9d0)*
     :dum3(k,n) + 0.9d0)*ww(k,n)
      if(iss(k,n).eq. 1 . and . iss(k-1,n).eq.0)ww(k+1,n)=((1.d0-0.9d0)*
     :dum3(k,n)+ 0.9d0)*ww(k+1,n)
      if(iss(k,n).eq. 1 . and . iss(k-1,n).eq.0)ww(k+2,n)=((1.d0-0.9d0)*
     :dum3(k,n)+0.9d0)*ww(k+2,n)
      end if
      
      end if	! I have put everything under the condition that we are inside the stripe; that solves the k-2 issue as well.
             
      end do
      end do       

!!!!!! Create some Antidots below !!!!!!
****************************************

C       do k=50,nn-50,20            
C        do n=10,nn-10,20
C         do i=-3,4,1
C          do j=-3,4,1
C           iss(k+i,n+j)= 0.
C           end do
C         end do
C        end do
c      end do

	
!      if(iread.eq.1 .and. idd.eq.0)then
!      open(1,file='input.dat')
!      do k=1,nn
!      do n=1,nn
!      read(1,*)qf(k,n)
!      read(1,*)ex(k,n)
!      read(1,*)ey(k,n)
!      end do
!      end do
!      close(1)
!      end if !!of iread=1

!      if(iread.eq.2)then !***********************
      
!      i_vort=1	!how many vortices you want in the sample to begin with

!      do nv=1,i_vort				!giving the initial position of vortices, in the center of the stripe, random in x
!      rndm=RAND48(SEED)
!      xvp(nv)=0.5d0*size_x*(1.d0-2.d0*rndm)
!      yvp(nv)=0.d0
!      end do

!      do k=1,nn
!      do n=1,nn
      
!      ampl=1.d0
!      sum=0.d0
!      do nv=1,i_vort
!      bx=x(n)-xvp(nv)           ! REAL VORTEX
!      by=y(k)-yvp(nv)
!      a=sqrt(bx**2+by**2)+1.d-10
!      b=dacos(bx/a)
!      if(by.le.0.d0)b=2*pi-b
!      sum=sum+b
!      if(a.le.0.5d0)ampl=ampl*a/(1+a)
!      end do
!      if(iss(k,n).eq.1)qf(k,n)=ampl*dcmplx(cos(sum),sin(sum))
!      pp(k,n)=dreal(qf(k,n)*conjg(qf(k,n)))
!      end do
!      end do

!      end if !of iread=2

       
c     test what has been prepared as input (structiral, pinning, whatever) 
      call iso_write('test_input.dat',nn,ww)
!      stop



!!!! WE just finished setting up the geometry of the stripe !!!!!!
*******************************************************************


 
!!!!!! The First Job !!!!!!
***************************
C The first job can be used to try to find the lowest energy minimum by
C looping the first job and rewriting data that corresponds to a lower
C energy state. This is useful if the wavefunction for sequential jobs
C should start from the previous state. Useful for critical current 
C calculations when we increase the constant vector potential during
C each job. 
C
C Note: We do not need to loop so many times if we are trying to find 
C the lowest energy for different fields for example. The first job
C and subsequent jobs should be looped the same number of times.


!      if(iread.eq.0 .and. idd.eq.0)then !***********************
      
      if(idd.eq.0)then !***********************	THIS SHOULD BE EXECUTED IN THE FIRST ITERATION REGARDLESS OF IREAD!
      
      print *, "Magnetic field is", field
      print *, "addx field is", addx

      text2='000.dat'
      if(idd.le.9 )write(text2(3:3),'(i1)')idd
      if(idd.gt.9.and.idd.le.99)write(text2(2:3),'(i2)')idd
      if(idd.gt.99)write(text2(1:3),'(i3)')idd
      
!!!!!! LOOP HERE TO FIND THE LOWEST MINIMUM ENERGY !!!!!!
*********************************************************
      
      do iii=1,nloopfirst    ! CHANGE FOR VALUE FOR NUMBER IF LOOPS
      
      do k=1,nn
      do n=1,nn
      call random_number(rnd(k,n))
      ex(k,n)=-0.5d0*field*y(k)+addx
      ey(k,n)=0.5d0*field*x(n)+addy		!!!defining the initial vector potential
      end do
      end do
      
      ampl=1.d0
      sum=0.d0
      if(iss(k,n).eq.1)qf(k,n)=ampl*dcmplx(.1d0*rnd(k,n),0.d0)  !I COMMENTED THIS OUT, OTHERWISE YOU OVERWRITE ANY INITIAL STATE PREPARED BY IREAD 
      
      if(iread.eq.1)then
      open(1,file='cont002.dat')
      do k=1,nn
      do n=1,nn
      read(1,*)qf(k,n)
c      read(1,*)ex(k,n)
c      read(1,*)ey(k,n)
      end do
      end do
      close(1)
      end if !!of iread=1  

      if(iread.eq.2)then !***********************
      
      i_vort=1	!how many vortices you want in the sample to begin with

      do nv=1,i_vort				!giving the initial position of vortices, in the center of the stripe, random in x
      rndm=RAND48(SEED)
      xvp(nv)=0.5d0*size_x*(1.d0-2.d0*rndm)
      yvp(nv)=0.d0
      end do

      do k=1,nn
      do n=1,nn
      
      ampl=1.d0
      sum=0.d0
      do nv=1,i_vort
      bx=x(n)-xvp(nv)           ! REAL VORTEX
      by=y(k)-yvp(nv)
      a=sqrt(bx**2+by**2)+1.d-10
      b=dacos(bx/a)
      if(by.le.0.d0)b=2*pi-b
      sum=sum+b
      if(a.le.0.5d0)ampl=ampl*a/(1+a)
      end do
      if(iss(k,n).eq.1)qf(k,n)=ampl*dcmplx(cos(sum),sin(sum))
      pp(k,n)=dreal(qf(k,n)*conjg(qf(k,n)))
      end do
      end do

      end if !of iread=2
      
      

      do it=1,iter   !----------iterations---------


      text3='000000.dat'

      if(it.le.9)write(text3(6:6),'(i1)')it
      if(it.gt.9.and.it.le.99)write(text3(5:6),'(i2)')it
      if(it.gt.99.and.it.le.999)write(text3(4:6),'(i3)')it
      if(it.gt.999.and.it.le.9999)write(text3(3:6),'(i4)')it
      if(it.gt.9999.and.it.le.99999)write(text3(2:6),'(i5)')it
      if(it.gt.99999.and.it.le.999999)write(text3(1:6),'(i6)')it

      do n=1,nn
      dpx(n)=0.5d0*field*size_x*y(n)	!gauge at the boundary
      dpy(n)=-0.5d0*field*size_y*x(n)
      end do

      delt=0.d0

      do k=1,nn
      do n=1,nn
      l=nn*(k-1)+n
      f(l)=qf(k,n)
      aac(l)=dcmplx(0.d0,0.d0)
      bbc(l)=dcmplx(0.d0,0.d0)
      ccc(l)=dcmplx(0.d0,0.d0)
      ddc(l)=dcmplx(0.d0,0.d0)
      fc(l)=dcmplx(0.d0,0.d0)
      if(iss(k,n).eq.1)then
      cen(l)=-1.d0

      aim=dcmplx(0.d0,0.d0)

      if(k.ne. 1)then ! AAAA
      if(iss(k-1,n).eq.1)then
      aac(l)=aky*exp(dcmplx(0.d0, hmin_y*ey(k,n)))
      aim=aim-aac(l)*qf(k-1,n)
      cen(l)=cen(l)+aky
      end if
      else
      if(iss(nn,n).eq.1)then
      aac(l)=aky*exp(dcmplx(0.d0, hmin_y*ey(k,n)-dpy(n)))
      aim=aim-aac(l)*qf(nn,n)
      cen(l)=cen(l)+aky
      end if
      end if

      if(k.ne.nn)then ! CCCC
      if(iss(k+1,n).eq.1)then
      ccc(l)=aky*exp(dcmplx(0.d0,-hmin_y*ey(k+1,n)))
     :*ww(k+1,n)/ww(k,n)
      aim=aim-ccc(l)*qf(k+1,n)
      cen(l)=cen(l)+aky*ww(k+1,n)/ww(k,n)
      end if
      else
      if(iss(1,n).eq.1)then
      ccc(l)=aky*exp(dcmplx(0.d0,-hmin_y*ey(1,n)+dpy(n)))
     :*ww(1,n)/ww(k,n)
      aim=aim-ccc(l)*qf(1,n)
      cen(l)=cen(l)+aky*ww(1,n)/ww(k,n)
      end if
      end if

      if(n.ne.1)then   ! BBBB
      if(iss(k,n-1).eq.1)then
      bbc(l)=akx*exp(dcmplx(0.d0, hmin_x*ex(k,n)))
      aim=aim-bbc(l)*qf(k,n-1)
      cen(l)=cen(l)+akx
      end if
      else
      if(iss(k,nn).eq.1)then
      bbc(l)=akx*exp(dcmplx(0.d0, hmin_x*ex(k,n)-dpx(k)))
      aim=aim-bbc(l)*qf(k,nn)
      cen(l)=cen(l)+akx
      end if
      end if

      if(n.ne.nn)then   ! DDDD
      if(iss(k,n+1).eq.1)then
      ddc(l)=akx*exp(dcmplx(0.d0,-hmin_x*ex(k,n+1)))
     :*ww(k,n+1)/ww(k,n)
      aim=aim-ddc(l)*qf(k,n+1)
      cen(l)=cen(l)+akx*ww(k,n+1)/ww(k,n)
      end if
      else
      if(iss(k,1).eq.1)then
      ddc(l)=akx*exp(dcmplx(0.d0,-hmin_x*ex(k,1)+dpx(k)))
     :*ww(k,1)/ww(k,n)
      aim=aim-ddc(l)*qf(k,1)
      cen(l)=cen(l)+akx*ww(k,1)/ww(k,n)
      end if
      end if

      a=dreal(f(l)*conjg(f(l)))
      aim=(cen(l)+a)*qf(k,n)+aim
      delt=dmax1(delt,dreal(aim*conjg(aim)))

      eec(l)=f(l)*f(l)
      fc(l)=demp*f(l)+2*a*f(l)
      cen(l)=cen(l)+demp+2*a
      else
      cen(l)=1.d0
      end if
      end do
      end do

      delt=sqrt(dabs(delt))
      nit1=nit
      call relax_o(cen,eec,aac,bbc,ccc,ddc,fc,f,nn,er,err,nit1,as,bs)
	
!%%%%%%%%%%*THE WAVE FUNCTION*%%%%%%%%%%%
      delta=0.d0

      do k=1,nn
      do n=1,nn
      l=nn*(k-1)+n
      aim=qf(k,n)-f(l)
      qf(k,n)=f(l)
      delta=dmax1(delta,dreal(aim*conjg(aim)))
      pp(k,n)=dreal(qf(k,n)*conjg(qf(k,n)))
      end do
      end do
      delta=sqrt(delta)

!%%%%%%%%%%%*CURRENT DENSITIES*%%%%%%%%%%
            
      do k=1,nn
      do n=1,nn
      vx(k,n)=0.d0
      vy(k,n)=0.d0

      if(k.ne.nn)then
      aim=exp(-hmin_y*dcmplx(0.d0,ey(k,n)))
      vy(k,n)=dreal((conjg(qf(k,n))*(aim*qf(k+1,n)-qf(k,n))
     :-conjg(aim*qf(k+1,n)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_y))
      else
      aim=exp(dcmplx(0.d0,-hmin_y*ey(1,n)+dpy(n)))
      vy(k,n)=dreal((conjg(qf(k,n))*(aim*qf(1,n)-qf(k,n))
     :-conjg(aim*qf(1,n)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_y))
      end if

          
      if(n.ne.nn)then
      aim=exp(-hmin_x*dcmplx(0.d0,ex(k,n)))
      vx(k,n)=dreal((conjg(qf(k,n))*(aim*qf(k,n+1)-qf(k,n))
     :-conjg(aim*qf(k,n+1)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_x))
      else
      aim=exp(dcmplx(0.d0,-hmin_x*ex(k,n)+dpx(k)))
      vx(k,n)=dreal((conjg(qf(k,n))*(aim*qf(k,1)-qf(k,n))
     :-conjg(aim*qf(k,1)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_x))
      end if

      sx(k,n)=vx(k,n)
      sy(k,n)=vy(k,n)
        
      end do
      end do
   
           
      do k=1,nn
      do n=1,nn
      curr(k,n)=sqrt(vx(k,n)*vx(k,n)+vy(k,n)*vy(k,n))
      end do
      end do
	
      if(iself.eq.1)then ! 2ND EQUATION

!%%%%%%%%%%%%*DIRECT FOURIER TRANSFORM*%%%%%%%%%%%%%%
      do k=1,nn
      do n=1,nn
      ft1(n)=sx(k,n)
      ft2(n)=sy(k,n)
      end do
      call realft(ft1,nn/2,1)
      call realft(ft2,nn/2,1)
      do n=1,nn
      sx(k,n)=2*ft1(n)/nn
      sy(k,n)=2*ft2(n)/nn
      end do
      end do
      do n=1,nn
      do k=1,nn
      ft1(k)=sx(k,n)
      ft2(k)=sy(k,n)
      end do
      call realft(ft1,nn/2,1)
      call realft(ft2,nn/2,1)
      do k=1,nn
      sx(k,n)=2*ft1(k)/nn
      sy(k,n)=2*ft2(k)/nn
      end do
      end do
!%%%%%%%%%%%%%%*BACK FOURIER TRANSFORM*%%%%%%%%%%%%%%%%%%
      
      do k=1,nn
      do n=1,nn
      if(k.eq.2)then
      qy=2*pi*int(nn/2   )/size_y
      else
      qy=2*pi*int((k-1)/2)/size_y
      end if
      if(n.eq.2)then
      qx=2*pi*int(nn/2   )/size_x
      else
      qx=2*pi*int((n-1)/2)/size_x
      end if
      q=sqrt(qx**2+qy**2)
      
      if(q.gt.1.d-6)then
      sx(k,n)=sx(k,n)*(1.d0-exp(-ww(k,n)*q/2.d0))/(q**2)/ak1**2
      sy(k,n)=sy(k,n)*(1.d0-exp(-ww(k,n)*q/2.d0))/(q**2)/ak1**2
      else
      sx(k,n)=0.d0
      sy(k,n)=0.d0
      end if
      end do
      end do

      do k=1,nn
      do n=1,nn
      ft1(n)=sx(k,n)
      ft2(n)=sy(k,n)
      end do
      call realft(ft1,nn/2,-1)
      call realft(ft2,nn/2,-1)
      do n=1,nn
      sx(k,n)=ft1(n)
      sy(k,n)=ft2(n)
      end do
      end do

      do n=1,nn
      do k=1,nn
      ft1(k)=sx(k,n)
      ft2(k)=sy(k,n)
      end do
      call realft(ft1,nn/2,-1)
      call realft(ft2,nn/2,-1)
      do k=1,nn
      sx(k,n)=ft1(k)
      sy(k,n)=ft2(k)
      end do
      end do

      deltae=0.d0
      do k=1,nn
      do n=1,nn
      a=sx(k,n)-0.5d0*field*y(k)+addx
      deltae=dmax1(deltae,dabs(a-ex(k,n)))
      ex(k,n)=(1-dempv)*ex(k,n)+dempv*a
      
      a=sy(k,n)+0.5d0*field*x(n)+addy
      deltae=dmax1(deltae,dabs(a-ey(k,n)))
      ey(k,n)=(1-dempv)*ey(k,n)+dempv*a
      end do
      end do
      else
      do k=1,nn
      do n=1,nn
      sx(k,n)=0.d0
      sy(k,n)=0.d0
      ex(k,n)=-0.5d0*field*y(k)+addx
      ey(k,n)=0.5d0*field*x(n)+addy
      end do
      end do
      end if 

      dens_f=0.d0
      dens_h=0.d0
      volume=0.d0
      do k=1,nn
      do n=1,nn
      if(iss(k,n).eq.1)volume=volume+ww(k,n)*hmin_x*hmin_y !! counts the volume of superconducting material
      dens_f=dens_f+ww(k,n)*hmin_x*hmin_y*pp(k,n)**2 !!counts free energy from sc
      dens_h=dens_h+ww(k,n)*hmin_x*hmin_y*
     :(sx(k,n)*vx(k,n)+sy(k,n)*vy(k,n))  !!counts energy from magnetic field
      end do
      end do
 
      free=(-dens_f/2+dens_h)/volume    
      delta_f=dabs(free-free_0)
      free_0=free
	
      if(mod(it,100).eq.0)then
      print*,it,dmax1(delt,delta),free
      end if
		
!       if(mod(it,1000).eq.0)then
 !      call iso_write('dens_tussen_'//text3,nn,pp)
 !      end if
       
       

      if(dmax1(delt,delta).le.error.and.it.gt.10) go to 1
    

      end do !iteration

1     continue
      if(iii.eq.1)then
        free_old=free
        vvx=vx
        vvy=vy
        ppp = pp
        currr=curr
        sxx=sx
        syy=sy
        qqf =qf
      end if
      if(iii.gt.1)then
        if(free.lt.free_old)then
            free_old=free
            vvx=vx
            vvy=vy
            ppp = pp
            currr=curr
            sxx=sx
            syy=sy
            qqf=qf
      end  if
      end  if
      
      free_t(idd)=free_old
	free_new=free
      
      
      call iso_write('dens_'//text2,nn,ppp)
		
      do k=1,nn
      do n=1,nn
      if(iss(k,n).eq.1)then
      aim=qqf(k,n)
      a=dreal(aim)
      b=dimag(aim)
      c=sqrt(a**2+b**2+1.d-10)
      c=dacos(a/c)
      if(b.le.0.d0)c=2*pi-c
      php(k,n)=c/(2*pi)
      end if
      end do
      end do
      
      call iso_write1('phase_'//text2,nn,php)
      call iso_write('curramp_'//text2,nn,currr)
C      call iso_write('curxcom_'//text2,nn,vvx)
C      call iso_write('curycom_'//text2,nn,vvy)
      
      
      if(mod(idd,icont).eq.0) then
      open(1,file='cont'//text2)
      do k=1,nn
      do n=1,nn
      write(1,*)qqf(k,n)
      end do
      end do
      close(1) 
      end if


      do k=1,nn
      do n=2,nn
      hz(k,n)=(vvy(k,n)-vvy(k,n-1))/hmin_x ! MAGNETIC FIELD
      end do
      hz(k,1)=(vvy(k,1)-vvy(k, nn))/hmin_x ! MAGNETIC FIELD
      end do
      do n=1,nn
      do k=2,nn
      hz(k,n)=hz(k,n)-(vvx(k,n)-vvx(k-1,n))/hmin_y ! MAGNETIC FIELD
      end do
      hz(1,n)=hz(1,n)-(vvx(1,n)-vvx(nn ,n))/hmin_y ! MAGNETIC FIELD
      end do
      
      
c	call iso_write1('magn_'//text2,nn,hz)
c	call iso_write2('magn_'//text2,nn,x,y,hz)
c	call vector('curre'//text2,ng,x,vx,vy)
        
        
       a=0.d0
       vol_n=0.d0
       do k=0,nn
       do n=0,nn
       vol_n=vol_n+ww(k,n)*hmin_x*hmin_y
       a=a+hz(k,n)*ww(k,n)*hmin_x*hmin_y
       end do
       end do
       amag(idd)=-a/(4*pi*vol_n)
       
       
       open(1,file='free_en.dat')
       open(2,file='magnet.dat')
       do n=0,idd
       write(1,*)n,free_t(n)
       write(2,*)n,amag(n)
       end do
       close(1)
       close(2)
              
      end do
      end if

!!!!!! END OF FIRST JOB !!!!!!
******************************

***********************************************************************

***************************************
!!!!!! NOW GO TO SUBSEQUENT JOBS !!!!!!
*************************************** 
333   continue
!!!!!! INCREASE JOB STEP !!!!!!
*******************************
      idd=idd+1      
!!!!!! INCREASE PARAMETER FOR JOB !!!!!!
****************************************      
      if(idd.lt.iddm)then
C      addx = addx + daddx            !!!This increases the field from its minimum up to a value depeneding on the field step and iddm(= no. of jobs)
      field=field+dfield



      print *, "Magnetic field is", field
      print *, "addx field is", addx

      text2='000.dat'
      if(idd.le.9 )write(text2(3:3),'(i1)')idd
      if(idd.gt.9.and.idd.le.99)write(text2(2:3),'(i2)')idd
      if(idd.gt.99)write(text2(1:3),'(i3)')idd
      
      
!!!!!! STARTS SUBSEQUENT JOBS !!!!!!      
      if(idd.gt.0)then ! THIS CHECK TO SEE IF WE ARE AT JOBS AFTER THE FIRST
      
!!!!!! LOOP FOR SUBSEQUENT JOBS !!!!!!
**************************************
C If the simulations uses previous states then do not use random 
C wavefunctions (qf), and loop only ONCE!
C
C For jobs locating free energy minimum at every job then loop a number
C of times and use random wave function.
      
      do iii=1,nloopafter   !!! Loop to find lower free energy, ends after the "1  continue" around line 430
        
      do k=1,nn
      do n=1,nn
      call random_number(rnd2(k,n))
      ex(k,n)=-0.5d0*field*y(k)+addx
      ey(k,n)=0.5d0*field*x(n)+addy
      
c      qf(k,n)=qqf(k,n)                                          !!! THIS USES PREVIOUS MINIMIZED WAVEFUNCTION.
      
      qf(k,n) = dcmplx(0.1d0*rnd2(k,n),0.0d0)        !!! Start the iteration with some random wavefunction.

!!!!!! START WITH FULLY SUPERCONDUCTING STATE !!!!!!
C      qf(k,n)=dcmplx(0.d0,0.d0)                                !!! Matrix size KxN of complex type full of zeros.
C      if(iss(k,n).eq.1)qf(k,n)=dcmplx(1.d0,0.d0)               !!! if in the superconducting state we apply a "vector" in the k,n (x,y) direction

      enddo
      enddo
      
      if(iread.eq.2)then !***********************
      
      i_vort=1	!how many vortices you want in the sample to begin with

      do nv=1,i_vort				!giving the initial position of vortices, in the center of the stripe, random in x
      rndm=RAND48(SEED)
      xvp(nv)=0.5d0*size_x*(1.d0-2.d0*rndm)
      yvp(nv)=0.d0
      end do

      do k=1,nn
      do n=1,nn
      
      ampl=1.d0
      sum=0.d0
      do nv=1,i_vort
      bx=x(n)-xvp(nv)           ! REAL VORTEX
      by=y(k)-yvp(nv)
      a=sqrt(bx**2+by**2)+1.d-10
      b=dacos(bx/a)
      if(by.le.0.d0)b=2*pi-b
      sum=sum+b
      if(a.le.0.5d0)ampl=ampl*a/(1+a)
      end do
      if(iss(k,n).eq.1)qf(k,n)=ampl*dcmplx(cos(sum),sin(sum))
      pp(k,n)=dreal(qf(k,n)*conjg(qf(k,n)))
      end do
      end do

      end if !of iread=2
      
      
      
    
      do it=1,iter   !----------iterations---------


      text3='000000.dat'

      if(it.le.9)write(text3(6:6),'(i1)')it
      if(it.gt.9.and.it.le.99)write(text3(5:6),'(i2)')it
      if(it.gt.99.and.it.le.999)write(text3(4:6),'(i3)')it
      if(it.gt.999.and.it.le.9999)write(text3(3:6),'(i4)')it
      if(it.gt.9999.and.it.le.99999)write(text3(2:6),'(i5)')it
      if(it.gt.99999.and.it.le.999999)write(text3(1:6),'(i6)')it

      do n=1,nn
      dpx(n)=0.5d0*field*size_x*y(n)	!gauge at the boundary
      dpy(n)=-0.5d0*field*size_y*x(n)
      end do

      delt=0.d0

      do k=1,nn
      do n=1,nn
      l=nn*(k-1)+n
      f(l)=qf(k,n)
      aac(l)=dcmplx(0.d0,0.d0)
      bbc(l)=dcmplx(0.d0,0.d0)
      ccc(l)=dcmplx(0.d0,0.d0)
      ddc(l)=dcmplx(0.d0,0.d0)
      fc(l)=dcmplx(0.d0,0.d0)
      if(iss(k,n).eq.1)then
      cen(l)=-1.d0

      aim=dcmplx(0.d0,0.d0)

      if(k.ne. 1)then ! AAAA
      if(iss(k-1,n).eq.1)then
      aac(l)=aky*exp(dcmplx(0.d0, hmin_y*ey(k,n)))
      aim=aim-aac(l)*qf(k-1,n)
      cen(l)=cen(l)+aky
      end if
      else
      if(iss(nn,n).eq.1)then
      aac(l)=aky*exp(dcmplx(0.d0, hmin_y*ey(k,n)-dpy(n)))
      aim=aim-aac(l)*qf(nn,n)
      cen(l)=cen(l)+aky
      end if
      end if

      if(k.ne.nn)then ! CCCC
      if(iss(k+1,n).eq.1)then
      ccc(l)=aky*exp(dcmplx(0.d0,-hmin_y*ey(k+1,n)))
     :*ww(k+1,n)/ww(k,n)
      aim=aim-ccc(l)*qf(k+1,n)
      cen(l)=cen(l)+aky*ww(k+1,n)/ww(k,n)
      end if
      else
      if(iss(1,n).eq.1)then
      ccc(l)=aky*exp(dcmplx(0.d0,-hmin_y*ey(1,n)+dpy(n)))
     :*ww(1,n)/ww(k,n)
      aim=aim-ccc(l)*qf(1,n)
      cen(l)=cen(l)+aky*ww(1,n)/ww(k,n)
      end if
      end if

      if(n.ne.1)then   ! BBBB
      if(iss(k,n-1).eq.1)then
      bbc(l)=akx*exp(dcmplx(0.d0, hmin_x*ex(k,n)))
      aim=aim-bbc(l)*qf(k,n-1)
      cen(l)=cen(l)+akx
      end if
      else
      if(iss(k,nn).eq.1)then
      bbc(l)=akx*exp(dcmplx(0.d0, hmin_x*ex(k,n)-dpx(k)))
      aim=aim-bbc(l)*qf(k,nn)
      cen(l)=cen(l)+akx
      end if
      end if

      if(n.ne.nn)then   ! DDDD
      if(iss(k,n+1).eq.1)then
      ddc(l)=akx*exp(dcmplx(0.d0,-hmin_x*ex(k,n+1)))
     :*ww(k,n+1)/ww(k,n)
      aim=aim-ddc(l)*qf(k,n+1)
      cen(l)=cen(l)+akx*ww(k,n+1)/ww(k,n)
      end if
      else
      if(iss(k,1).eq.1)then
      ddc(l)=akx*exp(dcmplx(0.d0,-hmin_x*ex(k,1)+dpx(k)))
     :*ww(k,1)/ww(k,n)
      aim=aim-ddc(l)*qf(k,1)
      cen(l)=cen(l)+akx*ww(k,1)/ww(k,n)
      end if
      end if

      a=dreal(f(l)*conjg(f(l)))
      aim=(cen(l)+a)*qf(k,n)+aim
      delt=dmax1(delt,dreal(aim*conjg(aim)))

      eec(l)=f(l)*f(l)
      fc(l)=demp*f(l)+2*a*f(l)
      cen(l)=cen(l)+demp+2*a
      else
      cen(l)=1.d0
      end if
      end do
      end do

      delt=sqrt(dabs(delt))
      nit1=nit
      call relax_o(cen,eec,aac,bbc,ccc,ddc,fc,f,nn,er,err,nit1,as,bs)
	
!%%%%%%%%%%*THE WAVE FUNCTION*%%%%%%%%%%%
      delta=0.d0

      do k=1,nn
      do n=1,nn
      l=nn*(k-1)+n
      aim=qf(k,n)-f(l)
      qf(k,n)=f(l)
      delta=dmax1(delta,dreal(aim*conjg(aim)))
      pp(k,n)=dreal(qf(k,n)*conjg(qf(k,n)))
      end do
      end do
      delta=sqrt(delta)

!%%%%%%%%%%%*CURRENT DENSITIES*%%%%%%%%%%
            
      do k=1,nn
      do n=1,nn
      vx(k,n)=0.d0
      vy(k,n)=0.d0

      if(k.ne.nn)then
      aim=exp(-hmin_y*dcmplx(0.d0,ey(k,n)))
      vy(k,n)=dreal((conjg(qf(k,n))*(aim*qf(k+1,n)-qf(k,n))
     :-conjg(aim*qf(k+1,n)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_y))
      else
      aim=exp(dcmplx(0.d0,-hmin_y*ey(1,n)+dpy(n)))
      vy(k,n)=dreal((conjg(qf(k,n))*(aim*qf(1,n)-qf(k,n))
     :-conjg(aim*qf(1,n)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_y))
      end if

          
      if(n.ne.nn)then
      aim=exp(-hmin_x*dcmplx(0.d0,ex(k,n)))
      vx(k,n)=dreal((conjg(qf(k,n))*(aim*qf(k,n+1)-qf(k,n))
     :-conjg(aim*qf(k,n+1)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_x))
      else
      aim=exp(dcmplx(0.d0,-hmin_x*ex(k,n)+dpx(k)))
      vx(k,n)=dreal((conjg(qf(k,n))*(aim*qf(k,1)-qf(k,n))
     :-conjg(aim*qf(k,1)-qf(k,n))*qf(k,n))
     :/(4.*dcmplx(0.d0,1.d0)*hmin_x))
      end if

      sx(k,n)=vx(k,n)
      sy(k,n)=vy(k,n)
        
      end do
      end do
   
           
      do k=1,nn
      do n=1,nn
      curr(k,n)=sqrt(vx(k,n)*vx(k,n)+vy(k,n)*vy(k,n))
      end do
      end do
	
      if(iself.eq.1)then ! 2ND EQUATION

!%%%%%%%%%%%%*DIRECT FOURIER TRANSFORM*%%%%%%%%%%%%%%
      do k=1,nn
      do n=1,nn
      ft1(n)=sx(k,n)
      ft2(n)=sy(k,n)
      end do
      call realft(ft1,nn/2,1)
      call realft(ft2,nn/2,1)
      do n=1,nn
      sx(k,n)=2*ft1(n)/nn
      sy(k,n)=2*ft2(n)/nn
      end do
      end do
      do n=1,nn
      do k=1,nn
      ft1(k)=sx(k,n)
      ft2(k)=sy(k,n)
      end do
      call realft(ft1,nn/2,1)
      call realft(ft2,nn/2,1)
      do k=1,nn
      sx(k,n)=2*ft1(k)/nn
      sy(k,n)=2*ft2(k)/nn
      end do
      end do
!%%%%%%%%%%%%%%*BACK FOURIER TRANSFORM*%%%%%%%%%%%%%%%%%%
      
      do k=1,nn
      do n=1,nn
      if(k.eq.2)then
      qy=2*pi*int(nn/2   )/size_y
      else
      qy=2*pi*int((k-1)/2)/size_y
      end if
      if(n.eq.2)then
      qx=2*pi*int(nn/2   )/size_x
      else
      qx=2*pi*int((n-1)/2)/size_x
      end if
      q=sqrt(qx**2+qy**2)
      
      if(q.gt.1.d-6)then
      sx(k,n)=sx(k,n)*(1.d0-exp(-ww(k,n)*q/2.d0))/(q**2)/ak1**2
      sy(k,n)=sy(k,n)*(1.d0-exp(-ww(k,n)*q/2.d0))/(q**2)/ak1**2
      else
      sx(k,n)=0.d0
      sy(k,n)=0.d0
      end if
      end do
      end do

      do k=1,nn
      do n=1,nn
      ft1(n)=sx(k,n)
      ft2(n)=sy(k,n)
      end do
      call realft(ft1,nn/2,-1)
      call realft(ft2,nn/2,-1)
      do n=1,nn
      sx(k,n)=ft1(n)
      sy(k,n)=ft2(n)
      end do
      end do

      do n=1,nn
      do k=1,nn
      ft1(k)=sx(k,n)
      ft2(k)=sy(k,n)
      end do
      call realft(ft1,nn/2,-1)
      call realft(ft2,nn/2,-1)
      do k=1,nn
      sx(k,n)=ft1(k)
      sy(k,n)=ft2(k)
      end do
      end do

      deltae=0.d0
      do k=1,nn
      do n=1,nn
      a=sx(k,n)-0.5d0*field*y(k)+addx
      deltae=dmax1(deltae,dabs(a-ex(k,n)))
      ex(k,n)=(1-dempv)*ex(k,n)+dempv*a
      
      a=sy(k,n)+0.5d0*field*x(n)+addy
      deltae=dmax1(deltae,dabs(a-ey(k,n)))
      ey(k,n)=(1-dempv)*ey(k,n)+dempv*a
      end do
      end do
      else
      do k=1,nn
      do n=1,nn
      sx(k,n)=0.d0
      sy(k,n)=0.d0
      ex(k,n)=-0.5d0*field*y(k)+addx
      ey(k,n)=0.5d0*field*x(n)+addy
      end do
      end do
      end if 

      dens_f=0.d0
      dens_h=0.d0
      volume=0.d0
      do k=1,nn
      do n=1,nn
      if(iss(k,n).eq.1)volume=volume+ww(k,n)*hmin_x*hmin_y !! counts the volume of superconducting material
      dens_f=dens_f+ww(k,n)*hmin_x*hmin_y*pp(k,n)**2 !!counts free energy from sc
      dens_h=dens_h+ww(k,n)*hmin_x*hmin_y*
     :(sx(k,n)*vx(k,n)+sy(k,n)*vy(k,n))  !!counts energy from magnetic field
      end do
      end do
 
      free=(-dens_f/2+dens_h)/volume    
      delta_f=dabs(free-free_0)
      free_0=free
	
      if(mod(it,100).eq.0)then
      print*,it,dmax1(delt,delta),free
      end if
		
!       if(mod(it,1000).eq.0)then
 !      call iso_write('dens_tussen_'//text3,nn,pp)
 !      end if
       
       

      if(dmax1(delt,delta).le.error2.and.it.gt.10) go to 2
    

      end do !iteration

2     continue
      if(iii.eq.1)then
        free_old=free
        vvx=vx			! I'm not sure what is the purpose of this. To save the lowest energy state in a different variable? If so, you may not mess with ppp later on...
        vvy=vy
        ppp = pp
        currr=curr
        sxx=sx
        syy=sy
        qqf=qf
      end if
      if(iii.gt.1)then
        if(free.lt.free_old)then
            free_old=free
            vvx=vx
            vvy=vy
            ppp = pp
            currr=curr
            sxx=sx
            syy=sy
            qqf=qf
        end if
      end  if
      
      free_t(idd)=free_old
      free_new=free
      
      
      call iso_write('dens_'//text2,nn,ppp)
		
      do k=1,nn
      do n=1,nn
      if(iss(k,n).eq.1)then
      aim=qqf(k,n)
      a=dreal(aim)
      b=dimag(aim)
      c=sqrt(a**2+b**2+1.d-10)
      c=dacos(a/c)
      if(b.le.0.d0)c=2*pi-c
c      ppp(k,n)=c/(2*pi)			! it is plain s...illy to use ppp here. Why not introduce a new variable for phase??? 
      php(k,n)=c/(2*pi)
      end if
      end do
      end do
      
      call iso_write1('phase_'//text2,nn,php)
      call iso_write('curramp_'//text2,nn,currr)
C      call iso_write('curxcom_'//text2,nn,vvx)
C      call iso_write('curycom_'//text2,nn,vvy)
      
      
      if(mod(idd,icont).eq.0) then
      open(1,file='cont'//text2)
      do k=1,nn
      do n=1,nn
      write(1,*)qqf(k,n)
      end do
      end do
      close(1) 
      end if


      do k=1,nn
      do n=2,nn
      hz(k,n)=(vvy(k,n)-vvy(k,n-1))/hmin_x ! MAGNETIC FIELD
      end do
      hz(k,1)=(vvy(k,1)-vvy(k, nn))/hmin_x ! MAGNETIC FIELD
      end do
      do n=1,nn
      do k=2,nn
      hz(k,n)=hz(k,n)-(vvx(k,n)-vvx(k-1,n))/hmin_y ! MAGNETIC FIELD
      end do
      hz(1,n)=hz(1,n)-(vvx(1,n)-vvx(nn ,n))/hmin_y ! MAGNETIC FIELD
      end do
      
      
!       call iso_write1('magn_'//text2,nn,hz)
c       call iso_write2('magn_'//text2,nn,x,y,hz)
        !!call vector('curre'//text2,ng,x,vx,vy)
        
        
       a=0.d0
       vol_n=0.d0
       do k=0,nn
       do n=0,nn
       vol_n=vol_n+ww(k,n)*hmin_x*hmin_y
       a=a+hz(k,n)*ww(k,n)*hmin_x*hmin_y
       end do
       end do
       amag(idd)=-a/(4*pi*vol_n)
       
       
       open(11,file='free_en.dat')
       open(22,file='magnet.dat')
       do n=0,idd
       write(11,*)n,free_t(n)
       write(22,*)n,amag(n)
       end do
      close(11)
      close(22)
              
      end do
      endif      
      
       
       continue

       go to 333
       endif                            !!!! THIS IS THE END OF THE ITERATION LOOP CHECK LOWEST FREE ENERGY FOR THE JON STEP HERE!!!!
       
     
       continue
       
!!!!!! ENDS THE PROGRAM SIMULATION HERE !!!!!       
       END 





**************************PHASE************************
      subroutine vector(name_f,nn,x,fx,fy)
      implicit real*8 (a-h,o-y)
      dimension x(-512:512),fx(-257:257,-257:257),
     :fy(-257:257,-257:257)
      character*(*) name_f
      d=0.d0
      do k=-nn,nn
      do n=-nn,nn
      d=dmax1(d,(fx(k,n)**2+fy(k,n)**2))
      end do
      end do
c       if(d.gt.1.d-6)then
      open(3,file=name_f)
      pi=4*datan(1.d0)
      do k=-nn,nn,4
      do n=-nn,nn,4
      a=fx(k,n)**2+fy(k,n)**2+1.d-10
c       if(a.gt.1.d-4*d)then
      a=sqrt(a)
      b=dacos(fx(k,n)/a)+1.d-10
      if(fy(k,n).lt.0.d0)b=2*pi-b
      write(3,'(4(1x,1pe12.5))')x(n),x(k),b,1.d1*a/sqrt(d)
c       end if
      end do
      end do
      close(3)
c       end if
      end


********************  Relaxation subroutine *****************************
      subroutine relax_o(cen,e,a,b,c,d,fc,f,nn,er,err,nit,as,bs)
      implicit complex*16 (a-h,o-z)
      real*8 er,err,as,bs,cen
      dimension
     :a(100000),b(100000),c(100000),d(100000),fc(100000),e(100000),
     :f(100000),cen(100000)
      do i=1,nit
      bs=0.d0

      do k=1,nn
      do n=1,nn
      l=nn*(k-1)+n
      aa=fc(l)

      if(k.ne.1)then
      aa=aa+a(l)*f(l-nn)
      else
      aa=aa+a(l)*f(nn*(nn-1)+n)
      end if

      if(k.ne.nn)then
      aa=aa+c(l)*f(l+nn)
      else
      aa=aa+c(l)*f(n)
      end if

      if(n.ne.1)then
      aa=aa+b(l)*f(l-1)
      else
      aa=aa+b(l)*f(nn*(k-1)+nn)
      end if

      if(n.ne.nn)then
      aa=aa+d(l)*f(l+1)
      else
      aa=aa+d(l)*f(nn*(k-1)+1)
      end if

      bb=aa-cen(l)*f(l)-e(l)*conjg(f(l))
      bs=dmax1(bs,dreal(bb*conjg(bb)))
      f(l)=(aa-e(l)*conjg(aa)/cen(l))
     :/(cen(l)-e(l)*conjg(e(l))/cen(l))
      end do
      end do
      if(i.eq.1)as=dmax1(bs,1.d-20)
      if(bs.le.er.or.(bs/as).le.err.and.i.ne.1) go to 1
      end do
  1   nit=i
      end



*********************************************************
*************************ISO_WRITE***********************
*********************************************************
      subroutine iso_write(name_f1,nn,f)	!I removed x and y from this subroutine (and all its calls), as it is not needed
      implicit real*8 (a-h,o-y)
      dimension f(256,256),f1(0:256,0:256)
      character*(*) name_f1

      do k=1,nn
      do n=1,nn
      f1(k,n)=f(k,n)
      end do
      end do
      f1(0,0)=f(nn,nn)
      do n=1,nn
      f1(0,n)=f(nn,n)
      f1(n,0)=f(n,nn)
      end do

      open(3,file=name_f1)
      do k=0,nn
      write(3,'(257(1x,1pe11.4))')(f1(k,n),n=0,nn)
      end do
      close(3)

      end


      subroutine iso_write1(name_f1,nn,f)
      implicit real*8 (a-h,o-y)
      dimension f(256,256),f1(256,256)
      character*(*) name_f1

      open(3,file=name_f1)
      do k=1,nn
      write(3,'(256(1x,1pe11.4))')(f(k,n),n=1,nn)
      end do
      close(3)
      end

      subroutine iso_write2(name_f1,nn,x,y,f)
      implicit real*8 (a-h,o-y)
      dimension x(256),y(256),f(256,256),f1(256,256)
      character*(*) name_f1

      open(3,file=name_f1)
      do k=1,nn
      do n=1,nn
      write(3,'(3(1x,1pe11.4))')x(k),x(n),f(k,n)
      end do
      end do
      close(3)
      end



***********************Fourier subs*******************
        SUBROUTINE FOUR1(DATA,NN,ISIGN)
                          implicit real*8 (a-h,o-z)
        REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
        DIMENSION DATA(*)
        N=2*NN
        J=1
        DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11      CONTINUE
        MMAX=2
2       IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
        GO TO 2
        ENDIF
        RETURN
        END
        SUBROUTINE REALFT(DATA,N,ISIGN)
                          implicit real*8 (a-h,o-z)
        REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
        DIMENSION DATA(*)
        THETA=6.28318530717959D0/2.0D0/DBLE(N)
        C1=0.5
        IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
        ELSE
        C2=0.5
        THETA=-THETA
        ENDIF
        WPR=-2.0D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.0D0+WPR
        WI=WPI
        N2P3=2*N+3
        DO 11 I=2,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11      CONTINUE
        IF (ISIGN.EQ.1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
        ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N,-1)
        ENDIF
        RETURN
        END


!     *******************  This Function generates an aleatory number *********************
      FUNCTION rand48(seed)
c
c                                Fortran 77 function that implements
c                                the random number generator used
c                                by the C utility erand48. It produces
c                                the same values as erand48 for the
c                                three integers seed(3) and a sequence
c                                of random numbers that agree in all but
c                                the least significant digits. 
c                                (Very small discrepancies in the least 
c                                significant digits are produced by the 
c                                different truncations used by Fortran
c                                and C.)   
c
c                                        Claudio Rebbi - Boston University
c                                                        May 1992 
c
      INTEGER*2 seed(3)
      INTEGER*4 i1,i2,i3,i11,i21,i31,i12,i22,i13
      INTEGER*4 E66D,DEEC,FFFF
      PARAMETER(E66D=58989, DEEC=57068, FFFF=65535)
      REAL*8 rand48

      i1=seed(1)
      IF(i1.LT.0) i1=i1+65536
      i2=seed(2)
      IF(i2.LT.0) i2=i2+65536
      i3=seed(3)
      IF(i3.LT.0) i3=i3+65536

        i11=i1*E66D
        i21=i2*E66D
        i31=i3*E66D
        i12=i1*DEEC
        i22=i2*DEEC
        i13=i1*5
        i1=IAND(i11,FFFF)+11
        i11=ISHFT(i11,-16)+ISHFT(i1,-16)
        i1=IAND(i1,FFFF)
        i11=i11+IAND(i21,FFFF)+IAND(i12,FFFF)
        i2=IAND(i11,FFFF)
        i3=ISHFT(i11,-16)+ISHFT(i21,-16)+ISHFT(i12,-16)+
     &     IAND(i31,FFFF)+IAND(i22,FFFF)+IAND(i13,FFFF)
        i3=IAND(i3,FFFF)
c
c       rand48=i3*2**(-16)+i2*2**(-32)
c
        rand48=i3*1.52587890625D-05+i2*2.328306D-10

      IF(i1.GE.32768) i1=i1-65536
      seed(1)=i1
      IF(i2.GE.32768) i2=i2-65536
      seed(2)=i2
      IF(i3.GE.32768) i3=i3-65536
      seed(3)=i3

      RETURN
      END
