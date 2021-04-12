	  program TurbulentBetweenParallelPlates
c Taha Yaşar Demir / 1881978
c CE-580 - Homework #3	  
	  parameter (mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, rnu, tau
	  common/coef/ A(mx), B(mx), C(mx), D(mx)
	  common/error/u_old(mx),error(mx) 
c     sml = mixinglength , fu = viscous damping function, yp=y+, Ap=A+
c     us = u_star, rnu = kinematic viscosity, tau = wall shear

	  open(1,file='error.dat')
	  open(2,file='yplus.dat')
	  open(3,file='velpr.dat')
	  open(4,file='loglaw.dat')
	  open(6,file='vist.dat')


	  call init 
	  do l=1,100
		call stress(l) ! Evaluate stresses
		call coefficients ! Get Three diagonal system coefficents
        call THOMAS(2,N-1,A,B,C,D) ! Returns the solution in D vector
        u(1) = 0. ! no-slip Boundary Condition
        u(N) = D(N-1) ! Neumann Boundary at centerline
        do k = 2,N-1 ! extract the solution
        	u(k) = D(k)
        enddo 
        if (i.gt.1) call error_cal(l) ! After first iteration calculate error
        do k=1,N
        	u_old(k) = u(k) 
        enddo
        call output(l) ! write solution out in relative files
	  enddo 

	  print*, Cp, u(N) ! get max velocity at centerline

	  close(1)
	  close(2)
	  close(3)
	  close(4)
	  close(6)


	  stop
	  end

	  subroutine init
	  parameter (mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, rnu, tau
	  common/coef/ A(mx), B(mx), C(mx), D(mx)

	  H   = 0.02 ! m
	  vis = 0.001 ! N.s/m^2
	  rho = 1000.0 ! kg/m^3
	  N   = mx
	  print*, "Enter the pressure coef. Cp :"
	  read(*,*) Cp
	  if (Cp.eq.-100.) then
	  	Beta = 0.96
	  elseif(Cp.eq.-1000.) then
	  	Beta = 0.94
	  else
	  	Beta = 0.93
	  endif
	  call makegrid
	  call analytic

	  return
	  end

	  subroutine makegrid
	  parameter (mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, rnu, tau
	  common/coef/ A(mx), B(mx), C(mx), D(mx)
	  real sum

	  sum = 0.0
	  do i=0,N-2
	  	sum = sum+ Beta**i
	  enddo
	  dy(N) = H/sum
	  y(N)  = H

	  do i=N,2,-1
	  	y(i-1)  = y(i)-dy(i)
	  	dy(i-1) = Beta*dy(i)
	  	yc(i)   = (y(i)+y(i-1))/2
	  enddo



	  return
	  end

	  subroutine analytic
	  parameter (mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N

	  do i=1,N
	  	ua(i) = -500*y(i)**2 + 200*y(i)
	  enddo
	  ua(1) = 0.

	  return
	  end

	  subroutine stress(iter)
	  parameter (mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, rnu, tau

	  Ap  = 26
	  rnu = vis/rho
	  tau = -Cp*H	 
	  us  = sqrt(tau/rho)
	  do k=2,N
	   yp(k)  = yc(k)*us/rnu
	   fu(k)  = 1-exp(-yp(k)/Ap)
	   sml(k) = H*(0.14-0.08*(1-(yc(k)/H))**2
     &  	   -0.06*(1-(yc(k)/H))**4)*fu(k)
	   if (iter.eq.1) then
	    vist(k) = 0.5*(rho*sml(k)**2)*(ua(k)-ua(k-1))/dy(k) ! half of the laminar case
	   else
	    vist(k) = 0.5*(vist(k)+(rho*sml(k)**2)*(u(k)-u(k-1))/dy(k)) 
	   endif 
	   vise(k) = (vist(k)+vis) 
	  enddo


	  return
	  end

	  subroutine coefficients
	  parameter (mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, rnu, tau
	  common/coef/ A(mx), B(mx), C(mx), D(mx)

	  do i=2,N-1
	  	if (i.eq.2) then ! no slip Boundary Conditions
	  		A(i) =  0
	  		C(i) = -vise(i+1)/(dy(i+1)*(dy(i+1)+dy(i))/2)
	  		B(i) = (vise(i)/(dy(i)*(dy(i+1)+dy(i))/2)) - C(i)
	  		D(i) = -Cp
	  	elseif (i.eq.N-1) then ! Neumann Boundary Condition u_n=u_n-1
	  		A(i) = -vise(i)/(dy(i)*(dy(i+1)+dy(i))/2)
	  		C(i) =  0
	  		B(i) = -A(i)
	  		D(i) = -Cp 
	  	else
	  		A(i) = -vise(i)/(dy(i)*(dy(i+1)+dy(i))/2)
	  		C(i) = -vise(i+1)/(dy(i+1)*(dy(i+1)+dy(i))/2)
	  		B(i) = -A(i) - C(i)
	  		D(i) = -Cp
	  	endif
	  enddo

	  return
	  end

c-------------------------------------------------------------------
	  subroutine error_cal(ite)
	  parameter (mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, rnu, tau
	  common/coef/ A(mx), B(mx), C(mx), D(mx)
	  common/error/u_old(mx),error(mx)

	  error(ite) = 0.
	  do i=1,N
	  error(ite) = error(ite) + (1/((N-1)*u(N)))*abs(u_old(i)-u(i))
	  enddo

	  return
	  end
c-------------------------------------------------------------------
	  subroutine output(it)
	  parameter(mx=101)
	  common/flow/ rho, H, Cp, vis, vist(mx), vise(mx), u(mx), ua(mx)
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),N
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, rnu, tau
	  common/error/u_old(mx),error(mx)
	  real law(mx),up(mx)

	  write(1,*) it,error(it)
	  if(it.eq.1) write(3,*) u(1),y(1),ua(1)
	  if(it.eq.100) then
	  	do i=2,N
 	  	 law(i) = (1/0.41)*log(yp(i))+5.1
 	  	 up(i)  = 0.5*(u(i)+u(i-1))/us	
	  	 write(2,*) yp(i),up(i)
	  	 write(3,*) u(i),y(i),ua(i)
	  	 write(4,*) yp(i),law(i)
		 write(6,*) vist(i),y(i)
	  	enddo
	  endif	

	  return
	  end

c-------------------------------------------------------------------
      subroutine THOMAS(il,iu,aa,bb,cc,ff)
c............................................................
c Solution of a tridiagonal system of n equations of the form
c  A(i)*x(i-1) + Bb(i)*x(i) + C(i)*x(i+1) = R(i)  for i=il,iu
c  the solution X(i) is stored in F(i)
c  A(il-1) and C(iu+1) are not used
c  A,Bb,C,R are arrays to bbe provided bby the user
c............................................................
      parameter (mx=101)
      dimension  aa(mx),bb(mx),cc(mx),ff(mx),tmp(mx)

      tmp(il)=cc(il)/bb(il)
      ff(il)=ff(il)/bb(il)
      ilp1 = il+1
      do i=ilp1,iu
         z=1./(bb(i)-aa(i)*tmp(i-1))
         tmp(i)=cc(i)*Z
         ff(i)=(ff(i)-aa(i)*ff(i-1))*z
      enddo
      iupil=iu+il
      do ii=ilp1,iu
         i=iupil-ii
         ff(i)=ff(i)-tmp(i)*ff(i+1)
      enddo
      return
      end
