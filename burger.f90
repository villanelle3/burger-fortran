program burger
use constantes
implicit none
 
open (unit=01,file='dados.txt',STATUS='unknown',ACTION='READ')

read(1,*)N
read(1,*)cfl
close(01)

allocate(x(N),u(N),kx(N),ua(N),uu(N),rhs(N),tnl(N),B(N),AUX(N),ufis(N))
tf=pi/8.0d0

!Definindo o espaco fisico, dx e x referentes a malha base
L=2.0d0*pi
dx=L/N
ni=0.2d0
c=4.0d0

do i=1,N
	x(i)=dx*(i-1) -pi
end do

!Escrevendo os numeros de onda
do i=1,N/2+1
	kx(i)=2.0d0*pi*(i-1)/(N*dx)
end do
  
do i=N/2+2,N
 	Kx(i)=2.0d0*pi*(i-1-n)/(N*dx)
end do

!Criando a condicao inicial
t=0.0d0
call CI(u,x,t,N,ni,c)

dt = cfl*min(2.0d0/ni*1.0d0/dx**2, dx/c)

open (unit=02,file='CI.dat',STATUS='unknown',ACTION='write')
write(02,*)'#ci'
do i=1,N
	write(02,*)x(i),real(u(i))
end do
close(02)

!inicando LOOP temporal
call CPU_time(inicio)
CALL ZFFT1D(ufis,N,0,B)
CALL ZFFT1D(u,N,0,B)
CALL ZFFT1D(uu,N,0,B)
CALL ZFFT1D(rhs,N,0,B)
CALL ZFFT1D(tnl,N,0,B)

nn=0
do while (t<=(tf+dt)) 	

	!dt = cfl*min(2.0d0/ni*1/dx**2, dx/maxval(real(u)))
	call avanco_temporal
	!t= t+dt
	nn=nn+1
	t=nn*dt
end do

open (unit=03,file='saida.dat',STATUS='unknown',ACTION='write')
write(03,*)'#num'
do i=1,N
	write(03,*)x(i),real(u(i))
end do

!write(*,*)('dt',dt)
!write(*,*)('t',t)
!write(*,*)('tf',tf)
call sa(ua,x,t,N,ni,c)
write(03,*)'#cont'
do i=1,N
	write(03,*)x(i),real(ua(i))
end do

close(03)

L2=0.0d0
do i=1,n
	L2=(real(u(i))-real(ua(i)))**2 + L2
end do
L2 = dsqrt(L2/real(n))
write(*,*)'L2= ',L2


Call CPU_time(final)
 
write(*,*)'Custo(s):' 
write(*,'(f6.3)')final-inicio

contains


subroutine calculo_rhs
use constantes
implicit none

uu=u*u
CALL ZFFT1D(uu,N,-1,B)	
uu=-0.25d0*img*kx*uu


CALL ZFFT1D(u,N,-1,B)	
tnl =img*kx*u

CALL ZFFT1D(u,N,1,B)	
CALL ZFFT1D(tnl,N,1,B)
tnl = -0.5d0*tnl*u

CALL ZFFT1D(tnl,N,-1,B)
CALL ZFFT1D(u,N,-1,B)
tnl = tnl + uu
rhs= tnl + ni*(img**2)*(kx**2)*u


end subroutine calculo_rhs

subroutine avanco_temporal
use constantes
implicit none
complex*16,dimension(n)::k1,k2,k3,k4
real*8,parameter,dimension(6)::alfa=(/ 0.0d0, -0.691750960670d0, -1.727127405211d0, -0.694890150986d0,  -1.039942756197d0, -1.531977447611d0 /), beta=(/ 0.122d0, 0.477263056358d0, 0.381941220320d0, 0.447757195744d0, 0.498614246822d0, 0.186648570846d0 /), gama=(/ 0.0d0, 0.122d0, 0.269115878630d0, 0.447717183551d0, 0.749979795490d0, 0.898555413085d0 /) !coeficientes RK46

	!!do i=1,6
!		call calculo_rhs
	!	AUX=alfa(i)*AUX+dt*rhs
	!	u=u+beta(i)*AUX    
	!	CALL ZFFT1D(u,N,+1,B)
	!end do 

	CALL ZFFT1D(u,N,-1,B)
	call RK(K1,dt,u)
	call RK(K2,t+0.5d0*dt,u+0.5d0*dt*K1)
	call RK(K3,t+0.5d0*dt,u+0.5d0*dt*K2)
	call RK(K4,t+dt,u+dt*K3)

	!CALL ZFFT1D(u,N,-1,B)
	u=u+dt/6.0d0*(k1+2.0d0*k2+2.0d0*k3+k4)
	!u=u+dt*k1
	CALL ZFFT1D(u,N,+1,B)

end subroutine avanco_temporal

subroutine RK(k,dt,y)
use constantes
implicit none 
complex*16,dimension(n),intent(out)::k
complex*16,dimension(n),intent(in)::y
real*8,intent(in)::dt
!1. Calculate the multiplicative of the velocity at physical space, as uu;
!2. Transform the product to spectral space, as uu;
!3. Calculate the derivative and multiply by 1/2 after calculating the advective form,
!4. Transform the field u to spectral space;
!5. Calculate the derivatives, iku^;
!6. Calculate the inverse transform of the velocity derivatives multiplied by the velocity in the physical space.
!7. Transform the derivative from the physical space to spectral space,

	CALL ZFFT1D(y,N,0,B)	
	CALL ZFFT1D(y,N,+1,B)	
	uu=y*y
	CALL ZFFT1D(uu,N,-1,B)	
	uu=0.25d0*img*kx*uu


	CALL ZFFT1D(y,N,-1,B)	
	tnl =img*kx*y

	CALL ZFFT1D(y,N,1,B)	
	CALL ZFFT1D(tnl,N,1,B)
	tnl = 0.5d0*tnl*y

	CALL ZFFT1D(tnl,N,-1,B)
	CALL ZFFT1D(y,N,-1,B)
	tnl = tnl + uu
	k= - tnl - ni*kx**2*y

end subroutine RK

end program burger
