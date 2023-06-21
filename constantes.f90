module constantes
implicit none
!integer,parameter::pr=selected_real_kind(p=8,r=200),pint=selected_int_kind(r=8)  
integer*4::N,i,nn 
real*8::pi=acos(-1.0d0),cfl,tf,L,inicio,final,dx,dt,t,ni,c,L2
complex*16::img=cmplx(0.0d0,1.0d0)
real*8, allocatable, dimension(:)::x
complex*16, allocatable,dimension (:)::u,kx,ua,uu,rhs,tnl,B,AUX,ufis
character(9)::filename  


end module constantes
