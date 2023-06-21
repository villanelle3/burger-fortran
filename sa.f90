subroutine sa (ua,x,t,N, ni, c)
implicit none
integer*4, intent (in):: N
real*8, intent (in)::t,ni,c
real*8, dimension(N), intent (in)::x
complex*16, dimension(N), intent (inout)::ua
real*8::pi=acos(-1.0d0),e,der
integer*4::i,j
real*8, dimension(N):: ub


e = 0.0d0
der = 0.0d0
ub = 0.0d0
ua=0.0d0

do i=1,N
	e=0
	der=0
	do j=-200,200,1
		e=exp(-(((x(i)-c*t)-(2.0d0*j+1.0d0)*pi)**2.0d0)/(4.0d0*ni*(t+1.0d0)))+e
		der=exp(-(((x(i)-c*t)-(2.0d0*j+1.0d0)*pi)**2.0d0)/(4.0d0*ni*(t+1.0d0)))* (-(2.0d0*x(i)-2.0d0*c*t-4.0d0*j*pi-2.0d0*pi))/(4.0d0*ni*(t+1.0d0)) + der
	end do
	ua(i)=(c-2.0d0*ni*der/e)
enddo          
             
   
end subroutine
