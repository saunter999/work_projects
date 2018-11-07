module params
implicit none
real*8::t=1.0,U=4.0,beta=2.0,mu=0.0
integer*4 :: N=6,Ns,Natom,Ntau=500,Nh=10
integer*4,allocatable::rls(:,:)
real*8,allocatable::hls(:),RBZ(:,:),ekls(:),tauls(:)
real*8,allocatable::omg(:,:)
real*8::hb=0.5
real*8::pi=3.14159265358979
complex*16,allocatable::Grtau_aup(:,:,:),Grtau_adn(:,:,:),Grtau_bup(:,:,:),Grtau_bdn(:,:,:)
end module


Program threedHubAFM
use params
implicit none
integer*4::Nexp=2,i,j,k
real*8 temp
!******functions declaration*******!
real*8 Gqp_ktau,Gckup,GcAup

Ns=N**3
Natom=Ns*2
allocate(rls(Ns,3))
allocate(hls(Nh),tauls(Ntau),omg(Nexp,Nh))
allocate(RBZ(Ns,3),ekls(Ns))
allocate(Grtau_aup(Ns,2,Ntau),Grtau_adn(Ns,2,Ntau),Grtau_bup(Ns,2,Ntau),Grtau_bdn(Ns,2,Ntau))
open(unit=90,file="F0.txt")
open(unit=91,file="F1.txt")
open(unit=92,file="F2.txt")
call hmesh()
call tauls_gen()
call RBZ_gen()
call ekdisp()

!!evaluating free energy order by order
call F0(omg(1,:)) 
call F1(omg(2,:)) 
do i=1,Nh
write(90,*)hls(i),omg(1,i)
write(91,*)hls(i),omg(2,i)+omg(1,i)
enddo

do i=1,Nh
call Grtau_gen(hls(i))
!write(92,*)hls(i),Grtau_aup(1,1,Ntau),Grtau_adn(1,1,Ntau)
enddo



end program threedHubAFM


subroutine hmesh()
use params
implicit none
integer*4 i
do i=1,Nh
   hls(i)=-hb+2.0*hb/(Nh-1)*(i-1)
enddo
return
end subroutine

subroutine tauls_gen()
use params
implicit none
integer*4 i
do i=1,Ntau
   tauls(i)=beta*(i-1.)/(Ntau-1)
enddo
return
end subroutine


subroutine RBZ_gen()
!Three unit vector for the primitive cell in real space with AFM spin structure
!b1=pi*array([1,1,-1]),b2=pi*array([1,-1,1]),b3=pi*array([-1,1,1]) 
use params
implicit none
integer*4 i,j,k,ik
ik=1
do i=1,N
   do j=1,N
      do k=1,N
      RBZ(ik,1)=pi*( (i-1.0)/N+(j-1.0)/N-(k-1.0)/N)
      RBZ(ik,2)=pi*( (i-1.0)/N-(j-1.0)/N+(k-1.0)/N)
      RBZ(ik,3)=pi*( -(i-1.0)/N+(j-1.0)/N+(k-1.0)/N)
      ik=ik+1
      enddo
   enddo
enddo
return
end subroutine

subroutine ekdisp()
use params
implicit none
integer*4 i
do i=1,Ns
    ekls(i)=-2.0*t*(cos(RBZ(i,1))+cos(RBZ(i,2))+cos(RBZ(i,3)))
enddo
end subroutine


subroutine Grtau_gen(h)
use params
implicit none
integer*4 i,j,k,m,l,ir
real*8 tau,h
real*8 Gckup,GcAup
complex*16 t1,t2,p1,p2
real*8 f1,f2
Grtau_aup=dcmplx(0.0,0.0)
Grtau_adn=dcmplx(0.0,0.0)
Grtau_bup=dcmplx(0.0,0.0)
Grtau_bdn=dcmplx(0.0,0.0)
ir=1

do i=1,N
   do j=1,N
      do k=1,N
        write(*,*)i,j,k
        rls(ir,1)=i
        rls(ir,2)=j
        rls(ir,3)=k
        do l=1,Ntau
            tau=tauls(l)
            do m=1,Ns
              f1=Gckup(ekls(m),h,tau)
              f2=GcAup(ekls(m),h,tau)
              p1=exp(dcmplx(0.0,(i-1+j-1)*RBZ(m,1)+(i-1+k-1)*RBZ(m,2)+(j-1+k-1)*RBZ(m,3)))
              p2=exp(dcmplx(0.0,RBZ(m,1)+RBZ(m,2)+RBZ(m,3)))
              t1=(f1+f2)*p1
              t2=(f1-f2)*p1
              Grtau_aup(ir,1,l)=Grtau_aup(ir,1,l)+t1
              Grtau_aup(ir,2,l)=Grtau_aup(ir,2,l)+t1*p2
              Grtau_bup(ir,1,l)=Grtau_bup(ir,1,l)+t2/p2
              Grtau_bup(ir,2,l)=Grtau_bup(ir,2,l)+t2

              Grtau_adn(ir,1,l)=Grtau_adn(ir,1,l)+t2
              Grtau_adn(ir,2,l)=Grtau_adn(ir,2,l)+t2*p2
              Grtau_bdn(ir,1,l)=Grtau_bdn(ir,1,l)+t1/p2
              Grtau_bdn(ir,2,l)=Grtau_bdn(ir,2,l)+t1
            enddo
            Grtau_aup(ir,:,l)=Grtau_aup(ir,:,l)/Ns
            Grtau_adn(ir,:,l)=Grtau_adn(ir,:,l)/Ns
            Grtau_bup(ir,:,l)=Grtau_bup(ir,:,l)/Ns
            Grtau_bdn(ir,:,l)=Grtau_bdn(ir,:,l)/Ns
         enddo
         ir=ir+1
      enddo
   enddo
enddo

return
end subroutine

subroutine F0(omg0)
use params
implicit none
integer*4 i,j
real*8 x,y,h
real*8 omg0(Nh)

do i=1,Nh
    omg0(i)=0.d0
    h=hls(i)
    !print *,h,Ns,mu,beta,Natom,U
    do j=1,Ns
       x=sqrt(ekls(j)**2+h**2)-mu
!       print *,ek(j),h
       y=-sqrt(ekls(j)**2+h**2)-mu
       if(x>0) then
          omg0(i)=omg0(i)+log(1.0+exp(-beta*x))
       else
          omg0(i)=omg0(i)-beta*x+log(1.0+exp(beta*x))
       endif
       if (y>0) then
          omg0(i)=omg0(i)+log(1.0+exp(-beta*y))
       else
          omg0(i)=omg0(i)-beta*y+log(1.0+exp(beta*y))
       endif
    enddo
     
    omg0(i)=omg0(i)*(-2.0/beta)*(1./Natom) 
!    omg0(i)=omg0(i)+h**2/U
enddo
return
end subroutine

subroutine F1(omg1)
use params
implicit none
real*8 omg1(Nh),neff_0up,neff_0dn,h
real*8 Gckup,GcAup
integer*4 i,j

!open(unit=93,file="check.txt")
do i=1,Nh
   neff_0up=0.d0
   neff_0dn=0.d0
   h=hls(i)
   do j=1,Ns 
     neff_0up=neff_0up+Gckup(ekls(j),h,beta)+GcAup(ekls(j),h,beta)
     neff_0dn=neff_0dn+Gckup(ekls(j),h,beta)-GcAup(ekls(j),h,beta)
   enddo
   neff_0up=-1.0*neff_0up/Ns
   neff_0dn=-1.0*neff_0dn/Ns
   omg1(i)=U*(neff_0up-0.5+h/U)*(neff_0dn-0.5-h/U)
   omg1(i)=omg1(i)+h**2/U
!   write(93,*)h,neff_0up,neff_0dn
enddo 
return 
end subroutine


function Gqp_ktau(e,tau)
!***return the Green's function of quasiparticle with energy e in (k,tau) representation ***!
use params
implicit none
real*8 Gqp_ktau,e,tau
if (e>0) then
Gqp_ktau=-exp(-e*tau)/(1.0+exp(-e*beta))
else
Gqp_ktau=-exp(-e*(tau-beta))/(1.0+exp(e*beta))
endif

return
end function

function Gckup(ek,h,tau)
use params
implicit none
real*8 Gckup,ek,h,tau
real*8 ekh,Gqp_ktau
ekh=sqrt(ek**2+h**2)
Gckup=1./(2.0*ekh*(ekh+ek))*((ekh+ek)**2*Gqp_ktau(ekh-mu,tau)+h**2*Gqp_ktau(-ekh-mu,tau))

return
end function


function GcAup(ek,h,tau)
use params
implicit none
real*8 GcAup,ek,h,tau
real*8 ekh,Gqp_ktau
ekh=sqrt(ek**2+h**2)
GcAup=h/(2.0*ekh)*(Gqp_ktau(ekh-mu,tau)-Gqp_ktau(-ekh-mu,tau))

return
end function
