!-----------------------------------------------------------------------
    subroutine getEcoils()
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      
      integer(ki) :: i,c
      real(kr)    :: m,r,Rc,z,Zc,E,K,tmp
      
      Ecoils = zero
      DO i=1,nbrNodes
         z = node(i,1)
         r = node(i,2)
         ! Ec is null on the axis
         if (abs(r).le.eps) cycle
         
         ! Loop on every coils in the domain
         DO c=1,size(coils,1)
            Zc = coils(c,1)
            Rc = coils(c,2)
            
            IF (abs(r-Rc).gt.eps .and. abs(z-Zc).gt.eps) THEN
               m = 4.0d00*r*Rc/((r+Rc)*(r+Rc) + (z-Zc)*(z-Zc))
               CALL CElliptic(m,K,E)
               tmp = omega*mu0*Icoil*sqrt(Rc/r)*((1.0d00-m*0.5d00)*K-E)/(pi*sqrt(m))
               Ecoils(i) = Ecoils(i) + tmp
            ELSE
               Ecoils(i) = 1.0E16
            ENDIF
            
         ENDDO
      ENDDO
            
!-----------------------------------------------------------------------
    end subroutine getEcoils
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    subroutine CElliptic(m,K,E)  
!-----------------------------------------------------------------------
      use module_icp,only:kr,ki,pi
      implicit none
      
      real(kr)  m,alpha,E,K,A,B,C,A_p,B_p,C_0,suma
      integer(ki) j,N   
      
      N=100
      alpha=asin(sqrt(m))
      A_p=1.0d00
      B_p=cos(alpha)
      C_0=sin(alpha)
      suma=0.0d00
      do j=1,N
          A=(A_p+B_p)*0.5d00
          B=dsqrt(A_p*B_p)
          C=(A_p-B_p)*0.5d00
          suma=suma+2**(j)*C**2
          A_p=A
          B_p=B
      end do
      K=pi/(2*A)
      E=(1.0d00 - 0.5d00*(C_0**2+suma))*K
!-----------------------------------------------------------------------
    end Subroutine CElliptic
!-----------------------------------------------------------------------