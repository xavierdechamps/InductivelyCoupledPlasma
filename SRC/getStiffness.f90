!-----------------------------------------------------------------------
      subroutine getStiffness()
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      integer(ki) :: k,l,i,j,inod,jnod,m,n
      real(kr)    :: stiffmat(1:nbvar*nodelm,1:nbvar*nodelm)
      real(kr)    :: matK5   (1:nbvar*nodelm,1:nbvar*nodelm)
            
      call getEcoils()
          
      do k=1,nbrElem
         call getStiffLoc(k,stiffmat,matK5)
                  
         do i=1,nodelm
            inod = elem(k,i)
            do j=1,nodelm
               jnod = elem(k,j)
               do m=1,nbvar
                  do n=1,nbvar                  
                     call addvalue_csr(inod,jnod,m,n,stiffmat((i-1)*nbvar+m,(j-1)*nbvar+n))
                  enddo
               enddo
               
               ! add contribution of the element to the RHS
               ! Ecoils is purely imaginary, thus pick only corresponding 
               ! terms in the matK5 matrix
               rhs((inod-1)*nbvar+1) = rhs((inod-1)*nbvar+1) - &
     &             matK5((i-1)*nbvar+1,(j-1)*nbvar+2) * Ecoils(jnod)
               
            end do
            
         end do
      end do
      
!-----------------------------------------------------------------------
      end subroutine getStiffness
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine getStiffLoc(ielm,stiff,matK5)
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      
      integer(ki),intent(in)  :: ielm
      real(kr),intent(out)    :: stiff(1:nodelm*nbvar,1:nodelm*nbvar)
      real(kr)    :: xa,xb,xc,ya,yb,yc,surf,tau,ri,rj
      real(kr)    :: n1r,n2r,n3r,n1z,n2z,n3z,r123,det,h,ksi
      real(kr)    :: nr(3),nz(3),na(3),nb(3),nc(3),nd(3)
      real(kr)    :: a,b,c,d,r1,r2,r3,sigmaloc
      real(kr)    :: matK1(1:nodelm*nbvar,1:nodelm*nbvar)
      real(kr)    :: matK3(1:nodelm*nbvar,1:nodelm*nbvar)
      real(kr)    :: matK4(1:nodelm*nbvar,1:nodelm*nbvar)
      real(kr)    :: matK5(1:nodelm*nbvar,1:nodelm*nbvar)
      integer(ki) :: ii,jj,plasma
      real(kr)    :: delta(3,3)
      
      stiff = zero
      matK1 = zero
      matK3 = zero
      matK4 = zero
      matK5 = zero
      delta = zero
      
      ! CALL isPlasma(ielm,plasma)
      ! CALL getSigma(sigmaloc,plasma)
      sigmaloc = sigma_in(ielm)
      
      r1 = node(elem(ielm,1),2)
      r2 = node(elem(ielm,2),2)
      r3 = node(elem(ielm,3),2)
      n1z = r3 - r2! -r2+r3
      n2z = r1 - r3 ! -r3+r1
      n3z = r2 - r1 ! -r1+r2
      n1r = -node(elem(ielm,3),1) + node(elem(ielm,2),1) ! -z3+z2
      n2r = -node(elem(ielm,1),1) + node(elem(ielm,3),1) ! -z1+z3
      n3r = -node(elem(ielm,2),1) + node(elem(ielm,1),1) ! -z2+z1 
      nz(1:3) = (/n1z,n2z,n3z/)
      nr(1:3) = (/n1r,n2r,n3r/)
      det = (n3r*n2z)-(n2r*n3z)
      surf = abs(det)*0.5
      r123 = r1 + r2 + r3
      
      ! Shape functions at the quadrature points
      na(1:3) = (/1.0d00 , 1.0d00, 1.0d00/) / 3.0d00
      nb(1:3) = (/0.2d00 , 0.2d00, 0.6d00/) 
      nc(1:3) = (/0.2d00 , 0.6d00, 0.2d00/) 
      nd(1:3) = (/0.6d00 , 0.2d00, 0.2d00/) 
      
      ! Weights for the quadrature of Ni Nj / r
      a = - 81.0d00/(96.0d00*r123)
      b = 25.0d00/(96.0d00*(0.20d00*r1+0.20d00*r2+0.60d00*r3))
      c = 25.0d00/(96.0d00*(0.20d00*r1+0.60d00*r2+0.20d00*r3))
      d = 25.0d00/(96.0d00*(0.60d00*r1+0.20d00*r2+0.20d00*r3))
      
      DO ii=1,3
        delta(ii,ii) = 1.0d00
        ri = node(elem(ielm,ii),2)
        DO jj=1,3
          rj = node(elem(ielm,jj),2)
          
          matK1(2*ii-1,2*jj-1) = nr(ii)*nr(jj)
          matK1(2*ii  ,2*jj  ) = matK1(2*ii-1,2*jj-1)
          
          matK3(2*ii-1,2*jj-1) = nz(ii)*nz(jj)
          matK3(2*ii  ,2*jj  ) = matK3(2*ii-1,2*jj-1)
          
          matK4(2*ii-1,2*jj-1) = a*na(ii)*na(jj)+b*nb(ii)*nb(jj)+ &
     &                           c*nc(ii)*nc(jj)+d*nd(ii)*nd(jj)
          matK4(2*ii  ,2*jj  ) = matK4(2*ii-1,2*jj-1)
          
          matK5(2*ii  ,2*jj-1) = (r123 + ri + rj)*(1.0d00 + delta(ii,jj))
          matK5(2*ii-1,2*jj  ) = - matK5(2*ii  ,2*jj-1)
        ENDDO 
      ENDDO
      matK1 = - matK1 * r123 / (12.0d00 * surf)
      matK3 = - matK3 * r123 / (12.0d00 * surf)
      matK4 = - matK4 * abs(det)
      matK5 = - matK5 * omega * mu0 * sigmaloc * surf / 60.0d00
      
      ! write(*,*) "k1+k3"
      ! write(*,*) matk1+matk3
      
      stiff = matK1 + matK3 + matK4 + matK5

!-----------------------------------------------------------------------
      end subroutine getStiffLoc
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine isPlasma(ielm,plasma)
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      
      integer(ki) :: ielm,plasma
      
      plasma = 0
      if (elem(ielm,5) .eq. 2) plasma = 1
!-----------------------------------------------------------------------
      end subroutine isPlasma
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine getSigma(sig,plasma)
!-----------------------------------------------------------------------
      use module_icp, only :kr,ki,sigma,zero
      implicit none
      
      real(kr) :: sig
      integer(ki) :: plasma
      
      sig = zero
            
      if (plasma.eq.1) sig = sigma

!-----------------------------------------------------------------------
      end subroutine getSigma
!-----------------------------------------------------------------------