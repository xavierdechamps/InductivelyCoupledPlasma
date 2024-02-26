!-----------------------------------------------------------------------
      subroutine setBC()
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      integer(ki) :: i,idof,inode
      
      ! Impose Eind     = 0     on axis
      !      dEind / dn = 0     on far field (nothing special to do)
      !      dEind / dn = value on far field
      do i=1,nbrFront
            
         if (front(i,3).eq.CLTable(1,1)) then ! far field dEdn imposed
           call setBCfarfield(i)
           
!         if (front(i,3).eq.CLTable(1,2) .or. front(i,3).eq.CLTable(1,1)) then ! E=0 on all boundary edges
         else if (front(i,3).eq.CLTable(1,2)) then ! edge on axis
          ! First node of the edge element
           inode = front(i,1)
           ! real part
           idof = ( inode - 1 )*nbvar+1
           call setdiag_csr(idof)
           rhs(idof) = zero
           
           ! imaginary part
           idof = idof+1
           call setdiag_csr(idof)
           rhs(idof) = zero
           
           ! Second node of the edge element
           inode = front(i,2)
           ! real part
           idof = ( inode - 1 )*nbvar+1
           call setdiag_csr(idof)
           rhs(idof) = zero
           
           ! imaginary part
           idof = idof+1
           call setdiag_csr(idof)
           rhs(idof) = zero
           
         endif
                  
      end do
            
!-----------------------------------------------------------------------
      end subroutine setBC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine setBCfarfield(indFront)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      integer(ki) :: indFront
      integer(ki) :: node1,node2,i,j,m
      real(kr)    :: rc,zc,length,l2,r1,r2,z1,z2,nr,nz,sumquad
      real(kr)    :: rd(4),zd(4) ! for the numerical quadrature
      real(kr)    :: ksi(4),wd(4),shapefunc(2,4) ! for the numerical quadrature
      real(kr)    :: fac1(4),fac2(4)
      
      node1 = front(indFront,1)
      node2 = front(indFront,2)
      
      r1 = node(node1,2)
      r2 = node(node2,2)
      z1 = node(node1,1)
      z2 = node(node2,1)
      nr = z1 - z2
      nz = r2 - r1
      length = sqrt( nr*nr + nz*nz )
      nr = nr / length
      nz = nz / length
      
      ksi = (/ -0.8611363116 , -0.3399810436 , 0.3399810436 , 0.8611363116 /)
      wd  = (/  0.3478548451 ,  0.6521451548 , 0.6521451548 , 0.3478548451 /)
      
      shapefunc(1,1:4)  = 0.5d00 * (1.0d00 - ksi)
      shapefunc(2,1:4)  = 0.5d00 * (1.0d00 + ksi)
            
      rd = r1 + 0.5d00 * (1.0d00+ksi)*(r2-r1)
      zd = z1 + 0.5d00 * (1.0d00+ksi)*(z2-z1)
      
      zc = coils(2,1)
      zd = zd - zc ! dipole centered on middle coil
      
      ! nz contribution factor
      fac1 = -1.5d00 * wd * length * nz * rd * zd / ( rd*rd + zd*zd )
      ! nr contribution factor
      fac2 = 0.5d00 * wd * length * nr * ( zd*zd - 2.0d00*rd*rd ) / ( rd*rd + zd*zd )
      
      do i=1,2
        node1 = front(indFront,i)
        do j=1,2
          node2 = front(indFront,j)
          sumquad = sum( shapefunc(i,1:4) * shapefunc(j,1:4) * (fac1+fac2) )
          
          do m=1,nbvar
            call addvalue_csr(node1,node2,m,m,sumquad)
          enddo
        enddo
      enddo
      
!-----------------------------------------------------------------------
      end subroutine setBCfarfield
!-----------------------------------------------------------------------