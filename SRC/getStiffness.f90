!-----------------------------------------------------------------------
      subroutine setCL
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      integer(ki) :: i,noeud

      do i=1,nbrFront
         noeud = front(1,i)
         ! 8 = gauche, 9 = haut, 10 = droite, 11 = bas
         select case (front(3,i))
         case (8) ! gauche
            call setdiag_csr(noeud)
            rhs(noeud) = 1.
         case (9) ! haut
            call setdiag_csr(noeud)
            rhs(noeud) = 1.
         case (10) ! droite
            ! if (typeCLdroite == 2) then
               call setdiag_csr(noeud)
               rhs(noeud) = node(2,noeud)*2.
            ! end if
         case (11) ! bas
            call setdiag_csr(noeud)
            rhs(noeud) = 0.
         case default
            write(*,*) "Unknown type of BC"
         end select

         noeud = front(2,i)
         ! 8 = gauche, 9 = haut, 10 = droite, 11 = bas
         select case (front(3,i))
         case (8) ! gauche
            call setdiag_csr(noeud)
            rhs(noeud) = 1.
         case (9) ! haut
            call setdiag_csr(noeud)
            rhs(noeud) = 1.
         case (10) ! droite
            ! if (typeCLdroite == 2) then
               call setdiag_csr(noeud)
               rhs(noeud) = node(2,noeud)*2.
            ! end if
         case (11) ! bas
            call setdiag_csr(noeud)
            rhs(noeud) = 0.
         case default
            write(*,*) "Unknown type of BC"
         end select
      end do

!-----------------------------------------------------------------------
      end subroutine setCL
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine getStiffness
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      integer(ki) :: k,l,i,j,inod,jnod
      real(kr)    :: stiffmat(1:nbvar*nodelm,1:nbvar*nodelm)

      do k=1,nbrElem
         call getStiffLoc(k,stiffmat)
         l = 1
         do i=1,nodelm
            inod = elem(i,k)
            do j=1,nodelm
               jnod = elem(j,k)
               call addvalue_csr(inod,jnod,ii(l),jj(l),stiffmat(ii(l),jj(l)))
               l = l+1
            end do
         end do

      end do
!-----------------------------------------------------------------------
      end subroutine getStiffness
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine getStiffLoc(i,stiff)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      integer(ki),intent(in)  :: i
      real(kr),intent(out)    :: stiff(1:nodelm*nbvar,1:nodelm*nbvar)
      real(kr)    :: xa,xb,xc,ya,yb,yc,surf,tau,a1,a2,a3,b1,b2,b3,det,h,ksi
      real(kr)    :: matAa(1:nodelm*nbvar,1:nodelm*nbvar)
      real(kr)    :: matB(1:nodelm*nbvar,1:nodelm*nbvar)
      real(kr)    :: matS(1:nodelm*nbvar,1:nodelm*nbvar)

      a1 = node(2,elem(2,i)) - node(2,elem(3,i))
      a2 = node(2,elem(3,i)) - node(2,elem(1,i))
      a3 = node(2,elem(1,i)) - node(2,elem(2,i))
      b1 = node(1,elem(3,i)) - node(1,elem(2,i))
      b2 = node(1,elem(1,i)) - node(1,elem(3,i))
      b3 = node(1,elem(2,i)) - node(1,elem(1,i))
      det = (b3*a2)-(b2*a3)
      h = 2*abs(det)/(abs(a1)+abs(a2)+abs(a3))
      surf = abs(det)*0.5

      matAa(1,1:3) = (/a1*a1+b1*b1,a1*a2+b1*b2,a1*a3+b1*b3/)
      matAa(2,1:3) = (/a1*a2+b1*b2,a2*a2+b2*b2,a2*a3+b2*b3/)
      matAa(3,1:3) = (/a1*a3+b1*b3,a2*a3+b2*b3,a3*a3+b3*b3/)
      matAa = matAa/(2*abs(det))

      matB(1,1:3) = (/a1,a2,a3/)
      matB(2,1:3) = (/a1,a2,a3/)
      matB(3,1:3) = (/a1,a2,a3/)
      matB = matB/6

      tau = zero
      matS = zero
      stiff = zero
      ! if (type_reso /= 1) then
         matS(1,1:3) = (/a1*a1,a1*a2,a1*a3/)
         matS(2,1:3) = (/a2*a1,a2*a2,a2*a3/)
         matS(3,1:3) = (/a3*a1,a3*a2,a3*a3/)
         matS = matS/(2*abs(det))

         ! if (type_reso==2) then
            ! h  =  2*sqrt(surf/pi)
         ! elseif(type_reso==3) then
            ! h = 2*abs(det)/(abs(a1)+abs(a2)+abs(a3));
         ! end if

         ! if (u*h/(6*alpha)>1) then
            ! ksi = 1.
         ! else
            ! ksi = u*h/(6*alpha)
         ! end if

         ! tau = h*ksi/(2*u)

      ! end if

      ! stiff = alpha*matAa + u*matB + tau*u*u*matS

!-----------------------------------------------------------------------
      end subroutine getStiffLoc
!-----------------------------------------------------------------------
