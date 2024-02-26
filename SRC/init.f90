!-----------------------------------------------------------------------
      subroutine stifforme
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      allocate (stiffness(1:nbvar*nodelm,1:nbvar*nodelm))
      allocate (ii(1:nodelm**2))
      allocate (jj(1:nodelm**2))
      stiffness = 1
      nstencil  = nodelm**2
      ii = (/1,1,1,2,2,2,3,3,3/)
      jj = (/1,2,3,1,2,3,1,2,3/)

!-----------------------------------------------------------------------
      end subroutine stifforme
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine init
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      call stifforme

      allocate (U0(1:nbvar*nbrtotNodes))
      allocate (rhs(1:nbvar*nummaxnodes))
      allocate (stencil(1:nbrNodes,1:maxstencil))
      allocate (ia(1:nbrNodes*nbvar+1))
      allocate (ja(1:nbrNodes*maxstencil*nbvar))
      allocate (mat(1:nbrNodes*maxstencil*nbvar))

      !allocate (elem(1:4,1:nbrElem))
      !allocate (node(1:2,1:nbrnodes))
      !allocate (front(1:3,1:nbrfront))
      mat = 0.;   ia = 0;    ja = 0

      rhs = 0.
      U0  = 1.

!-----------------------------------------------------------------------
      end subroutine init
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine deallocate_end
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      if (allocated(U0))   deallocate(U0)
      if (allocated(node)) deallocate(node)
      if (allocated(elem)) deallocate(elem)
      if (allocated(stencil)) deallocate(stencil)
      if (allocated(rhs)) deallocate(rhs)
      if (allocated(ia)) deallocate(ia)
      if (allocated(ja)) deallocate(ja)
      if (allocated(stiffness)) deallocate(stiffness)
      if (allocated(mat)) deallocate(mat)
      if (allocated(nodCL)) deallocate(nodCL)
      if (allocated(nodeglob)) deallocate(nodeglob)
      if (allocated(frontglob)) deallocate(frontglob)
      if (allocated(elemglob)) deallocate(elemglob)
      if (allocated(node)) deallocate(node)
      if (allocated(front)) deallocate(front)
      if (allocated(elem)) deallocate(elem)
      if (allocated(loc2glob)) deallocate(loc2glob)

!-----------------------------------------------------------------------
      end subroutine deallocate_end
!-----------------------------------------------------------------------
