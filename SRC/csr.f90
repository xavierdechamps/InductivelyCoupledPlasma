!-----------------------------------------------------------------------
      subroutine createcsr(sstencil,lnk,stiffcsr,iacsr,jacsr,nbrnodes,nbrelem,numnz)
!-----------------------------------------------------------------------
      use module_icp,only : maxstencil,nbvar,irank,nodelm,ki,kr
      implicit none
      
      ! Variables passed through header
      integer(ki) :: numnz,nbrnodes,nbrelem
      integer(ki) :: sstencil(1:nbrNodes,1:maxstencil)
      integer(ki) :: lnk(1:3,1:nbrElem)
      integer(ki) :: stiffcsr(1:nbvar,1:nbvar)
      integer(ki) :: iacsr(*)
      integer(ki) :: jacsr(*)
      
      ! Variables declared locally
      integer(ki) :: ielm
      integer(ki) :: i,j,k,l,inod,jnod,index
      integer(ki) :: stenciltmp(1:maxstencil-1)
      integer(ki), allocatable:: nbneighbours(:)

      ! Generate ia and ja matrices
      sstencil = -1

      do inod = 1,nbrNodes
        sstencil(inod,1) = inod
      enddo

      do ielm=1,nbrElem !loop over elements
        do i=1,nodelm  !loop over nodes in each element
          inod=lnk(i,ielm)
          do j=1,nodelm
            jnod=lnk(j,ielm)
            !inod has jnod in its stencil. 
            !store if not yet included in stencil
            do k = 1,maxstencil
              if (sstencil(inod,k).eq.jnod) then
                exit
                !no need to store, already there.
              else if (sstencil(inod,k).eq.-1) then
                sstencil(inod,k) = jnod
                !store in first available place
                exit
              endif
            enddo
          enddo
        enddo
      enddo

      ! Reorder nodes within stencil
      do inod = 1,nbrNodes
        stenciltmp = sstencil(inod,2:maxstencil)
        call sortcsr(stenciltmp,maxstencil-1)
        sstencil(inod,2:maxstencil) = stenciltmp
      enddo    

      allocate (nbneighbours(1:nbrNodes)); nbneighbours = 0

      do i=1,nbrNodes
        do k =1, maxstencil
          if (sstencil(i,k).eq.-1) exit
        enddo
        nbneighbours(i) = k-1
      enddo

      index   = 0 
      do i = 1, nbrNodes
        do l = 1, nbvar !loops over a single row in the system matrix
          iacsr((i-1)*nbvar+l) = index+1
          do j = 1, nbneighbours(i)
            do k = 1, nbvar
              if (stiffcsr(l,k) .ne. 0) then
                jacsr(index+1) = sstencil(i,j)*nbvar+k-1
                index = index + 1
              endif
            enddo
          enddo
        enddo
      enddo

      iacsr(nbrNodes*nbvar+1) = index+1
      numnz = index
      deallocate(nbneighbours)

!-----------------------------------------------------------------------
      end subroutine createcsr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine sortcsr(a,n)
!-----------------------------------------------------------------------
      use module_icp, only : kr,ki,nbrNodes,nbvar
        implicit none
      
      ! Variables passed through header
      integer(ki) :: n
      integer(ki) :: a(1:n)
      integer(ki) :: b(1:n)
      
      ! Variables declared locally
      integer(ki) :: i,j,jmax,amax

      b = -1

      do i=1,n
        amax = -1
        jmax = 0
        do j=1,n
          if (a(j).gt.amax) then
            amax = a(j)
            jmax = j
          endif
        enddo
        if (amax.gt.-1) then
          b(i) = amax
          a(jmax) = -1
        endif
      enddo

      a = b

!-----------------------------------------------------------------------
      end subroutine sortcsr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine addvalue_csr(rnod,cnod,lrow,lcol,value)
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      
      ! Variables passed through header
      real(kr) :: value
      integer(ki) :: rnod,cnod,lrow,lcol
      
      ! Variables declared locally
      integer(ki) :: i,irow,icol

      ! Remember: arrays start at position 1
      irow=(rnod*nbvar)
      icol=(cnod*nbvar)
      do i=ia(irow),(ia(irow+1)-1)
        if(ja(i).eq.icol) then
          mat(i)=mat(i)+value
          exit
        endif
      enddo

      if(i==ia(irow+1)) then
        write(*,*) 'addvalue_csr: location not found in csr matrix'
        write(*,*) 'irow,icol',irow,icol
      endif

!-----------------------------------------------------------------------
      end subroutine addvalue_csr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine setdiag_csr(idof)
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      
      ! Variables passed through header
      integer(ki),intent(in) :: idof
      
      ! Variables declared locally
      integer :: i,icol,irow

      irow=idof
      do i=ia(irow),ia(irow+1)-1
        icol=ja(i)
        if(irow.eq.icol) then
          mat(i)=1.
        else
          mat(i)=0.
        endif
      enddo
      
!-----------------------------------------------------------------------
      end subroutine setdiag_csr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine rhsset_csr(idof,dval)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      real(kr),intent(in)    :: dval
      integer(ki),intent(in) :: idof
      integer(ki) :: irow

      irow=idof
      rhs(irow) = dval
      
!-----------------------------------------------------------------------
      end subroutine rhsset_csr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine csrcoo (nrow,nzmax,ia,ir,ierr)
!-----------------------------------------------------------------------
      use module_icp, only : kr, ki
      implicit none
 
      integer(ki) :: ir(*),ia(nrow+1)
      integer(ki) :: ierr,nrow,nzmax,nnz,k,i,k1,k2
!----------------------------------------------------------------------- 
!  Compressed Sparse Row      to      Coordinate 
!----------------------------------------------------------------------- 
! converts a matrix that is stored in coordinate format
!  ir into a row general sparse ia format.
!
! on entry: 
!---------
! nrow	= dimension of the matrix.

! ia    = matrix in compressed sparse row format.
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
! 
! on return:
!----------- 
! ir = matrix in coordinate format.
!
! ierr       = integer error indicator.
!         ierr .eq. 0 means normal retur
!         ierr .eq. 1 means that the the code stopped 
!         because there was no space in ao, ir, jc 
!         (according to the value of  nzmax).
!------------------------------------------------------------------------

      ierr = 0
      nnz = ia(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = 1
         write(*,*) "   Erreur lors de la conversion CSR -> COO",nnz,nzmax
         return
      endif

      do i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do k=k1,k2,-1
            ir(k) = i
         enddo
      enddo
      return
!----------------------------------------------------------------------- 
      end subroutine csrcoo
!-----------------------------------------------------------------------
