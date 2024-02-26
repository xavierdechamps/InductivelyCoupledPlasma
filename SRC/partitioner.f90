!-----------------------------------------------------------------------
      subroutine partitioner
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      include "mpif.h"
    
      real(kr), allocatable  :: CoordOfPartall(:,:,:)
      integer :: i,j,k,l,nameLength,iElm
      integer :: iNod,jNod
      integer :: iSharedNod
      integer :: index,iDummy,APART,BPART,ANOD,BNOD
      integer :: iNode1,iNode2,iBC1,iBC2,iNeighbElm,iCount 
      integer :: iLines,numFModes,iBType,numDOF,SEG1Line
      integer :: iAdded, iPart,jPart, argc,SEG1,SEG2,TEMP,SEG2Line
      integer :: iStart, nBFPart
      integer :: elmPart
      integer :: iLen
      integer :: iVar,iDOF,numTotDOF
      integer :: iverb,ierr
      integer, allocatable :: nodePart(:)
      integer, allocatable :: nodeFlags(:),iNodOfPart(:)
      integer, allocatable :: iNodOfPartM1(:,:),iElmOfPart(:)
      integer, allocatable :: iPartitioning(:,:)
      integer, allocatable :: invperm(:), perm(:)
      integer, allocatable :: elemallloc(:,:,:),frontallloc(:,:,:)

      integer(ki), allocatable :: numpartelem(:), numpartnod(:),numpartfront(:)
      character(len=3)  :: int2ch
      external :: int2ch
!-----------------------------------------------------------------------

      iverb = 0
      
      write(*,'(a,a,a,i3,a)') ' Reading ',trim(maillage),' and splitting into ',nproc,' partitions'

      !-----------------------------------------------------------------
      ! Node reordering
      !-----------------------------------------------------------------

      open(unit=11,file='grid.mts',status='unknown')
      write(11,*) nbrtotElem, ' 1'
      do ielm=1,nbrtotElem
        write(11,'(i7,1x,i7,1x,i7)') elemglob(1,ielm),elemglob(2,ielm),elemglob(3,ielm)
      enddo
      close(11)

      call mesh2nodal(iverb)

      !-----------------------------------------------------------------
      ! Partition graph
      !-----------------------------------------------------------------
      
      if (nproc.gt.1) then
        call partnmesh(nproc,iverb)
      else
        open(unit=11,file='grid.mts.epart.'//int2ch(nproc),status='unknown')
        do ielm = 1, nbrtotelem
          write(11,*) '0 ' 
        enddo
        close(11)
      endif

      ! Consider partitioning of elements

      allocate(nodeflags(1:nbrtotNodes))
      allocate(inodofpart(1:nbrtotNodes))
      allocate(inodofpartm1(1:nbrtotNodes,0:nproc-1))
      allocate(ielmofpart(1:nbrtotElem))
      allocate(ipartitioning(1:nbrtotNodes,0:nproc-1))
      allocate(numpartnod(0:nproc-1))
      allocate(numpartelem(0:nproc-1))
      allocate(numpartfront(0:nproc-1))

      allocate(elemallloc(1:nbrtotelem,1:3,0:nproc-1))
      allocate(frontallloc(1:nbrtotfront,1:3,0:nproc-1))
      allocate(coordofpartall(1:nbrtotNodes,1:2,0:nproc-1))

      numpartnod   = 0
      numpartelem  = 0
      numpartfront = 0

      do ipart = 0, nproc-1
        open(unit=11,file='grid.mts.epart.'//int2ch(nproc),status='unknown')

        nodeflags(1:nbrtotNodes)            = 0
        inodofpart(1:nbrtotNodes)           = 0
        inodofpartm1(1:nbrtotNodes,ipart)   = 0
        ielmofpart(1:nbrtotElem)            = 0
        ipartitioning(1:nbrtotNodes,ipart)  = 0

        elemallloc(1:nbrtotelem,1:3,ipart)    = 0
        frontallloc(1:nbrtotfront,1:3,ipart)  = 0
        coordofpartall(1:nbrtotNodes,1:2,ipart) = 0.d0        

        i = 0
        do ielm = 1,nbrtotElem
          read(11,*) elmpart
          if (elmpart.eq.ipart) then
            i = i + 1
            ielmofpart(i)                = ielm
            nodeflags(elemglob(1,ielm))  = 1
            nodeflags(elemglob(2,ielm))  = 1
            nodeflags(elemglob(3,ielm))  = 1
            numpartelem(ipart)           = numpartelem(ipart) + 1
          endif
        enddo
        close(11)

        i = 0
        do inod = 1, nbrtotNodes
          if (nodeflags(inod).eq.1) then
             i = i + 1
             inodofpart(i)                 = inod
             inodofpartm1(inod,ipart)      = i
             coordofpartall(i,1:2,ipart)   = nodeglob(1:2,inod)
             numpartnod(ipart)             = numpartnod(ipart) + 1
          endif
        enddo

        ! Register nodes in partitioning array
        ipartitioning(:,ipart) = nodeflags(:)

        do i = 1,numpartelem(ipart)
           elemallloc(i,1,ipart) = inodofpartm1(elemglob(1,ielmofpart(i)),ipart)
           elemallloc(i,2,ipart) = inodofpartm1(elemglob(2,ielmofpart(i)),ipart)
           elemallloc(i,3,ipart) = inodofpartm1(elemglob(3,ielmofpart(i)),ipart)
        enddo

        nbfpart = 0
        do j = 1,nbrtotfront
           if (nodeflags(frontglob(1,j)).eq.1 .and. nodeflags(frontglob(2,j)).eq.1) then
              numpartfront(ipart) = numpartfront(ipart) + 1
              nbfpart = nbfpart + 1
              frontallloc(nbfpart,1,ipart) = inodofpartm1(frontglob(1,j),ipart)
              frontallloc(nbfpart,2,ipart) = inodofpartm1(frontglob(2,j),ipart)
              frontallloc(nbfpart,3,ipart) = frontglob(3,j)
           endif
        enddo

        open(unit=11,file='aaa.elem'//int2ch(ipart),status='unknown')
        write(11,*) numpartelem(ipart), nbrtotelem
        do i = 1,numpartelem(ipart)
           write(11,*) elemallloc(i,1:3,ipart) ! local numerotation
        enddo
        close(11)
        
        open(unit=11,file='aaa.node'//int2ch(ipart),status='unknown')
        write(11,*) numpartnod(ipart), nbrtotnodes
        do i = 1,numpartnod(ipart)
           write(11,*) inodofpart(i), coordofpartall(i,1:2,ipart) !numglob + coord
        enddo
        close(11)

        open(unit=11,file='aaa.front'//int2ch(ipart),status='unknown')
        write(11,*) numpartfront(ipart), nbrtotfront
        do i = 1,numpartfront(ipart)
           write(11,*) frontallloc(i,1:3,ipart)
        enddo
        close(11)
      enddo
      
#ifdef WINDOWS
      call system('del grid.mts*')
#else
      call system('rm grid.mts*')
#endif
      deallocate(coordofpartall)
      deallocate(elemallloc)
      deallocate(nodeflags,inodofpart,inodofpartm1)
      deallocate(ielmofpart,ipartitioning)
      deallocate(numpartnod,numpartelem,numpartfront)
      deallocate(frontallloc)
!-----------------------------------------------------------------------
      end subroutine partitioner
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine readPartition
!-----------------------------------------------------------------------
        use module_icp
        implicit none

        include "mpif.h"

        integer :: ierr,i
        character(len=3)  :: int2ch
        external :: int2ch
        logical :: filethere

        do
           inquire(file='aaa.elem'//int2ch(irank),exist=filethere)
           if (filethere) exit
        end do

        ! surface elements for the local partition
        open(unit=11,file='aaa.elem'//int2ch(irank),status='unknown')

        read(11,*) nbrelem,nbrtotelem
        call mpi_allreduce(nbrelem,nummaxelm,1,mpi_integer,mpi_max,mpi_comm_world,ierr)

        allocate (elem(1:3,1:nummaxelm))
        elem = 0
        do i = 1,nbrelem
           read(11,*) elem(1:3,i) ! local numerotation
        enddo
        close(11)

        do
           inquire(file='aaa.node'//int2ch(irank),exist=filethere)
           if (filethere) exit
        end do

        ! node global numbering and coordinates 
        open(unit=11,file='aaa.node'//int2ch(irank),status='unknown')
        read(11,*) nbrnodes,nbrtotnodes

        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_allreduce(nbrnodes,nummaxnodes,1,mpi_integer,mpi_max,mpi_comm_world,ierr)

        allocate (node(1:2,1:nbrnodes))
        allocate (loc2glob(1:nummaxnodes))
        node = 0.d0
        loc2glob = 0
        do i = 1,nbrnodes
           read(11,*) loc2glob(i), node(1:2,i) !numglob + coord
        enddo
        close(11)

        do
           inquire(file='aaa.front'//int2ch(irank),exist=filethere)
           if (filethere) exit
        end do

        open(unit=11,file='aaa.front'//int2ch(irank),status='unknown')
        read(11,*) nbrfront,nbrtotfront

        allocate (front(1:3,1:nbrfront))

        do i = 1,nbrfront
           read(11,*) front(1:3,i)
        enddo
        close(11)
        
#ifdef WINDOWS
        call system('del aaa.elem'//int2ch(irank))
        call system('del aaa.node'//int2ch(irank))
        call system('del aaa.front'//int2ch(irank))
#else
        call system('rm aaa.elem'//int2ch(irank))
        call system('rm aaa.node'//int2ch(irank))
        call system('rm aaa.front'//int2ch(irank))
#endif
        
!-----------------------------------------------------------------------
      end subroutine readPartition
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      character(len=3) function int2ch(number)
!-----------------------------------------------------------------------
      implicit none
      
      ! Variables passed through header
      integer number
      
      ! Variables declared locally
      integer   tt, jj, kk, cero
      character cc(3)

!-----------------------------------------------------------------------

      cero  = ichar( '0' )  ! 48 is the ascii code of zero
      tt = number
      kk = 100
      do jj = 1, 3
        cc(jj) = char(cero + tt/kk)
        tt = mod(tt,kk)
        kk = kk / 10
      enddo
      if (cc(1).ne.'0') then
      int2ch = cc(1)//cc(2)//cc(3)
      else if (cc(2).ne.'0') then
      int2ch = cc(2)//cc(3)//' '
      else if (cc(3).ne.'0') then
      int2ch = cc(3)//'  '
      else
      int2ch = '0  '
      endif

!-----------------------------------------------------------------------
      end function int2ch
!-----------------------------------------------------------------------
