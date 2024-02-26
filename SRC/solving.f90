!-----------------------------------------------------------------------
      subroutine resolution(ok)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      INCLUDE 'dmumps_struc.h'
      include "mpif.h"

      TYPE (DMUMPS_STRUC) mumps_par

      integer(ki) :: n,numnzfonc,ierr,k,l
      integer(ki), allocatable :: iaa(:),jaa(:) ! vecteur contenant les indices de ligne pour le format coordinate
      integer     :: ok
      integer(ki), allocatable :: loc2globall(:,:),numofnodes(:)
      real(kr), allocatable    :: rhsall(:,:)

      ok = 1
      n = nbrtotNodes*nbvar

      mumps_par%COMM = mpi_comm_world

      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 1
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 100
      endif

      mumps_par%ICNTL(3)  = -1 ! verbosity level : uncomment this line if no verbosity is desired
      mumps_par%ICNTL(7)  = 5  ! metis
      mumps_par%ICNTL(5)  = 0  ! assembled format
      mumps_par%ICNTL(18) = 3  ! the distributed matrix is given by the user

      if (irank.eq.0) write(*,*) "Resolution of the system with MUMPS-parallel, version ",mumps_par%VERSION_NUMBER
      if (irank.eq.0) write(*,*) "The matrix is kept distributed"

      allocate(iaa(1:numnz))
      allocate(jaa(1:numnz))

      call sampletime(time1)
      if (irank.eq.0) write(*,*) "    Transforming CSR to coordinate format"
      call csrcoo(nbrnodes*nbvar,numnz,ia,iaa,ierr)
      call sampletime(time2)
      !call time_display

      do k=1,numnz
         ! iaa(k) = loc2glob(iaa(k))
         ! jaa(k) = loc2glob(ja(k))
      enddo

      if (irank.eq.0) mumps_PAR%N = n
      mumps_PAR%NZ_loc = numnz
      allocate( mumps_par%IRN_loc ( mumps_par%NZ_loc ) )
      allocate( mumps_par%JCN_loc ( mumps_par%NZ_loc ) )
      allocate( mumps_par%A_loc   ( mumps_par%NZ_loc ) )
      if (irank.eq.0) allocate( mumps_par%RHS     ( mumps_par%N  ) )

      mumps_par%IRN_loc = iaa(1:numnz)
      mumps_par%JCN_loc = jaa(1:numnz)
      mumps_par%A_loc   = mat(1:numnz)

      allocate(loc2globall(nummaxnodes,0:nproc-1))
      allocate(rhsall(nummaxnodes,0:nproc-1))
      loc2globall = 0
      rhsall   = 0.d0

      ! call mpi_allgather(loc2glob(1),nummaxnodes,mpi_integer,loc2globall(1,0),nummaxnodes,mpi_integer,mpi_comm_world,ierr)

      if (kr==4) then
         call mpi_allgather(rhs(1),nummaxnodes,mpi_real,rhsall(1,0),nummaxnodes,mpi_real,mpi_comm_world,ierr)
      elseif (kr==8) then
         call mpi_allgather(rhs(1),nummaxnodes,mpi_double_precision,rhsall(1,0),nummaxnodes,mpi_double_precision,mpi_comm_world,ierr)
      else
         write(*,*) "Unrecognized real definition in transfert of data"
         ok = 0
         goto 50
      endif

      allocate(numofnodes(0:nproc-1))
      call mpi_allgather(nbrnodes,1,mpi_integer,numofnodes,1,mpi_integer,mpi_comm_world,ierr)

      if (irank.eq.0) then
         do k = 0,nproc-1
            do l = 1,numofnodes(k)
               mumps_par%RHS(loc2globall(l,k)) = mumps_par%RHS(loc2globall(l,k)) + rhsall(l,k)
            enddo
         enddo
      endif

      !//////////////////////////////////
      deallocate(iaa,jaa)
      deallocate(loc2globall,numofnodes)
      deallocate(rhsall)
      !//////////////////////////////////

      if (irank.eq.0) write(*,*) "    Analyzing the system"
      call sampletime(time1)
      mumps_par%JOB = 1
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 50
      endif
      call sampletime(time2)
      call time_display

      if (irank.eq.0) write(*,*) "    Factorizing the system"
      call sampletime(time1)
      mumps_par%JOB = 2
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 50
      endif
      call sampletime(time2)
      call time_display

      !//////////////////////////////////
      deallocate(mumps_par%IRN_loc)
      deallocate(mumps_par%JCN_loc)
      deallocate(mumps_par%A_loc)
      !//////////////////////////////////

      if (irank.eq.0) write(*,*) "    Solving the system"
      call sampletime(time1)
      mumps_par%JOB = 3
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 50
      endif
      call sampletime(time2)
      call time_display

      if (mumps_par%MYID.eq.0) then
         U0 = mumps_par%RHS
      endif

      !//////////////////
      if (irank==0) then
         write(*,*) "allocated",mumps_par%INFOG(18),mumps_par%INFOG(19)/nproc
         write(*,*) "effective",mumps_par%INFOG(21),mumps_par%INFOG(22)/nproc
         write(*,*) "solve",mumps_par%INFOG(30),mumps_par%INFOG(31)/nproc
      endif
      !/////////////////

50    continue
      if (irank.eq.0) deallocate(mumps_par%RHS)

100   continue
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)

!-----------------------------------------------------------------------
      end subroutine resolution
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine resolution2(ok)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      INCLUDE 'dmumps_struc.h'
      include "mpif.h"

      TYPE (DMUMPS_STRUC) mumps_par

      integer(ki) :: n,numnzfonc,ierr,k,l,iproc,ielm,inod,ilocalnod,iglobalnod,ifill
      integer(ki) :: ivar,j1,j2,jcol,jnod,irow,jglobalnod,lrow,lcol,iglobalrow,jglobalcol,i,j
      real(kr)    :: matelm
      integer(ki), allocatable :: iaa(:) 
      integer     :: ok
      integer(ki), allocatable :: loc2globall(:,:),numofnodes(:),numofelem(:)
      real(kr), allocatable    :: rhsall(:,:)
      integer(ki) :: numtotnz,nummaxnz,temp
      integer(ki), allocatable :: ia1(:,:),ja1(:,:),stencil2(:,:),ja2(:),ia2(:)
      real(kr), allocatable :: stiffmat1(:,:),stiffmat2(:)

      ok = 1
      n = nbrtotNodes*nbvar

      call mpi_allreduce(numnz,numtotnz,   1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(numnz,nummaxnz,   1,mpi_integer,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(nbrElem,nummaxelm,1,mpi_integer,mpi_max,mpi_comm_world,ierr)

      allocate(ia1      (1:nummaxnodes*nbvar+1,0:nproc-1))
      allocate(ja1      (1:nummaxnz,0:nproc-1))
      allocate(stiffmat1(1:nummaxnz,0:nproc-1))
      allocate(loc2globall(nummaxnodes,0:nproc-1))
      ia1 = 0; ja1 = 0; stiffmat1 = 0.d0; loc2globall = 0

      allocate(numofnodes(0:nproc-1))
      call mpi_allgather(nbrnodes,1,mpi_integer,numofnodes,1,mpi_integer,mpi_comm_world,ierr)
      allocate(numofelem(0:nproc-1))
      call mpi_allgather(nbrelem,1,mpi_integer,numofelem,1,mpi_integer,mpi_comm_world,ierr)

      call mpi_allgather(ia(1),nummaxnodes*nbvar+1,mpi_integer,ia1(1,0),nummaxnodes*nbvar+1,mpi_integer,mpi_comm_world,ierr)
      call mpi_allgather(ja(1),nummaxnz,mpi_integer,ja1(1,0),nummaxnz,mpi_integer,mpi_comm_world,ierr)
      ! call mpi_allgather(loc2glob(1),nummaxnodes,mpi_integer,loc2globall(1,0),nummaxnodes,mpi_integer,mpi_comm_world,ierr)
      if (kr==4) then
         call mpi_allgather(mat(1),nummaxnz,mpi_real,stiffmat1(1,0),nummaxnz,mpi_real,mpi_comm_world,ierr)
      elseif (kr==8) then
         call mpi_allgather(mat(1),nummaxnz,mpi_double_precision,stiffmat1(1,0),nummaxnz,mpi_double_precision,mpi_comm_world,ierr)
      else
         write(*,*) "Unrecognized real definition in transfert of data"
         ok = 0
         goto 150
      endif

      if (irank==0) then
         allocate(stencil2(1:nbrtotnodes,1:maxstencil))
         allocate(ia2(1:nbrtotnodes*nbvar+1))
         allocate(ja2(1:numtotnz))
         allocate(stiffmat2(1:numtotnz))
         stencil2 = 0; ia2 = 0; ja2 = 0; stiffmat2 = 0.d0

         call createcsr(stencil2,elemglob(1:nodelm,1:nbrtotelem),stiffness,ia2,ja2,nbrtotnodes,nbrtotelem,numtotnz)

         do iproc=0,nproc-1
            do inod=1,numofnodes(iproc)
               do ivar=1,nbvar
                  irow = inod*nbvar+nbvar-1
                  j1 = ia1(irow,iproc)
                  j2 = ia1(irow+1,iproc)
                  do j = j1,j2-1
                     jcol = ja1(j,iproc)
                     jnod = jcol/nbvar
                     iglobalnod = loc2globall(inod,iproc)
                     jglobalnod = loc2globall(jnod,iproc)
                     lrow = irow - inod*nbvar
                     lcol = jcol - jnod*nbvar
                     iglobalrow = iglobalnod*nbvar + lrow
                     jglobalcol = jglobalnod*nbvar + lcol
                     matelm = stiffmat1(j,iproc)
                     do i=ia2(iglobalrow),ia2(iglobalrow+1)-1
                        if (ja2(i).eq.jglobalcol) then
                           stiffmat2(i) = stiffmat2(i) + matelm
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      call mpi_barrier(mpi_comm_world,ierr)

      mumps_par%COMM = mpi_comm_world

      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 1
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 100
      endif

      mumps_par%ICNTL(3) = -1  ! verbosity level : decomment this line if no verbosity is desired
      mumps_par%ICNTL(7) = 5 ! metis

      if (irank.eq.0) write(*,*) "Resolution of the system with MUMPS-parallel, version ",mumps_par%VERSION_NUMBER
      if (irank.eq.0) write(*,*) "The matrix is assembled on proc 0"

      allocate(iaa(1:numtotnz))

      call sampletime(time1)
      if (irank.eq.0) then
         write(*,*) "    Transforming CSR to coordinate format"
         call csrcoo(nbrtotnodes*nbvar,numtotnz,ia2,iaa,ierr)
      endif
      call sampletime(time2)
      !call time_display

      if (irank.eq.0) then
         mumps_PAR%N = n
         mumps_PAR%NZ = numtotnz
         allocate( mumps_par%IRN ( mumps_par%NZ ) )
         allocate( mumps_par%JCN ( mumps_par%NZ ) )
         allocate( mumps_par%A   ( mumps_par%NZ ) )
         allocate( mumps_par%RHS ( mumps_par%N  ) )

         mumps_par%IRN = iaa(1:numtotnz)
         mumps_par%JCN = ja2(1:numtotnz)
         mumps_par%A   = stiffmat2(1:numtotnz)
      endif

      !////////////////////////////////////
      if (irank==0) then
         deallocate(ja2,ia2,stencil2)
         deallocate(stiffmat2)
      endif
      !////////////////////////////////////

      allocate(rhsall(1:nummaxnodes*nbvar,0:nproc-1))
      rhsall = 0.d0

      if (kr==4) then
         call mpi_allgather(rhs(1),nummaxnodes*nbvar,mpi_real,rhsall(1,0),nummaxnodes*nbvar,mpi_real,mpi_comm_world,ierr)
      elseif (kr==8) then
         call mpi_allgather(rhs(1),nummaxnodes*nbvar,mpi_double_precision,rhsall(1,0),nummaxnodes*nbvar,mpi_double_precision,mpi_comm_world,ierr)
      else
         write(*,*) "Unrecognized real definition in transfert of data"
         ok = 0
         goto 50
      endif

      if (irank.eq.0) then
         do k = 0,nproc-1
            do l = 1,numofnodes(k)
               mumps_par%RHS(loc2globall(l,k)) = mumps_par%RHS(loc2globall(l,k)) + rhsall(l,k)
            enddo
         enddo
      endif

      if (irank.eq.0) write(*,*) "    Analyzing the system"
      call sampletime(time1)
      mumps_par%JOB = 1
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 50
      endif
      call sampletime(time2)
      call time_display

      if (irank.eq.0) write(*,*) "    Factorizing the system"
      call sampletime(time1)
      mumps_par%JOB = 2
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 50
      endif
      call sampletime(time2)
      call time_display

      !////////////////////////////////////
      if (irank.eq.0) then
         deallocate(mumps_par%A)
         deallocate(mumps_par%IRN,mumps_par%JCN)
      endif
      !////////////////////////////////////

      if (irank.eq.0) write(*,*) "    Solving the system"
      call sampletime(time1)
      mumps_par%JOB = 3
      CALL DMUMPS(mumps_par)
      if (mumps_par%INFO(1).ne.0) then
         ok = 0
         write(*,'(a35,I2)') "Problem calling MUMPS : INFO(1) = ",mumps_par%INFO(1)
         goto 50
      endif
      call sampletime(time2)
      call time_display

      if (mumps_par%MYID.eq.0) then
         U0 = mumps_par%RHS
      endif

      !//////////////////
      if (irank==0) then
         write(*,*) "allocated",mumps_par%INFOG(18),mumps_par%INFOG(19)/nproc
         write(*,*) "effective",mumps_par%INFOG(21),mumps_par%INFOG(22)/nproc
         write(*,*) "solve",mumps_par%INFOG(30),mumps_par%INFOG(31)/nproc
      endif
      !/////////////////

50    continue

      if (irank.eq.0) then
         deallocate(mumps_par%RHS)
      endif

100   continue
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)

150   continue
      deallocate(stiffmat1)
      deallocate(ja1,ia1,loc2globall,numofelem,numofnodes)
      deallocate(iaa)
      deallocate(rhsall)

!-----------------------------------------------------------------------
      end subroutine resolution2
!-----------------------------------------------------------------------
