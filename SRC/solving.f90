!-----------------------------------------------------------------------
      subroutine solve(ok)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      INCLUDE 'dmumps_struc.h'
      include "mpif.h"

      TYPE (DMUMPS_STRUC) mumps_par

      integer(ki) :: ierr,k,l
      integer(ki), allocatable :: iaa(:),jaa(:) ! vecteur contenant les indices de ligne pour le format coordinate
      integer     :: ok

      ok = 1

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

      mumps_par%ICNTL(3)  = -1 ! verbosity level : comment this line if no verbosity is desired
      mumps_par%ICNTL(7)  = 5  ! metis
      mumps_par%ICNTL(5)  = 0  ! assembled format
      mumps_par%ICNTL(18) = 3  ! the distributed matrix is given by the user

      ! if (irank.eq.0) write(*,*) "Resolution of the system with MUMPS-parallel, version ",mumps_par%VERSION_NUMBER
      ! if (irank.eq.0) write(*,*) "The matrix is kept distributed"

      allocate(iaa(1:numnz))
      ! allocate(jaa(1:numnz))

      call sampletime(time1)
      if (irank.eq.0) write(*,*) "    Transforming CSR to coordinate format"
      call csrcoo(nbrnodes*nbvar,numnz,ia,iaa,ierr)
      call sampletime(time2)
      call time_display

      if (irank.eq.0) mumps_PAR%N = nbrNodes*nbvar
      mumps_PAR%NZ_loc = numnz
      allocate( mumps_par%IRN_loc ( mumps_par%NZ_loc ) )
      allocate( mumps_par%JCN_loc ( mumps_par%NZ_loc ) )
      allocate( mumps_par%A_loc   ( mumps_par%NZ_loc ) )
      if (irank.eq.0) allocate( mumps_par%RHS     ( mumps_par%N  ) )
            
      mumps_par%IRN_loc = iaa(1:numnz)
      mumps_par%JCN_loc = ja(1:numnz)
      mumps_par%A_loc   = mat(1:numnz)
      mumps_par%RHS     = rhs(1:nbrNodes*nbvar)
            
      ! //////////////////////////////////
      deallocate(iaa)
      ! //////////////////////////////////

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

      ! //////////////////////////////////
      deallocate(mumps_par%IRN_loc)
      deallocate(mumps_par%JCN_loc)
      deallocate(mumps_par%A_loc)
      ! //////////////////////////////////

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

      ! //////////////////
      ! if (irank==0) then
         ! write(*,*) "allocated",mumps_par%INFOG(18),mumps_par%INFOG(19)/nproc
         ! write(*,*) "effective",mumps_par%INFOG(21),mumps_par%INFOG(22)/nproc
         ! write(*,*) "solve",mumps_par%INFOG(30),mumps_par%INFOG(31)/nproc
      ! endif
      ! /////////////////

50    continue
      if (irank.eq.0) deallocate(mumps_par%RHS)

100   continue
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)

!-----------------------------------------------------------------------
      end subroutine solve
!-----------------------------------------------------------------------