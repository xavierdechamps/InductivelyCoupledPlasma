PROGRAM main
    USE module_icp
    IMPLICIT NONE
    INCLUDE "mpif.h"

    ! Local parameters
    INTEGER(ki) :: ok,ierr
    
    write(*,*) "inside the program" 
    
    IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
        WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
        WRITE(*,*) "icp name_of_parameter_file"
        GOTO 200
    ENDIF 
    CALL get_command_argument(1,maillage)
    
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    ! ------- Données du problème ----------
    ! maillage     = "StructRafVertHor_fine.msh"
    file_gmsh    = "essai.msh"
    file_dat     = "essai.dat"
    file_restart = "essai.dat"
    nbvar        = 2
    
    eps = 1.E-12
    ok = 1
    ! --------------------------------------
    ! if (irank.eq.0) then
       ! call lecture_gmsh(ok)
    ! endif

! if (ok == 0) then
  ! write(*,*) "The program hasn't started because of a problem during the reading of the mesh"
  ! goto 200
! else
   ! if (irank.eq.0) call partitioner
   
   ! call mpi_barrier(mpi_comm_world,ierr)
   ! call readPartition

   ! call init
   ! if (irank.eq.0) write(*,*) "Constructing the CSR structure"
   ! call createcsr(stencil,elem,stiffness,ia,ja,nbrnodes,nbrelem,numnz) 

   ! if (irank.eq.0) write(*,*) "Constructing the stiffness matrix"
   ! call getStiffness

   ! if (irank.eq.0) write(*,*) "Imposing the boundary conditions"
   ! call setCL

   ! if (type_solver==1) then
      ! call resolution(ok)
   ! else
      ! call resolution2(ok)
   ! endif

   ! if (ok.eq.1 .and. irank.eq.0) call ecriture_gmsh
   
! endif

200 continue
    if (irank.eq.0) write(*,*) "End of the simulation"

    call deallocate_end
    call mpi_finalize( ierr )

end program


! *******************************************************************
subroutine sampletime(counter)
  use module_icp, only : kr,ki
      implicit none
      
      ! variables passed through header
      integer(ki) ::counter
      
      ! variables declared locally
      integer(ki) ::rate, contmax
 
      ! Determine CPU time
      call system_clock(counter, rate, contmax )
!-----------------------------------------------------------------------
      end subroutine sampletime
!-----------------------------------------------------------------------

! *******************************************************************
SUBROUTINE time_display
use module_icp
implicit none

real(kr) :: time
integer(ki) :: job, cont3
integer(ki) :: rate, contmax, itime

call system_clock(cont3, rate, contmax )
if (time2 .ge. time1) then
   itime=time2-time1
else
   itime=(contmax - time1) + (time2 + 1)
endif

time = dfloat(itime) / dfloat(rate)

if (irank==0) write(*,'(a,f10.4)') "        Time needed (s) :     ",time

end subroutine time_display
