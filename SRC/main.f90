PROGRAM main
    USE module_icp
    IMPLICIT NONE
    INCLUDE "mpif.h"

    ! Local parameters
    INTEGER(ki) :: ok,ierr
    CHARACTER(LEN=100)::param_file
        
    IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
        WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
        WRITE(*,*) "icp name_of_parameter_file"
        GOTO 200
    ENDIF 
    CALL get_command_argument(1,param_file)
        
    ! Get the flow parameters
    CALL read_parameters(param_file,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "Problem reading the parameters"
      GOTO 200
    ENDIF
    
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    ! ------- Données du problème ----------
    omega = 2.d00 * pi * frequency
    mu0   = 4e-7*pi
    eps   = 1.E-12
    
    ok = 1
    ! --------------------------------------
    CALL browse_gmsh(mesh_file,length_names,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,ok)
    IF (ok == 0) THEN
       WRITE(*,*) "The program hasn't started because of a problem during the browsing of the mesh"
       GOTO 200
    ENDIF
    
    ! Allocate the memory for the arrays
    CALL mem_allocate()
    
    CALL init()
    
    ! Read the mesh and the value of sigma in the torch
    CALL read_gmsh(mesh_file,length_names,node,elem,nbr_nodes_per_elem,front,stencilElem,&
&                  nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,0,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ENDIF
    
    write(*,*) "Building the CSR structure"
    call createcsr(stencil,elem(:,1:3),stiffness,ia,ja,nbrnodes,nbrelem,numnz) 
    
    write(*,*) "Building the stiffness matrix"
    call getStiffness()

   write(*,*) "Imposing the boundary conditions"
   call setBC()

   CALL solve(ok)
   
   if (ok.eq.1) THEN
      CALL postpro()
      CALL write_gmsh(ok)
   end if
   
200 continue
    if (irank.eq.0) write(*,*) "End of the simulation"
    
    ! Deallocate the memory for the arrays
    CALL mem_deallocate()
    ! call deallocate_end
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
  write(*,'(a,f10.4)') "        Time needed (s) :     ",time

end subroutine time_display
