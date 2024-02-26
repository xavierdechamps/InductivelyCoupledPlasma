PROGRAM main
    USE module_icp
    IMPLICIT NONE
    INCLUDE "mpif.h"

    ! Local parameters
    INTEGER(ki) :: ok,ierr
        
    IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
        WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
        WRITE(*,*) "icp name_of_parameter_file"
        GOTO 200
    ENDIF 
    CALL get_command_argument(1,mesh_file)
    
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,irank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    ! ------- Données du problème ----------
    ! mesh_file     = "StructRafVertHor_fine.msh"
    file_gmsh    = "essai.msh"
    file_dat     = "essai.dat"
    
    Icoil = 1.0d00
!    sigma = 100.d00
    omega = 2.d00 * pi * 27.6e6
    mu0   = 4e-7*pi
    eps   = 1.E-12
    
    coils(1,1:2) = (/0.04d00  , 0.019d00/)
    coils(2,1:2) = (/0.053d00 , 0.019d00/)
    coils(3,1:2) = (/0.066d00 , 0.019d00/)
    
    ! coils(4,1:2) = (/0.027d00  , 0.019d00/)
    ! coils(5,1:2) = (/0.079d00  , 0.019d00/)
    
    ! coils(6,1:2) = (/0.014d00  , 0.019d00/)
    ! coils(7,1:2) = (/0.092d00  , 0.019d00/)
    
    ! coils(8,1:2) = (/0.001d00  , 0.019d00/)
    ! coils(9,1:2) = (/0.105d00  , 0.019d00/)
    
    ! coils(10,1:2) = (/0.0075d00  , 0.023d00/)
    ! coils(11,1:2) = (/0.0205d00  , 0.023d00/)
    ! coils(12,1:2) = (/0.0335d00  , 0.023d00/)
    ! coils(13,1:2) = (/0.0465d00  , 0.023d00/)
    ! coils(14,1:2) = (/0.0595d00  , 0.023d00/)
    ! coils(15,1:2) = (/0.0725d00  , 0.023d00/)
    ! coils(16,1:2) = (/0.0855d00  , 0.023d00/)
    ! coils(17,1:2) = (/0.0985d00  , 0.023d00/)
    
    ok = 1
    ! --------------------------------------
    ! if (irank.eq.0) then
       CALL browse_gmsh(mesh_file,length_names,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,ok)
       IF (ok == 0) THEN
          WRITE(*,*) "The program hasn't started because of a problem during the browsing of the mesh"
          GOTO 200
       ENDIF
    ! endif
    
    ! Allocate the memory for the arrays
    CALL mem_allocate(node,front,elem,nbr_nodes_per_elem,U0,BoundCond,&
&                     rhs,stencil,stiffness,ia,ja,mat,&
&                     nbrNodes,nbrElem,nbrFront,0)
    
    CALL init()
    
    ! Read the mesh and the initial solution / boundary conditions
    CALL read_gmsh(U0,nbvar*nbrElem,mesh_file,length_names,node,elem,nbr_nodes_per_elem,front,&
&                  BoundCond,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,0,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ENDIF
    
    ! CALL partitioner()
    if (irank.eq.0) write(*,*) "Constructing the CSR structure"
    call createcsr(stencil,elem(:,1:3),stiffness,ia,ja,nbrnodes,nbrelem,numnz) 
    
    if (irank.eq.0) write(*,*) "Constructing the stiffness matrix"
    call getStiffness()

   ! if (irank.eq.0) write(*,*) "Imposing the boundary conditions"
   call setBC()

   CALL solve(ok)
      ! call resolution(ok)
   
   CALL write_gmsh(ok)
   
   ! if (ok.eq.1 .and. irank.eq.0) call ecriture_gmsh
   

200 continue
    if (irank.eq.0) write(*,*) "End of the simulation"
    
    ! Deallocate the memory for the arrays
    CALL mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,BoundCond,&
&                       rhs,stencil,stiffness,ia,ja,mat)
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

if (irank==0) write(*,'(a,f10.4)') "        Time needed (s) :     ",time

end subroutine time_display
