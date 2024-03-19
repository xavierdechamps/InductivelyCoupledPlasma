PROGRAM BUILD_INITIAL_SOLUTION
    USE module_icp
    IMPLICIT NONE
    
    INTEGER(ki) :: ok
    CHARACTER(LEN=100)::param_file
    
    ! Get the number of arguments
    IF(COMMAND_ARGUMENT_COUNT().NE.2)THEN
      WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
      WRITE(*,*) "build_initial_condition   mesh_file.msh   initial_condition_file.msh"
      GOTO 200
    ENDIF 
    CALL get_command_argument(1,mesh_file)
    CALL get_command_argument(2,file_gmsh)
        
    ! Browse the mesh to get the size of the arrays
    CALL browse_gmsh(mesh_file,length_names,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the browsing of the mesh"
      GOTO 200
    endif
    
    ! Allocate the memory for the arrays
    CALL mem_allocate()

    ! Read the mesh and the initial solution / boundary conditions
    CALL read_gmsh(mesh_file,length_names,node,elem,nbr_nodes_per_elem,front,stencilElem,&
&                  nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,1,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ELSE
      CALL write_initial_condition_gmsh()
    ENDIF
    
    ! Deallocate the memory for the arrays
    CALL mem_deallocate()
    
200 CONTINUE
    WRITE(*,*) "End of the program"
    
END PROGRAM

! ******************************************************************************
SUBROUTINE write_initial_condition_gmsh()
    USE module_icp
    IMPLICIT NONE

    INTEGER(ki) :: ierr, inod, numdigits,k
    CHARACTER(LEN=2) :: numdig
    CHARACTER(LEN=9) :: formatreal
    REAL(kr) :: z,r,z1,z2,r1,r2,rad,rad2
    REAL(kr) :: z3,r3,z4,r4,slope
    LOGICAL  :: cond1,cond2,cond3,cond4
    
    CALL get_number_digits_integer(nbrNodes,numdigits)
    WRITE(numdig,'(A,I1)') 'I',numdigits+1

    formatreal = 'ES24.15E3'

    !************************************* SIGMA  
    sigma_in = zero
    DO inod=1,nbrNodes
      IF (elem(ielm,5).eq.2) THEN ! Inside the torch
        sigma_in(inod) = 1000.0d00
      ENDIF
    ENDDO
    CALL write_gmsh_initial_solution(sigma_in)

END SUBROUTINE write_initial_condition_gmsh

! *******************************************************************
SUBROUTINE sampletime(counter)
    USE module_icp, ONLY : kr,ki
    IMPLICIT NONE
          
    ! variables passed through header
    INTEGER(ki) ::counter
      
    ! variables declared locally
    INTEGER(ki) ::rate, contmax
 
    ! Determine CPU time
    CALL SYSTEM_CLOCK(counter, rate, contmax )
!-----------------------------------------------------------------------
END SUBROUTINE sampletime
!-----------------------------------------------------------------------

! *******************************************************************
SUBROUTINE time_display
    USE module_icp
    IMPLICIT NONE

    REAL(kr) :: time
    INTEGER(ki) :: job, cont3
    INTEGER(ki) :: rate, contmax, itime

    CALL SYSTEM_CLOCK(cont3, rate, contmax )
    IF (time2 .ge. time1) THEN
       itime=time2-time1
    ELSE
       itime=(contmax - time1) + (time2 + 1)
    ENDIF

    time = DFLOAT(itime) / DFLOAT(rate)

    WRITE(*,'(a,f10.4)') " Time needed (s) :            ",time
    
END SUBROUTINE time_display