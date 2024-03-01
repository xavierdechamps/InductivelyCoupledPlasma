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
    CALL read_gmsh(mesh_file,length_names,node,elem,nbr_nodes_per_elem,front,&
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

    INTEGER(ki) :: ierr, ielm, numdigits,k
    CHARACTER(LEN=2) :: numdig
    CHARACTER(LEN=9) :: formatreal
    REAL(kr) :: sigma_init(nbrElem),z,r,z1,z2,r1,r2,rad,rad2
    REAL(kr) :: z3,r3,z4,r4,slope
    LOGICAL  :: cond1,cond2,cond3,cond4,cond5,cond6
    
    CALL get_number_digits_integer(nbrNodes,numdigits)
    WRITE(numdig,'(A,I1)') 'I',numdigits+1

    formatreal = 'ES24.15E3'

    !************************************* SIGMA  
    sigma_init = zero
    DO ielm=1,nbrElem
    
      z=zero
      r=zero
      do k=1,nbr_nodes_per_elem(ielm)
        z=z+node(elem(ielm,k),1)
        r=r+node(elem(ielm,k),2)
      enddo
      z = z/nbr_nodes_per_elem(ielm)
      r = r/nbr_nodes_per_elem(ielm)
    
      IF (elem(ielm,5).eq.2) THEN ! Inside the torch
        
      ! First zone sigma = 49
        z1 = 0.056; r1 = zero; z2 = 0.0526; r2 = 0.0037;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond1 = r3.le.zero
        z1 = 0.0517; r1 = 0.0045; z2 = 0.0461; r2 = 0.01;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond2 = r3.le.zero
        z1 = 0.0461; r1 = 0.01; z2 = 0.0576; r2 = 0.0166
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond3 = r3.le.zero
        z1 = 0.0576; r1 = 0.0166; z2 = 0.1176; r2 = 0.0207
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond4= r3.le.zero
        z1 = 0.1176; r1 = 0.0207; z2 = 0.2; r2 = 0.0182
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond5= r3.le.zero
        z1 = 0.0526; r1 = 0.0037; z2 = 0.0517; r2 = 0.0045;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond6= r3.le.zero
        if (cond1 .and. cond2 .and. cond3 .and. cond4 .and. cond5 .and. cond6) sigma_init(ielm) = 49.d00
        
      ! Second zone sigma = 932
        z1 = 0.0566; r1 = zero; z2 = 0.0533; r2 = 0.0054;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond1 = r3.le.zero
        z1 = 0.0533; r1 = 0.0054; z2 = 0.0505; r2 = 0.01;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond2 = r3.le.zero
        z1 = 0.0505; r1 = 0.01; z2 = 0.0575; r2 = 0.0141
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond3 = r3.le.zero
        z1 = 0.0575; r1 = 0.0141; z2 = 0.1182; r2 = 0.019
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond4= r3.le.zero
        z1 = 0.1182; r1 = 0.019; z2 = 0.2; r2 = 0.0115
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond5= r3.le.zero
        if (cond1 .and. cond2 .and. cond3 .and. cond4 .and. cond5 ) sigma_init(ielm) = 932.d00
        
      ! Third zone sigma = 1836
        z1 = 0.057; r1 = zero; z2 = 0.0564; r2 = 0.004;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond1 = r3.le.zero
        z1 = 0.0564; r1 = 0.004; z2 = 0.0543; r2 = 0.0092;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond2 = r3.le.zero
        z1 = 0.0543; r1 = 0.0092; z2 = 0.06; r2 = 0.0125
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond3 = r3.le.zero
        z1 = 0.06; r1 = 0.0125; z2 = 0.120; r2 = 0.0175
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond4= r3.le.zero
        z1 = 0.12; r1 = 0.0175; z2 = 0.2; r2 = 0.0
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond5= r3.le.zero
        if (cond1 .and. cond2 .and. cond3 .and. cond4 .and. cond5 ) sigma_init(ielm) = 1836.d00
        
      ! Fourth zone sigma = 2358
        z1 = 0.0592; r1 = zero; z2 = 0.0576; r2 = 0.0087;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond1 = r3.le.zero
        z1 = 0.0576; r1 = 0.0087; z2 = 0.0854; r2 = 0.014;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond2 = r3.le.zero
        z1 = 0.0854; r1 = 0.014; z2 = 0.1186; r2 = 0.0163
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond3 = r3.le.zero
        z1 = 0.1186; r1 = 0.0163; z2 = 0.169; r2 = 0.0
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond4= r3.le.zero
        if (cond1 .and. cond2 .and. cond3 .and. cond4 ) sigma_init(ielm) = 2358.d00
        
      ! Fifth zone sigma = 2570
        z1 = 0.0679; r1 = 0.002; z2 = 0.0611; r2 = 0.0083;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond1 = r3.le.zero
        z1 = 0.0611; r1 = 0.0083; z2 = 0.086; r2 = 0.0135;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond2 = r3.le.zero
        z1 = 0.086; r1 = 0.0135; z2 = 0.117; r2 = 0.0153
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond3 = r3.le.zero
        z1 = 0.117; r1 = 0.0153; z2 = 0.1355; r2 = 0.0055
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond4= r3.le.zero
        z1 = 0.1355; r1 = 0.0055; z2 = 0.0679; r2 = 0.002
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond5= r3.le.zero
        if (cond1 .and. cond2 .and. cond3 .and. cond4 .and. cond5 ) sigma_init(ielm) = 2570.d00
        
      ! Sixth zone sigma = 2679
        z1 = 0.0695; r1 = 0.0071; z2 = 0.0855; r2 = 0.0133;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond1 = r3.le.zero
        z1 = 0.0855; r1 = 0.0133; z2 = 0.1175; r2 = 0.0141;
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond2 = r3.le.zero
        z1 = 0.1175; r1 = 0.0141; z2 = 0.125; r2 = 0.0089
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond3 = r3.le.zero
        z1 = 0.125; r1 = 0.0089; z2 = 0.0695; r2 = 0.0071
        r3 = (z2-z1)*(r-r1) - (r2-r1)*(z-z1)
        cond4= r3.le.zero
        if (cond1 .and. cond2 .and. cond3 .and. cond4 ) sigma_init(ielm) = 2679.d00
        
      ! Seventh zone sigma = 2756
        z1 = 0.0875; r1 = 0.0107; z2 = 0.0875; r2 = 0.0125
        rad = sqrt((z2-z1)**2+(r2-r1)**2)
        rad2 = sqrt((z-z1)**2+(r-r1)**2)
        if (rad2.le.rad) sigma_init(ielm) = 2756.d00
        
      ENDIF
    ENDDO
    
    CALL write_gmsh_initial_solution(sigma_init)

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