
!##########################################################
! SUBROUTINE write_solution
!  Goal: write the restart solution in the dat file
!##########################################################
SUBROUTINE write_solution(U0,lengU0,file_dat,lengch,ok)
    USE module_icp, only : kr,ki
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT) :: ok
    INTEGER(ki), INTENT(IN)  :: lengU0,lengch
    REAL(kr), INTENT(IN)     :: U0(lengU0)
    CHARACTER(LEN=lengch), INTENT(IN) :: file_dat
    
    ! Local parameters
    INTEGER(ki) :: i, ierr

    ok = 1
    OPEN(UNIT=20,FILE=file_dat,STATUS="replace",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error writing the solution in .dat format"
      ok = 0
    END IF

    WRITE(20,'(ES15.7E3)') (U0(i),i=1,lengU0)

    CLOSE(UNIT=20)

END SUBROUTINE write_solution

!##########################################################
! SUBROUTINE read_solution
!  Goal: read the restart solution in the dat file
!##########################################################
SUBROUTINE read_solution(U0,lengU0,file_restart,lengch,ok)
    USE module_icp, only : kr,ki,zero
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT):: ok
    INTEGER(ki), INTENT(IN) :: lengU0,lengch
    REAL(kr), INTENT(OUT)   :: U0(lengU0)
    CHARACTER(LEN=lengch), INTENT(IN) :: file_restart
    
    ! Local parameters
    INTEGER(ki) :: i, ierr, count
    REAL(kr) :: a

    ok = 1
    U0 = zero
    OPEN(UNIT=20,FILE=file_restart,STATUS="old",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error occuring : file '",trim(file_restart),"' not found "
      ok = 0
      GOTO 200
    END IF
    count = 0

    READ(20,IOSTAT=ierr,FMT='(ES15.7E3)') a
    DO 
       count = count + 1
       READ(UNIT=20,IOSTAT=ierr,FMT='(ES15.7E3)') a
       IF (ierr .LT. 0) EXIT
    END DO

    IF (count .EQ. lengU0) THEN
       REWIND(20)
       READ(20,*) (U0(i), i=1,count)
    ELSE
       ok = 0
       WRITE(*,*) "Error occuring during the reading of the former solution : wrong vector length"
       WRITE(*,*) "Length of the data in the restart file: ",count
       WRITE(*,*) "Expected length of the data in the restart file: ",lengU0
    END IF
    CLOSE(UNIT=20)
    
200 CONTINUE  
END SUBROUTINE read_solution

!##########################################################
!SUBROUTINE write_gmsh
!  Goal: write the solution in the Gmsh .msh format
!##########################################################
! SUBROUTINE write_gmsh(U0,lengU0,file_gmsh,lengch,node,elem,front,nbrNodes,nbrElem,nbrTris,nbrQuads,&
! &                     nbrFront,nbr_nodes_per_elem,k,count)
SUBROUTINE write_gmsh()
    USE module_icp
    IMPLICIT NONE
    
    ! Local parameters
    INTEGER(ki) :: ierr, i, fin, numdigits
    CHARACTER(LEN=2) :: numdig
    CHARACTER(LEN=9) :: formatreal
    
    fin = 0
    formatreal = 'ES24.15E3'
    
    CALL get_number_digits_integer(nbrNodes,numdigits)
    WRITE(numdig,'(A,I1)') 'I',numdigits+1
        
    OPEN(UNIT=10,FILE=file_gmsh,STATUS="replace",IOSTAT=ierr,FORM='formatted')
    WRITE(10,'(T1,A11)') "$MeshFormat"
    WRITE(10,'(T1,A7)') "2.2 0 8"
    WRITE(10,'(T1,A14)') "$EndMeshFormat"
    WRITE(10,'(T1,A6)') "$Nodes"
    
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//',2'//formatreal//',I2)') (i, node(i,:),0,i=1,nbrNodes)
    WRITE(10,'(T1,A9)') "$EndNodes"
    WRITE(10,'(T1,A9)') "$Elements"
    WRITE(10,'(T1,'//numdig//')') nbrElem+nbrFront
    WRITE(10,'(T1,'//numdig//',2I2,'//numdig//',I2,2'//numdig//')') (i,1,2,front(i,3),1,front(i,1:2),i=1,nbrFront)    
    DO i=1,nbrElem
      IF (nbr_nodes_per_elem(i) .EQ. 3) THEN 
        write(10,'(T1,'//numdig//',2I2,'//numdig//',I2,3'//numdig//')') i+nbrFront,2,2,elem(i,5),1,elem(i,1:3)
      ELSE IF (nbr_nodes_per_elem(i) .EQ. 4) THEN 
        write(10,'(T1,'//numdig//',2I2,'//numdig//',I2,4'//numdig//')') i+nbrFront,3,2,elem(i,5),1,elem(i,1:4)
      END IF 
    ENDDO  
    WRITE(10,'(T1,A12)') "$EndElements"
    
    !*************************************
    ! Write the electric field due to the coils
    ! WRITE(10,'(T1,A9)') "$NodeData"
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,A9)') '"E_coils"'
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "3"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') nbrNodes
    ! WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, Ecoils(i),i=1,nbrNodes)
    ! WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the induced electric field
    WRITE(10,'(T1,A9)') "$NodeData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A16)') '"E_induced_real"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, U0(2*i-1),i=1,nbrNodes)
    WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the induced electric field
    ! WRITE(10,'(T1,A9)') "$NodeData"
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,A16)') '"E_induced_imag"'
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "3"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') nbrNodes
    ! WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, U0(2*i),i=1,nbrNodes)
    ! WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the total electric field
    WRITE(10,'(T1,A9)') "$NodeData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A12)') '"E_tot_imag"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, (U0(2*i)+Ecoils(i)),i=1,nbrNodes)
    WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the norm of the total electric field
    ! WRITE(10,'(T1,A9)') "$NodeData"
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,A12)') '"E_tot_norm"'
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "3"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') nbrNodes
    ! WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, sqrt(U0(2*i-1)**2 + (U0(2*i)+Ecoils(i))**2),i=1,nbrNodes)
    ! WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the phase of the total electric field
    ! WRITE(10,'(T1,A9)') "$NodeData"
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,A13)') '"E_tot_phase"'
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "3"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') nbrNodes
    ! WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, 180.0d00/pi*atan2(U0(2*i)+Ecoils(i),U0(2*i-1)),i=1,nbrNodes)
    ! WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the Joule heating
    ! WRITE(10,'(T1,A9)') "$NodeData"
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,A16)') '"Joule heating"'
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "3"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') nbrNodes
    ! WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, dissJoule(i),i=1,nbrNodes)
    ! WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the Lorentz force
    WRITE(10,'(T1,A9)') "$NodeData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A21)') '"Lorentz force axial"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, Fax(i),i=1,nbrNodes)
    WRITE(10,'(T1,A12)') "$EndNodeData"
    
    WRITE(10,'(T1,A9)') "$NodeData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A22)') '"Lorentz force radial"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, Frad(i),i=1,nbrNodes)
    WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the magnetic field, axial component, real part
    WRITE(10,'(T1,A9)') "$NodeData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A16)') '"Magnetic_axial"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,'//numdig//')') 0   ! 1st time step = real part
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, Bzr(i),i=1,nbrNodes)
    WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the magnetic field, axial component, imaginary part
    ! WRITE(10,'(T1,A9)') "$NodeData"
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,A16)') '"Magnetic_axial"'
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "3"
    ! WRITE(10,'(T1,'//numdig//')') 1   ! 2nd time step = imaginary part
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') nbrNodes
    ! WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, Bzi(i),i=1,nbrNodes)
    ! WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the magnetic field, radial component, real part
    WRITE(10,'(T1,A9)') "$NodeData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A17)') '"Magnetic_radial"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') 0
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,'//numdig//')') 0   ! 1st time step = real part
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, Brr(i),i=1,nbrNodes)
    WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    ! Write the magnetic field, radial component, imaginary part
    ! WRITE(10,'(T1,A9)') "$NodeData"
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,A17)') '"Magnetic_radial"'
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') 0
    ! WRITE(10,'(T1,A1)') "3"
    ! WRITE(10,'(T1,'//numdig//')') 1   ! 2nd time step = imaginary part
    ! WRITE(10,'(T1,A1)') "1"
    ! WRITE(10,'(T1,'//numdig//')') nbrNodes
    ! WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, Bri(i),i=1,nbrNodes)
    ! WRITE(10,'(T1,A12)') "$EndNodeData"
    
    !*************************************
    CLOSE(UNIT=10)

END SUBROUTINE write_gmsh

!##########################################################
! SUBROUTINE read_gmsh
!  Goal: read the mesh and the initial solution in the Gmsh .msh format
!##########################################################
SUBROUTINE read_gmsh(mesh_file,lengch,node,elem,nbr_nodes_per_elem,front,stencilElem,&
&                    nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,skip_data,ok)
    USE module_icp, only : kr,ki,nbvar,maxstencil,zero,time_begin, time_end,sigma_in
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT) :: ok
    INTEGER(ki), INTENT(IN)  :: nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront
    INTEGER(ki), INTENT(IN)  :: lengch,skip_data
    REAL(kr), INTENT(OUT)    :: node(nbrNodes,2)
    INTEGER(ki), INTENT(OUT)    :: elem(nbrElem,5)
    INTEGER(ki), INTENT(OUT)    :: front(nbrFront,4)
    INTEGER(ki), INTENT(OUT)    :: stencilElem(nbrNodes,maxstencil+1)
    INTEGER(ki), INTENT(OUT)    :: nbr_nodes_per_elem(nbrElem)
    CHARACTER(LEN=lengch), INTENT(IN) :: mesh_file

    ! Local parameters
    CHARACTER(len=256)    :: line,dataname
    INTEGER(ki) :: i, j, ierr, a, b, c, d, e, nbrElemTot,nbrNodeEdge
    INTEGER(ki) :: istep,entnump,entnumc,entnums
    INTEGER(ki) :: nodesnument,surfnument,tag
    INTEGER(ki) :: ind,indf,inde
    REAL(kr) :: header,tmp
    EXTERNAL :: time_display
    INTEGER(ki), ALLOCATABLE :: entities_curves(:,:)
    INTEGER(ki), ALLOCATABLE :: entities_surfs(:,:)
    
    ind = 0
    indf = 0
    inde = 0
    nbrElemTot = 0
    nbrNodeEdge = 0
    
    CALL sampletime(time_begin)
    ok = 1
    WRITE(*,*) "-------------------------------------------------------"
    WRITE(*,*) "Name of the mesh data : ", TRIM(ADJUSTL(mesh_file))
    
    istep = 1
    
    OPEN(UNIT=10,FILE=mesh_file,STATUS="old",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error occuring : file '",TRIM(mesh_file),"' not found "
      ok = 0
      GOTO 100
    END IF

    READ(10,*) line, header, b, c
    
    IF (header .LT. 2.) THEN
      WRITE(*,*) "Error : wrong format for gmsh data"
      ok = 0
      GOTO 100
    END IF
    
    IF (header .GE. 2. .AND. header .LT. 4.0) THEN ! Format 2.1
      GOTO 10
    ELSE IF (header .GE. 4.0) THEN
      GOTO 20
    END IF
    
    istep = 2
    
    node = zero
    
10  DO WHILE (line .NE. "$Nodes")
      READ(10,*) line
    END DO

    READ(10,*) a

    DO i=1,nbrNodes ! Read the node coordinates
      READ(10,*) a, node(i,1), node(i,2), tmp
    END DO
    WRITE(*,101) nbrNodes
    
    istep = 3

    DO WHILE (line .NE. "$Elements")
      READ(10,*) line
    END DO

    READ(10,*) a
    
    ! front (:,1:2) : node IDs composing the boundary edge
    ! front (:,3)   : physical tag of the boundary edge
    ! front (:,4)   : ID of the 2D element linked to the boundary edge
    
    ! elem(:,1:4)   : node IDs composing the 2D element
    ! elem(:,5)     : physical tag of the 2D element
    
    ! stencilElem(:,1)     : number of 2D elements around each node
    ! stencilElem(:,2:end) : IDs of the 2D elements around each node
    nbr_nodes_per_elem = 0
    elem               = 0
    front              = 0
    stencilElem        = 0
    
    READ(10,*) (a,a,a,front(i,3),a,front(i,1),front(i,2),i=1,nbrFront)
    
    DO i=1,nbrElem
      READ(10,*) a,inde
      IF ( inde .EQ. 2 ) THEN
        nbr_nodes_per_elem(i) = 3
        BACKSPACE 10
        READ(10,*) a,a,a,elem(i,5),a,elem(i,1),elem(i,2),elem(i,3)
      ELSE IF (inde .EQ. 3) THEN
        nbr_nodes_per_elem(i) = 4
        BACKSPACE 10
        READ(10,*) a,a,a,elem(i,5),a,elem(i,1),elem(i,2),elem(i,3),elem(i,4)
      END IF
      
      DO j=1,nbr_nodes_per_elem(i)
        ! Increase counter of number of elements around each node
        a = stencilElem(elem(i,j),1) + 1
        stencilElem(elem(i,j),1) = a
        ! Set element ID around each node
        stencilElem(elem(i,j),a+1) = i
      ENDDO
            
    END DO
    
    GOTO 50
    
20  DO WHILE (line .NE. "$Entities")
      READ(10,*) line
    END DO
    
    READ(10,*) entnump,entnumc,entnums
    
    ALLOCATE( entities_curves( 1:entnumc , 2 ) )
    ALLOCATE( entities_surfs( 1:entnums , 2 ) )
    
    DO i=1,entnump
     READ(10,*)
    ENDDO
    DO i=1,entnumc
      READ(10,*) entities_curves(i,1),tmp,tmp,tmp,tmp,tmp,tmp,a,entities_curves(i,2)
    ENDDO
    DO i=1,entnums
      READ(10,*) entities_surfs(i,1),tmp,tmp,tmp,tmp,tmp,tmp,a,entities_surfs(i,2)
    ENDDO
    
    DO WHILE (line .NE. "$Nodes")
      READ(10,*) line
    END DO
    READ(10,*) nodesnument
    
    DO i=1,nodesnument
    ! Loop on the number of entities inside the group Nodes
      READ(10,*) a,a,a,b
      DO c=1,b
        READ(10,*) 
      END DO
      DO c=1,b
        READ(10,*) node(ind+c,1), node(ind+c,2)
      END DO
      ind = ind + b
    END DO
    WRITE(*,101) nbrNodes
    
    DO WHILE (line .NE. "$Elements")
      READ(10,*) line
    END DO
    READ(10,*) surfnument
      
    DO i=1,surfnument
    ! Loop on the number of entities inside the group Elements
      READ(10,*) a,tag,ind,b
      DO c=1,b
        IF (ind .EQ. 1) THEN ! Front
          indf = indf + 1
          READ(10,*) a,front(indf,1),front(indf,2)
          CALL find_vector(entities_curves(:,1),entnumc,tag,d,e)
          front(indf,3) = entities_curves(d,2)
          
        ELSE IF (ind .EQ. 2) THEN ! Triangle
          inde = inde + 1
          nbr_nodes_per_elem(inde) = 3
          READ(10,*) a,elem(inde,1),elem(inde,2),elem(inde,3)
          CALL find_vector(entities_surfs(:,1),entnums,tag,d,e)
          elem(inde,5) = entities_surfs(d,2)
                    
        ELSE IF (ind .EQ. 3) THEN ! Quadrangle
          inde = inde + 1
          nbr_nodes_per_elem(inde) = 4
          READ(10,*) a,elem(inde,1),elem(inde,2),elem(inde,3),elem(inde,4)
          CALL find_vector(entities_surfs(:,1),entnums,tag,d,e)
          elem(inde,5) = entities_surfs(d,2)
          
        END IF
        
      END DO
    END DO
    
    DEALLOCATE( entities_curves , entities_surfs )
    
    GOTO 50
    
50  WRITE(*,102) nbrFront
    WRITE(*,103) nbrTris
    WRITE(*,104) nbrQuads
    
    IF (skip_data.EQ.1) GOTO 100
    
    dataname = "Sigma"
    DO WHILE (line .NE. "$NodeData")
      READ(10,*,END=90) line
    END DO
    DO i=1,2
      READ(10,*,END=90) line
    ENDDO
!    IF (TRIM(line) .NE. '"Sigma"') GO TO 90
    DO i=1,5
      READ(10,*,END=90) line
    ENDDO
    READ(10,*,END=90) a    
    IF (a.ne.nbrNodes) THEN
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      WRITE(*,*) "Error reading the initial height."
      WRITE(*,*) "The number of elements with initial height is",a
      WRITE(*,*) "The number of 2D elements in the mesh is",nbrElem
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      ok = 0
      GOTO 100
    ENDIF
    READ(10,*,END=90) (a,sigma_in(i),i=1,nbrNodes)
    GOTO 100
    
 90 ok = 0
    WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    WRITE(*,*) "Error reading the following data: ", TRIM(dataname)
    WRITE(*,*) "Read the following data in the file: ", TRIM(line)
    WRITE(*,*) "Check that the name of the mesh file points to the correct file"
    WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    
100 CALL sampletime(time_end)
    CALL time_display
    WRITE(*,*) "-------------------------------------------------------"
    101 FORMAT(" Number of nodes : ",T41,I6)
    102 FORMAT(" Number of edge elements :", T41, I6)
    103 FORMAT(" Number of triangular surface elements :",T41,I6)
    104 FORMAT(" Number of rectangular surface elements :",T41,I6)
    CLOSE(UNIT=10)
    
END SUBROUTINE read_gmsh

!##########################################################
! SUBROUTINE write_gmsh_initial_solution
!  Goal: write the initial solution in the Gmsh .msh format
!##########################################################
SUBROUTINE write_gmsh_initial_solution()
    USE module_icp
    IMPLICIT NONE
    
    ! Local parameters
    INTEGER(ki) :: ierr, i, numdigits
    CHARACTER(LEN=2) :: numdig
    CHARACTER(LEN=9) :: formatreal
    
    CALL get_number_digits_integer(nbrNodes,numdigits)
    WRITE(numdig,'(A,I1)') 'I',numdigits+1

    formatreal = 'ES24.15E3'
    
    OPEN(UNIT=10,FILE=file_gmsh,STATUS="replace",IOSTAT=ierr,FORM='formatted')
    WRITE(10,'(T1,A11)') "$MeshFormat"
    WRITE(10,'(T1,A7)') "2.2 0 8"
    WRITE(10,'(T1,A14)') "$EndMeshFormat"
    WRITE(10,'(T1,A6)') "$Nodes"
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//',2'//formatreal//',I2)') (i, node(i,:),0,i=1,nbrNodes)
    WRITE(10,'(T1,A9)') "$EndNodes"
    WRITE(10,'(T1,A9)') "$Elements"
    WRITE(10,'(T1,'//numdig//')') nbrElem+nbrFront
    WRITE(10,'(T1,'//numdig//',2I2,'//numdig//',I2,2'//numdig//')') (i,1,2,front(i,3),1,front(i,1:2),i=1,nbrFront)    
    DO i=1,nbrElem
      ! IF (nbr_nodes_per_elem(i) .EQ. 3) THEN 
        WRITE(10,'(T1,'//numdig//',2I2,'//numdig//',I2,3'//numdig//')') i+nbrFront,2,2,elem(i,5),1,elem(i,1:3)
      ! ELSE IF (nbr_nodes_per_elem(i) .EQ. 4) THEN 
        ! WRITE(10,'(T1,'//numdig//',2I2,'//numdig//',I2,4'//numdig//')') i+nbrFront,3,2,elem(i,5),1,&
! &                                                        elem(i,1:4)
      ! END IF 
    ENDDO    
    WRITE(10,'(T1,A12)') "$EndElements"
    
    !************************************* ELECTRICAL CONDUCTIVITY
    WRITE(10,'(T1,A9)') "$NodeData"
    WRITE(10,'(T1,I1)') 1
    WRITE(10,'(T1,A7)') '"Sigma"'
    WRITE(10,'(T1,I1)') 1
    WRITE(10,'(T1,I1)') 0
    WRITE(10,'(T1,I1)') 3
    WRITE(10,'(T1,I1)') 0
    WRITE(10,'(T1,I1)') 1
    WRITE(10,'(T1,'//numdig//')') nbrNodes
    WRITE(10,'(T1,'//numdig//','//formatreal//')') (i, sigma_in(i),i=1,nbrNodes)
    WRITE(10,'(T1,A12)') "$EndNodeData"
    
    CLOSE(UNIT=10)
    
END SUBROUTINE write_gmsh_initial_solution

!##########################################################
! SUBROUTINE read_gmsh
!  Goal: read the mesh and the initial solution in the Gmsh .msh format
!##########################################################
SUBROUTINE browse_gmsh(mesh_file,lengch,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,ok)
    USE module_icp, only : ki,kr
    IMPLICIT NONE
    
    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT) :: ok
    INTEGER(ki), INTENT(OUT) :: nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront
    INTEGER(ki), INTENT(IN)  :: lengch
    CHARACTER(LEN=lengch), INTENT(IN) :: mesh_file

    ! Local parameters
    CHARACTER(len=256)    :: line
    INTEGER(ki) :: i, j, ierr, a, b, c
    INTEGER(ki) :: nbrEntities,nbrElemTot,nbrNodeEdge
    REAL(kr) :: header,tmp

    ok = 1
    nbrNodes = 0
    nbrElem = 0
    nbrTris= 0
    nbrQuads = 0
    nbrFront = 0
    nbrElemTot = 0
    nbrNodeEdge = 0
    
    OPEN(UNIT=10,FILE=mesh_file,STATUS="old",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error occuring : file '",TRIM(mesh_file),"' not found "
      ok = 0
      GOTO 200
    END IF

    READ(10,*) line, header, b, c
    
    IF (header .LT. 2.) THEN
      WRITE(*,*) "Error : wrong format for gmsh data"
      ok = 0
      GOTO 200
    END IF
    
    IF (header .GE. 2. .AND. header .LT. 4.0) THEN ! Format 2.1
      GOTO 110
    ELSE IF (header .GE. 4.0) THEN
      GOTO 120
    END IF
    
110 DO WHILE (line .NE. "$Nodes")
      READ(10,*) line
    END DO

    READ(10,*) nbrNodes
    
    DO WHILE (line .NE. "$Elements")
      READ(10,*) line
    END DO

    READ(10,*) nbrElemTot
    
    DO j=1,nbrElemTot
      READ(10,*) a,i
      IF (i .EQ. 1) THEN
        nbrFront = nbrFront + 1
      ELSE IF (i .EQ. 2) THEN
        nbrTris = nbrTris + 1
      ELSE IF (i .EQ. 3) THEN
        nbrQuads = nbrQuads + 1
      ELSE IF (i .EQ. 15) THEN
        nbrNodeEdge = nbrNodeEdge + 1
      END IF
    END DO
    GOTO 200
    
120 DO WHILE (line .NE. "$Nodes")
      READ(10,*) line
    END DO

    READ(10,*) a,nbrNodes
    
    DO WHILE (line .NE. "$Elements")
      READ(10,*) line
    END DO

    READ(10,*) nbrEntities,nbrElemTot
    
    DO j=1,nbrEntities
      READ(10,*) a,a,i,b
      
      IF (i .EQ. 1) THEN
        nbrFront = nbrFront + b
      ELSE IF (i .EQ. 2) THEN
        nbrTris = nbrTris + b
      ELSE IF (i .EQ. 3) THEN
        nbrQuads = nbrQuads + b
      ELSE IF (i .EQ. 15) THEN
        nbrNodeEdge = nbrNodeEdge + b
      END IF
      
      DO c=1,b
         READ(10,*) a
      ENDDO
      
    ENDDO
    
    GOTO 200

200 nbrElem = nbrTris + nbrQuads
    CLOSE(UNIT=10)
    
END SUBROUTINE browse_gmsh

!##########################################################
! SUBROUTINE get_number_digits_integer
!  Goal: Get the number of digits in an Integer number
!##########################################################
SUBROUTINE get_number_digits_integer(input,num_digits)
    USE module_icp, only : ki
    IMPLICIT NONE
    
    ! Subroutine parameters
    INTEGER(ki), INTENT(IN)  :: input
    INTEGER(ki), INTENT(OUT) :: num_digits

    ! Local parameters
    INTEGER(ki) :: modul,div
    
    num_digits = 0
    div = 1
    modul = input / div
    
    DO WHILE (modul.NE.0)
      div = div * 10
      modul = input / div
      num_digits = num_digits + 1      
    END DO
    
END SUBROUTINE get_number_digits_integer