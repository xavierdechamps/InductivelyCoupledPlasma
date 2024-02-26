!##########################################################
! SUBROUTINE read_parameters
!  Goal: read the parameters from the parameter file given 
!        by the user in the command window
!##########################################################
SUBROUTINE read_parameters(file_name,ok)
    USE module_icp
    IMPLICIT NONE
    
    ! Subroutine parameters
    CHARACTER(LEN=100),INTENT(IN)  :: file_name
    INTEGER(ki),       INTENT(OUT) :: ok
    
    ! Local parameters
    CHARACTER(LEN=200)::tmp 
    CHARACTER(LEN=256) :: my_iomsg
    INTEGER(ki) :: i, ierr,line

    ok = 1
    line = 1
    OPEN(unit=20,file=file_name,status="old",iostat=ierr,iomsg=my_iomsg,form='formatted')
    IF (ierr .NE. 0) THEN
        WRITE(*,*) my_iomsg
        ok = 0
        RETURN
    ENDIF
    
    READ(20,iostat=ierr,err=8,fmt='(a50,a100)') mesh_file,tmp
    line = line + 1
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)') nbrBC,tmp
    line = line + 1
    IF (nbrBC.le.0 .or. nbrBC.gt.5) THEN
      WRITE(*,*) "read_parameters: the number of boundary types is not in the range [1-5]"
      GOTO 50
    ENDIF
    CLTable = 0
    DO i=1,nbrBC
      READ(20,iostat=ierr,err=8,fmt='(I9,a100)') CLTable(1,i),tmp
      line = line + 1
      CLTable(2,i) = i-1
    ENDDO
    READ(20,iostat=ierr,err=8,fmt='(a50,a100)')       file_gmsh,tmp
    line = line + 1
    READ(20,iostat=ierr,err=8,fmt='(a50,a100)')       file_dat,tmp
    line = line + 1
    READ(20,iostat=ierr,err=8,fmt='(ES15.7E3,a100)')  Icoil,tmp
    line = line + 1
    READ(20,iostat=ierr,err=8,fmt='(ES15.7E3,a100)')  sigma,tmp
    line = line + 1
    READ(20,iostat=ierr,err=8,fmt='(ES15.7E3,a100)')  frequency,tmp
    line = line + 1
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)')        nbrCoils,tmp
    line = line + 1
    IF (nbrCoils.gt.0) ALLOCATE(coils(nbrCoils,2))
    DO i=1,nbrCoils
      READ(20,iostat=ierr,err=8,fmt='(2ES15.7E3,a100)') coils(i,1),coils(i,2),tmp
      line = line + 1
    ENDDO
        
    CLOSE (20)
! ! normal exit
    CALL print_parameters(file_name)
    RETURN

! ! exit on error during read operation 
 8   WRITE(*,*) "I/O error in parameters at line ", line
    
 50 CONTINUE   
    ok = 0
    call print_correct_parameters()

END SUBROUTINE read_parameters

!##########################################################
! SUBROUTINE print_correct_parameters
!  Goal: print on the screen the parameters that should 
!        appear in the parameter file
!##########################################################
SUBROUTINE print_correct_parameters()
    USE module_icp, ONLY : ki
    IMPLICIT NONE
    
    ! Local parameters
    INTEGER(ki) :: i
    
    WRITE(*,*) "The input file must contain the following lines"
    WRITE(*,*) ('*',i=1,79)
    WRITE(*,'(a12,18x,a)') "CHARACTER*50",": name of mesh file"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": number of boundary types as built in gmsh, see"
    WRITE(*,'(32x,a)')     "the following five physical tags"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": physical tag associated to the Far field"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": [optional] physical tag associated to the Axis"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": [optional] physical tag associated to the Torch wall"
    WRITE(*,'(a12,18x,a)') "CHARACTER*50",": name of the output gmsh file"
    WRITE(*,'(a12,18x,a)') "CHARACTER*50",": name of the output dat file"
    WRITE(*,'(a4,26x,a)')  "REAL",": Icoil, intensity of the electric current in the coil"
    WRITE(*,'(a4,26x,a)')  "REAL",": sigma, electric conductivity of the plasma"
    WRITE(*,'(a4,26x,a)')  "REAL",": excitation frequency [Herz] of the coil"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": number of coils"
    WRITE(*,'(a9,21x,a)')  "REAL REAL",": axial and radial coordinates of each coil"
    
    WRITE(*,*) ('*',i=1,79)
    
END SUBROUTINE print_correct_parameters

!##########################################################
! SUBROUTINE print_parameters
!  Goal: print on the screen the parameters that were read
!        in the parameter file
!##########################################################
SUBROUTINE print_parameters(file_name)
    USE module_icp
    IMPLICIT NONE
    
    ! Subroutine parameters
    CHARACTER(LEN=100),INTENT(IN)  :: file_name
    
    ! Local parameters
    INTEGER(ki) :: i
    
    WRITE(*,*) "The following parameters were read from parameter file '",trim(file_name),"'"
    WRITE(*,'(a50,a)')          mesh_file,": name of mesh file"
    WRITE(*,'(i9,41x,a)')       nbrBC,": number of boundary types as built in gmsh, see"
    DO i=1,nbrBC
      WRITE(*,'(i9,41x,a)')     CLTable(1,i),": boundary physical tag"
    ENDDO
    WRITE(*,'(a50,a)')          file_gmsh,": name of the output gmsh file"
    WRITE(*,'(a50,a)')          file_dat,": name of the output dat file"
    WRITE(*,'(ES15.7E3,35x,a)') Icoil,": Icoil, intensity of the electric current in the coil"
    WRITE(*,'(ES15.7E3,35x,a)') sigma,": sigma, electric conductivity of the plasma"
    WRITE(*,'(ES15.7E3,35x,a)') frequency,": excitation frequency [Herz] of the coil"
    WRITE(*,'(i9,41x,a)')       nbrCoils,": number of coils"
    DO i=1,nbrCoils
      WRITE(*,'(2ES15.7E3,20x,a)')     coils(i,1),coils(i,2),": axial and radial coordinates of each coil"
    ENDDO
    
END SUBROUTINE print_parameters
