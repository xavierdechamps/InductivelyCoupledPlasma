!##########################################################
! SUBROUTINE find_vector
!  Goal: find index 'element' inside array 'vector'
!##########################################################
SUBROUTINE find_vector(vector,lengvec,element,front1,front2)
    USE module_icp
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(IN) :: element,lengvec
    INTEGER(ki), DIMENSION(lengvec), INTENT(IN) :: vector
    INTEGER(ki), INTENT(OUT) :: front1,front2
    
    ! Local parameters
    INTEGER(ki) :: i, ind
    
    front1 = 0
    front2 = 0
    ind = 0

    loop : DO i=1,lengvec
      IF (vector(i) .EQ. element) THEN
        ind = ind + 1
        IF (ind .EQ. 1) THEN
           front1 = i
        ELSE
           front2 = i
           EXIT loop
        END IF
      END IF
    END DO loop

END SUBROUTINE find_vector

!##########################################################
! SUBROUTINE find_mat
!  Goal: find an index inside a 2D array
!##########################################################
SUBROUTINE find_mat(matrix,sizemat,sizemat2,element,ind,pos,lengthpos)
    USE module_icp
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(IN) :: element,sizemat,sizemat2,lengthpos
    INTEGER(ki), INTENT(IN) :: matrix(sizemat,sizemat2)
    INTEGER(ki), INTENT(OUT):: ind
    INTEGER(ki), INTENT(OUT):: pos(lengthpos)
    
    ! Local parameters
    INTEGER(ki) :: i,j
    
    ind = 0 ! number of times that element appears in matrix
    pos = 0

    DO i=1,sizemat
       DO j=1,sizemat2
          IF ( matrix(i,j) .eq. element ) THEN
             ind = ind + 1
             pos(ind) = i
          END IF
          IF (matrix(i,j).eq.0) EXIT
       END DO
    END DO

END SUBROUTINE find_mat
