!##########################################################
! SUBROUTINE mem_allocate
! 
!##########################################################
SUBROUTINE mem_allocate()
    USE module_icp
    IMPLICIT NONE
    
    WRITE(*,*) "Allocating the memory..."
        
    ! Requested in read_gmsh
    ALLOCATE(node(1:nbrNodes,1:2))
    ALLOCATE(front(1:nbrFront,1:4))
    ALLOCATE(elem(1:nbrElem,1:5))
    ALLOCATE(nbr_nodes_per_elem(1:nbrElem))
    ALLOCATE(U0(1:nbvar*nbrElem))
    ALLOCATE(Ecoils(1:nbrNodes))
    ALLOCATE(sigma_in(1:nbrElem))
    
    ! ALLOCATE(BoundCond(1:nbrFront,1:3))
    
    ! Solution and CSR format data
    ALLOCATE (rhs(1:nbvar*nbrNodes))
    ALLOCATE (stencil(1:nbrNodes,1:maxstencil))
    ALLOCATE (stiffness(1:4,1:4))
    ALLOCATE (ia(1:nbrNodes*nbvar+1))
    ALLOCATE (ja(1:nbrNodes*maxstencil*nbvar*nbvar))
    ALLOCATE (mat(1:nbrNodes*maxstencil*nbvar*nbvar))
    
END SUBROUTINE mem_allocate

!##########################################################
! SUBROUTINE mem_deallocate
! 
!##########################################################
SUBROUTINE mem_deallocate()
    USE module_icp
    IMPLICIT NONE
            
    WRITE(*,*) "Deallocating the memory..."
    
    IF (ALLOCATED(node))  DEALLOCATE(node)
    IF (ALLOCATED(front)) DEALLOCATE(front)
    IF (ALLOCATED(elem))  DEALLOCATE(elem)
    IF (ALLOCATED(nbr_nodes_per_elem))  DEALLOCATE(nbr_nodes_per_elem)
    IF (ALLOCATED(U0))    DEALLOCATE(U0)
    IF (ALLOCATED(Ecoils))    DEALLOCATE(Ecoils)
    IF (ALLOCATED(sigma_in))    DEALLOCATE(sigma_in)
    
    ! IF (ALLOCATED(BoundCond)) DEALLOCATE(BoundCond)
    IF (ALLOCATED(coils)) DEALLOCATE(coils)
    
    IF (ALLOCATED(stencil)) DEALLOCATE(stencil)
    IF (ALLOCATED(rhs)) DEALLOCATE(rhs)
    IF (ALLOCATED(ia)) DEALLOCATE(ia)
    IF (ALLOCATED(ja)) DEALLOCATE(ja)
    IF (ALLOCATED(mat)) DEALLOCATE(mat)
      
      ! if (allocated(stiffness)) deallocate(stiffness)
      ! if (allocated(nodCL)) deallocate(nodCL)
      ! if (allocated(loc2glob)) deallocate(loc2glob)
    
END SUBROUTINE mem_deallocate