!##########################################################
! SUBROUTINE mem_allocate
! 
!##########################################################
SUBROUTINE mem_allocate(node,front,elem,nbr_nodes_per_elem,U0,BoundCond,&
&                       nbrNodes,nbrElem,nbrFront,only_mesh)
    USE module_icp, only : kr,ki,nbvar
    IMPLICIT NONE

!   WARNING !!!!!!!!!!!!
! If you change the input/output parameters of this routine
! you also have to change them in the interface 
!      module_shallow.f90 / module_mem_allocate

    REAL(kr), ALLOCATABLE    :: U0(:)
    REAL(kr), ALLOCATABLE    :: node(:,:)
    REAL(kr), ALLOCATABLE    :: BoundCond(:,:)
    INTEGER(ki), ALLOCATABLE    :: elem(:,:)
    INTEGER(ki), ALLOCATABLE    :: front(:,:)
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elem(:)
    
    INTEGER(ki), INTENT(IN) :: nbrNodes,nbrElem,nbrFront,only_mesh
    
    WRITE(*,*) "Allocating the memory..."
    
    ! Requested in read_gmsh
    ALLOCATE(node(1:nbrNodes,1:2))
    ALLOCATE(front(1:nbrFront,1:4))
    ALLOCATE(elem(1:nbrElem,1:5))
    ALLOCATE(nbr_nodes_per_elem(1:nbrElem))
    ALLOCATE(U0(1:nbvar*nbrElem))
    ALLOCATE(BoundCond(1:nbrFront,1:3))
    
END SUBROUTINE mem_allocate

!##########################################################
! SUBROUTINE mem_deallocate
! 
!##########################################################
SUBROUTINE mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,BoundCond)
    USE module_icp, only : kr,ki
    IMPLICIT NONE

!   WARNING !!!!!!!!!!!!
! If you change the input/output parameters of this routine
! you also have to change them in the interface 
!      module_shallow.f90 / module_mem_allocate

    REAL(kr), ALLOCATABLE    :: U0(:)
    REAL(kr), ALLOCATABLE    :: node(:,:)
    REAL(kr), ALLOCATABLE    :: BoundCond(:,:)
    INTEGER(ki), ALLOCATABLE    :: elem(:,:)
    INTEGER(ki), ALLOCATABLE    :: front(:,:)
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elem(:)
            
    WRITE(*,*) "Deallocating the memory..."
    
    IF (ALLOCATED(node))  DEALLOCATE(node)
    IF (ALLOCATED(front)) DEALLOCATE(front)
    IF (ALLOCATED(elem))  DEALLOCATE(elem)
    IF (ALLOCATED(nbr_nodes_per_elem))  DEALLOCATE(nbr_nodes_per_elem)
    IF (ALLOCATED(U0))    DEALLOCATE(U0)
    IF (ALLOCATED(BoundCond)) DEALLOCATE(BoundCond)
    
    
      ! if (allocated(U0))   deallocate(U0)
      ! if (allocated(node)) deallocate(node)
      ! if (allocated(elem)) deallocate(elem)
      ! if (allocated(stencil)) deallocate(stencil)
      ! if (allocated(rhs)) deallocate(rhs)
      ! if (allocated(ia)) deallocate(ia)
      ! if (allocated(ja)) deallocate(ja)
      ! if (allocated(stiffness)) deallocate(stiffness)
      ! if (allocated(mat)) deallocate(mat)
      ! if (allocated(nodCL)) deallocate(nodCL)
      ! if (allocated(nodeglob)) deallocate(nodeglob)
      ! if (allocated(frontglob)) deallocate(frontglob)
      ! if (allocated(elemglob)) deallocate(elemglob)
      ! if (allocated(node)) deallocate(node)
      ! if (allocated(front)) deallocate(front)
      ! if (allocated(elem)) deallocate(elem)
      ! if (allocated(loc2glob)) deallocate(loc2glob)
    
END SUBROUTINE mem_deallocate