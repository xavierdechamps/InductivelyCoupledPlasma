!-----------------------------------------------------------------------
      MODULE module_icp
!-----------------------------------------------------------------------
      implicit none
      save

      integer, parameter :: kr = 8
      integer, parameter :: ki = 4
      
      INTEGER(ki), PARAMETER:: nbvar = 2
      real(kr), parameter:: pi = 3.1415926535897932384626434d00
      real(kr), parameter:: zero = 0.0d00
      
      real(kr), allocatable :: node(:,:)
      real(kr), allocatable :: nodeglob(:,:)
      real(kr), allocatable :: U0(:),rhs(:),mat(:)
      REAL(kr), ALLOCATABLE :: BoundCond(:,:)   ! Imposed boundary conditions
      INTEGER(ki), ALLOCATABLE :: nbr_nodes_per_elem(:) ! see read_gmsh
      integer(ki), allocatable :: frontglob(:,:),elemglob(:,:)
      integer(ki), allocatable :: front(:,:),elem(:,:),stencil(:,:),stiffness(:,:),nodCL(:,:)
      integer(ki), allocatable :: ia(:),ja(:),ii(:),jj(:),loc2glob(:)
      integer(ki) :: irank,nproc
      
      INTEGER(ki), PARAMETER :: length_names = 50
      CHARACTER(LEN=length_names) :: mesh_file, file_restart, file_gmsh, file_dat
      ! character(len=200) :: maillage, file_restart, file_gmsh, file_dat
      
      INTEGER(ki) :: nbrNodes, nbrElem, nbrTris, nbrQuads, nbrFront, nbrInt, nbrBC
    
      integer(ki) :: nbrtotnodes,nbrtotelem,nbrtotfront
      integer(ki) :: nummaxnodes,nummaxelm
      integer(ki) :: nstencil,maxstencil=10, nodelm=3
      integer(ki) :: time1,time2
      integer(ki) :: time_begin(8), time_end(8), time_diff(8), seconds, milliseconds,numnz
      real(kr)    :: eps

!-----------------------------------------------------------------------
      end module module_icp
!-----------------------------------------------------------------------


MODULE module_mem_allocate
  INTERFACE 
  
! -----------------------
! Subroutine mem_allocate
! -----------------------
  SUBROUTINE mem_allocate(node,front,elem,nbr_nodes_per_elem,U0,BoundCond,&
&                       nbrNodes,nbrElem,nbrFront,only_mesh)

    USE module_icp, only : kr,ki,nbvar
    IMPLICIT NONE

    REAL(kr), ALLOCATABLE    :: U0(:)
    REAL(kr), ALLOCATABLE    :: node(:,:)
    REAL(kr), ALLOCATABLE    :: BoundCond(:,:)
    INTEGER(ki), ALLOCATABLE    :: elem(:,:)
    INTEGER(ki), ALLOCATABLE    :: front(:,:)
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elem(:)
    
    INTEGER(ki), INTENT(IN) :: nbrNodes,nbrElem,nbrFront,only_mesh
        
  END SUBROUTINE mem_allocate
  
! -------------------------
! Subroutine mem_deallocate
! -------------------------
  SUBROUTINE mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,BoundCond)
    USE module_icp, only : kr,ki
    IMPLICIT NONE

    REAL(kr), ALLOCATABLE    :: U0(:)
    REAL(kr), ALLOCATABLE    :: node(:,:)
    REAL(kr), ALLOCATABLE    :: BoundCond(:,:)
    INTEGER(ki), ALLOCATABLE    :: elem(:,:)
    INTEGER(ki), ALLOCATABLE    :: front(:,:)
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elem(:)
        
  END SUBROUTINE mem_deallocate
  
  END INTERFACE 
  
END MODULE module_mem_allocate