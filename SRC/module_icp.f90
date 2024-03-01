!-----------------------------------------------------------------------
      MODULE module_icp
!-----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE

      INTEGER, PARAMETER :: kr = 8
      INTEGER, PARAMETER :: ki = 4
      
      INTEGER(ki), PARAMETER:: nbvar = 2
      REAL(kr), PARAMETER:: pi = 4.0d00*datan(1.0d00)
      REAL(kr), PARAMETER:: zero = 0.0d00
      
      REAL(kr), ALLOCATABLE :: node(:,:)
      REAL(kr), ALLOCATABLE :: nodeglob(:,:)
      REAL(kr), ALLOCATABLE :: U0(:),rhs(:),mat(:),Ecoils(:),sigma_in(:)
      ! REAL(kr), ALLOCATABLE :: BoundCond(:,:)   ! Imposed boundary conditions
      INTEGER(ki), ALLOCATABLE :: nbr_nodes_per_elem(:) ! see read_gmsh
      INTEGER(ki), ALLOCATABLE :: frontglob(:,:),elemglob(:,:)
      INTEGER(ki), ALLOCATABLE :: front(:,:),elem(:,:),stencil(:,:),stiffness(:,:),nodCL(:,:)
      INTEGER(ki), ALLOCATABLE :: ia(:),ja(:)
      ! INTEGER(ki), ALLOCATABLE :: ii(:),jj(:),loc2glob(:)
      INTEGER(ki) :: irank,nproc
      
      REAL(kr)    :: mu0,omega,Icoil,frequency
      REAL(kr), ALLOCATABLE    :: coils(:,:) ! coils(3,2)
      
      INTEGER(ki), PARAMETER :: length_names = 40
      CHARACTER(LEN=length_names) :: mesh_file, file_gmsh
      
      INTEGER(ki) :: nbrNodes, nbrElem, nbrTris, nbrQuads, nbrFront
      INTEGER(ki) :: nbrInt, nbrBC,nbrCoils
    
      INTEGER(ki), DIMENSION(2,5) :: CLTable
      INTEGER(ki) :: nbrtotnodes,nbrtotelem,nbrtotfront
      INTEGER(ki) :: nummaxnodes,nummaxelm
      INTEGER(ki) :: nstencil,maxstencil=12, nodelm=3
      INTEGER(ki) :: time1,time2
      INTEGER(ki) :: time_begin(8), time_end(8), time_diff(8), seconds, milliseconds,numnz
      REAL(kr)    :: eps

!-----------------------------------------------------------------------
      END MODULE module_icp
!-----------------------------------------------------------------------