!-----------------------------------------------------------------------
      MODULE module_icp
!-----------------------------------------------------------------------
      implicit none
      save

      integer, parameter :: kr = 8
      integer, parameter :: ki = 4
      real(kr), parameter:: pi = 3.1415926535897932384626434d00
      real(kr), parameter:: zero = 0.0d00
      real(kr), allocatable :: node(:,:)
      real(kr), allocatable :: nodeglob(:,:)
      real(kr), allocatable :: U0(:),rhs(:),mat(:)
      integer(ki), allocatable :: frontglob(:,:),elemglob(:,:)
      integer(ki), allocatable :: front(:,:),elem(:,:),stencil(:,:),stiffness(:,:),nodCL(:,:)
      integer(ki), allocatable :: ia(:),ja(:),ii(:),jj(:),loc2glob(:)
      integer(ki) :: irank,nproc
      character(len=200) :: maillage, file_restart, file_gmsh, file_dat
      integer(ki) :: nbrtotnodes,nbrtotelem,nbrtotfront
      integer(ki) :: nbrNodes, nbrElem, nbrFront,nummaxnodes,nummaxelm
      integer(ki) :: nbvar, nstencil,maxstencil=10, nodelm=3
      integer(ki) :: time1,time2
      integer(ki) :: time_begin(8), time_end(8), time_diff(8), seconds, milliseconds,numnz
      real(kr)    :: eps

!-----------------------------------------------------------------------
      end module module_icp
!-----------------------------------------------------------------------
