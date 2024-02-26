!-----------------------------------------------------------------------
      subroutine sendall(vect1,vect2,vect3)
!-----------------------------------------------------------------------
      use module_icp,only : nproc,nbrelem,nbrnodes,nbrfront,nbrtotnodes,nbrtotelem,nbrtotfront
      implicit none

      include "mpif.h"

      integer :: vect1(1:nproc),vect2(1:nproc),vect3(1:nproc)
      integer :: ierr

      call mpi_scatter(vect1(1:nproc),1,mpi_integer,nbrelem,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_scatter(vect2(1:nproc),1,mpi_integer,nbrnodes,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_scatter(vect3(1:nproc),1,mpi_integer,nbrfront,1,mpi_integer,0,mpi_comm_world,ierr)

      call mpi_bcast(nbrtotnodes,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(nbrtotelem,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(nbrtotfront,1,mpi_integer,0,mpi_comm_world,ierr)

!-----------------------------------------------------------------------
      end subroutine sendall
!-----------------------------------------------------------------------


! ------------------------------------------------------------------
      subroutine sendall2(vect4)
!-----------------------------------------------------------------------
      use module_icp,only : nproc,nbrelem,nbrnodes,nbrfront,nbrtotnodes,nbrtotelem,nbrtotfront,elem,node,front
      implicit none
      
      include "mpif.h"
      
      integer :: vect4(1:nbrtotelem,1:3,0:nproc-1)
      integer i

!!$do i = 0,nproc-1
!!$   if (irank == 0) then
!!$      call mpi_send()
!!$   elseif (irank == i) then
!!$      call mpi_recv()
!!$   endif
!!$enddo

!-----------------------------------------------------------------------
      end subroutine sendall2
!-----------------------------------------------------------------------
