!-----------------------------------------------------------------------
      SUBROUTINE ecriture_solution(ok)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      integer(ki), INTENT(out) :: ok
      integer(ki) :: i, ierr

      ok = 1
      open(unit=20,file=file_dat,status="replace",iostat=ierr,form='formatted')
      if (ierr /= 0) then
         write(*,*) "Error writing the solution in .dat format"
         ok = 0
      end if

      do i=1,nbrElem*nbvar
         write(20,'(ES15.8)') U0(i)
      end do

      close(unit=20)

!-----------------------------------------------------------------------
      end subroutine ecriture_solution
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE ecriture_gmsh
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      
      integer(ki) :: ierr, i,j
      write(*,*) "Writing the solution in gmsh format"
      open(unit=10,file=file_gmsh,status="replace",iostat=ierr,form='formatted')
      write(10,'(T1,A11)') "$MeshFormat"
      write(10,'(T1,A7)') "2.2 0 8"
      write(10,'(T1,A14)') "$EndMeshFormat"
      write(10,'(T1,A6)') "$Nodes"
      write(10,'(T1,I8)') nbrtotNodes
      do i=1,nbrtotNodes
         write(10,'(T1,I8,2F15.8,1X,A2)') i, nodeglob(1:2,i), "0."
      end do
      write(10,'(T1,A9)') "$EndNodes"
      write(10,'(T1,A9)') "$Elements"
      write(10,'(T1,I8)') nbrtotElem
      do i=1,nbrtotElem
         write(10,'(T1,I8,1X,4I2,1X,I8,1X,I8,1X,I8)') i,2,2,7,6,(elemglob(j,i),j=1,3)!elemglob(2,i),elemglob(3,i)
      end do
      write(10,'(T1,A12)') "$EndElements"

      write(10,'(T1,A9)') "$NodeData"
      write(10,'(T1,A1)') "1"
      write(10,'(T1,A13)') '"Temperature"'
      write(10,'(T1,A1)') "1"
      write(10,'(T1,A1)') "0"
      write(10,'(T1,A1)') "3"
      write(10,'(T1,A1)') "0"
      write(10,'(T1,A1)') "1"
      write(10,'(T1,I8)') nbrtotNodes
      do i=1,nbrtotNodes
         if (U0(i)<eps) U0(i)=eps
         write(10,'(T1,I8,ES15.8)') i, U0(i)
      end do
      write(10,'(T1,A12)') "$EndNodeData"
      close(unit=10)

!-----------------------------------------------------------------------
      end subroutine ecriture_gmsh
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE lecture_gmsh(ok)
!-----------------------------------------------------------------------
      use module_icp
      implicit none

      ! Parameters of the subroutine
      integer(ki), INTENT(out) :: ok

      ! Local parameters
      character(len=256)    :: line
      integer(ki) :: i, ierr, a, b, c, nbrElemTot=0,nbrnoeudfront=0
      real(kr) :: header
      logical :: type2 = .true.
      ! external :: find_elem_front

      call sampletime(time1)
      ok = 1
      write(*,*) "-------------------------------------------------------"
      write(*,*) "Name of the mesh data : ", trim(adjustl(maillage))

      open(unit=10,file=maillage,status="old",iostat=ierr,form='formatted')
      if (ierr /= 0) then
         write(*,*) "Error occuring : file '",trim(maillage),"' not found "
         ok = 0
      end if

      read(10,*) line, header, b, c

      if (header < 2.) then
         write(*,*) "Error : wrong format for gmsh data"
         ok = 0
      end if

      if (header >2.1) type2 = .false.

      do while (line /="$Nodes")
         read(10,*) line
      end do

      read(10,*) nbrtotNodes
      write(*,101) nbrtotNodes

      allocate(nodeglob(1:2,1:nbrtotNodes))

      do i=1,nbrtotNodes ! Lecture des coordonnées des noeuds
         read(10,*) a, nodeglob(1,i), nodeglob(2,i), b
         if (a/=i) write(*,*) "merde dans le maillage"
      end do

      do while (line /="$Elements")
         read(10,*) line
      end do

      read(10,*) nbrElemTot
      read(10,*) a,i
      do while (i/=2) ! parcourir une 1ere fois les elements de frontiere pour en determiner leur nombre
         if (i==1) then
            nbrtotFront = nbrtotFront + 1
         else if (i==15) then
            nbrnoeudfront = nbrnoeudfront + 1
         end if
         read(10,*) a,i
      end do
      nbrtotelem = nbrElemTot - nbrtotFront - nbrnoeudfront

      write(*,102) nbrtotFront
      write(*,103) nbrtotelem
      allocate(frontglob(1:3,1:nbrtotFront))
      allocate(elemglob(1:4,1:nbrtotElem))

      do i=1,nbrtotFront+1
         backspace(10)
      end do

      if (type2) then
         do i=1,nbrtotFront ! Lecture des elements de frontiere
            read(10,*) a,a,a,frontglob(3,i),a,a,frontglob(1,i),frontglob(2,i)
         end do
      else
         do i=1,nbrtotFront ! Lecture des elements de frontiere
            read(10,*) a,a,a,frontglob(3,i),a,frontglob(1,i),frontglob(2,i)
         end do
      end if

      if (type2) then
         do i=1,nbrtotElem ! Lecture des elements de surface
            read(10,*) a,a,a,elemglob(4,i),a,a,elemglob(1,i),elemglob(2,i),elemglob(3,i)
         end do
      else
         do i=1,nbrtotElem ! Lecture des elements de surface
            read(10,*) a,a,a,a,elemglob(4,i),elemglob(1,i),elemglob(2,i),elemglob(3,i)
         end do
      end if

      call setCLlink

      call sampletime(time2)
      call time_display
      write(*,*) "-------------------------------------------------------"
101   format(" Number of nodes : ",T35,I7)
102   format(" Number of border elements :", T35, I7)
103   format(" Number of surface elements :",T35,I7)
      close(unit=10)
!-----------------------------------------------------------------------
      END SUBROUTINE lecture_gmsh
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine setCLlink
!-----------------------------------------------------------------------
! Créer un lien entre noeud de frontière et tag de condition au limite
! Dans le maillage, la frontière est créée dans le sens anti-horlogique
      use module_icp
      implicit none

      integer(ki) ::i

      allocate(nodCL(1:nbrtotFront,1:2))

      do i=1,nbrtotFront
         nodCL(i,1) = frontglob(1,i)
         nodCL(i,2) = frontglob(3,i)
      end do
!-----------------------------------------------------------------------
      end subroutine setCLlink
!-----------------------------------------------------------------------
