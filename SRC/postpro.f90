!-----------------------------------------------------------------------
      subroutine postpro
!-----------------------------------------------------------------------
      use module_icp
      implicit none
      
    ! Local parameters
      INTEGER(ki) :: i,j,ielm
      REAL(kr)    :: tmp1,n1r,n2r,n3r,n1z,n2z,n3z,det,r123,surf,sumsurf(nbrNodes)
      REAL(kr)    :: Er1,Er2,Er3,Ei1,Ei2,Ei3,Er(nbrNodes),Ei(nbrNodes)
      REAL(kr)    :: dissPower,coilPower
      real(kr)    :: na(3),nb(3),nc(3),nd(3)
      real(kr)    :: a,b,c,d,r1,r2,r3
      
      ALLOCATE(Brr(1:nbrNodes)); Brr = zero
      ALLOCATE(Bri(1:nbrNodes)); Bri = zero
      ALLOCATE(Bzr(1:nbrNodes)); Bzr = zero
      ALLOCATE(Bzi(1:nbrNodes)); Bzi = zero
      ALLOCATE(Frad(1:nbrNodes)); Frad = zero
      ALLOCATE(Fax(1:nbrNodes)); Fax = zero
      ALLOCATE(dissJoule(1:nbrNodes)); dissJoule = zero
      
      ! Shape functions at the quadrature points
      na(1:3) = (/1.0d00 , 1.0d00, 1.0d00/) / 3.0d00
      nb(1:3) = (/0.2d00 , 0.2d00, 0.6d00/) 
      nc(1:3) = (/0.2d00 , 0.6d00, 0.2d00/) 
      nd(1:3) = (/0.6d00 , 0.2d00, 0.2d00/) 
      
      sumsurf = zero
        
      ! Loop on all nodes to compute magnetic field
      DO i=1,nbrNodes
        ! Loop on every elements surrounding each node
        DO j=1,stencilElem(i,1)
           ielm = stencilElem(i,j+1)
           r123 = ( node(elem(ielm,1),2) + node(elem(ielm,2),2) + node(elem(ielm,3),2) ) / 3.0d00
           n1z =  node(elem(ielm,3),2) - node(elem(ielm,2),2) ! -r2+r3
           n2z =  node(elem(ielm,1),2) - node(elem(ielm,3),2) ! -r3+r1
           n3z =  node(elem(ielm,2),2) - node(elem(ielm,1),2) ! -r1+r2
           n1r = -node(elem(ielm,3),1) + node(elem(ielm,2),1) ! -z3+z2
           n2r = -node(elem(ielm,1),1) + node(elem(ielm,3),1) ! -z1+z3
           n3r = -node(elem(ielm,2),1) + node(elem(ielm,1),1) ! -z2+z1 
           det = (n3r*n2z)-(n2r*n3z)
           surf = abs(det)*0.5d00
           sumsurf(i) = sumsurf(i) + 1.d00/surf
           
           Er1 = U0(2*elem(ielm,1)-1)
           Er2 = U0(2*elem(ielm,2)-1)
           Er3 = U0(2*elem(ielm,3)-1)
           Ei1 = U0(2*elem(ielm,1))+Ecoils(elem(ielm,1))
           Ei2 = U0(2*elem(ielm,2))+Ecoils(elem(ielm,2))
           Ei3 = U0(2*elem(ielm,3))+Ecoils(elem(ielm,3))
           
        ! dEi/dz
           Brr(i) = Brr(i) + ( n1z*Ei1 + n2z*Ei2 + n3z*Ei3 ) / (det*omega*surf)
     
        ! dEr/dz
           Bri(i) = Bri(i) - ( n1z*Er1 + n2z*Er2 + n3z*Er3 ) / (det*omega*surf)
     
        ! dEi/dr
           Bzr(i) = Bzr(i) - ( n1r*Ei1 + n2r*Ei2 + n3r*Ei3 ) / (det*omega*surf)
           
        ! dEr/dr
           Bzi(i) = Bzi(i) + ( n1r*Er1 + n2r*Er2 + n3r*Er3 ) / (det*omega*surf)
        ENDDO
        
        Er(i) = U0(2*i-1)
        Ei(i) = U0(2*i)
    ENDDO 
    
    ! Average the gradient over the elements
    Brr = Brr / sumsurf
    Bri = Bri / sumsurf
    Bzr = Bzr / sumsurf
    Bzi = Bzi / sumsurf
    ! Add the nodal contributions to axial magnetic field
    where (node(:,2) .GT. zero) Bzr = Bzr - Er / (node(:,2)*omega)
    where (node(:,2) .GT. zero) Bzi = Bzi + Er / (node(:,2)*omega)
    
    ! Local Lorentz  force
    Frad =   0.5d00 * sigma_in * ( Er*Bzr + Ei*Bzi )
    Fax  = - 0.5d00 * sigma_in * ( Er*Brr + Ei*Bri )
    
    ! Local Joule dissipation
    dissJoule = 0.5d00 * sigma_in * (Er*Er + Ei*Ei)
    
    ! dissPower = zero
    ! coilPower = zero
    ! DO ielm=1,nbrElem
       ! If outside the torch, then skip the element
       ! if (elem(ielm,5).ne.2) cycle 
       
       ! r1   = node(elem(ielm,1),2)
       ! r2   = node(elem(ielm,2),2)
       ! r3   = node(elem(ielm,3),2)
       ! r123 = ( r1 + r2 + r3 ) / 3.0d00
      
       ! n2z =  node(elem(ielm,1),2) - node(elem(ielm,3),2) ! -r3+r1
       ! n3z =  node(elem(ielm,2),2) - node(elem(ielm,1),2) ! -r1+r2
       ! n2r = -node(elem(ielm,1),1) + node(elem(ielm,3),1) ! -z1+z3
       ! n3r = -node(elem(ielm,2),1) + node(elem(ielm,1),1) ! -z2+z1 
       ! det = (n3r*n2z)-(n2r*n3z)
       
      ! Weights for the numerical quadrature
       ! a = -27.0d00 * r123 * abs(det) / 96.0d00
       ! b =  25.0d00 * (0.2d00*r1 + 0.2d00*r2 + 0.6d00*r3) * abs(det) / 96.0d00
       ! c =  25.0d00 * (0.2d00*r1 + 0.6d00*r2 + 0.2d00*r3) * abs(det) / 96.0d00
       ! d =  25.0d00 * (0.6d00*r1 + 0.2d00*r2 + 0.0d00*r3) * abs(det) / 96.0d00
      
       ! DO j=1,3
          ! dissPower = dissPower + 2.0d00 * pi * (a*na(j)+b*nb(j)+c*nc(j)+d*nd(j)) * dissJoule(elem(ielm,j))
          ! coilPower = coilPower + omega * pi *  (a*na(j)+b*nb(j)+c*nc(j)+d*nd(j)) * (Brr(elem(ielm,j))**2 + Bzr(elem(ielm,j))**2) / ( 2.d00 * mu0 )
       ! ENDDO
    ! ENDDO
    ! write(*,*) "---------------------------------------------"
    ! write(*,'(a,f9.3)') "  Input coil power              [kW]:",coilPower/1000.0d00
    ! write(*,'(a,f9.3)') "  Dissipated power in the torch [kW]:",dissPower/1000.0d00
    ! write(*,'(a,f9.3)') "  Conversion efficiency             :",100.d00*dissPower/coilPower
    ! write(*,*) "---------------------------------------------"
!-----------------------------------------------------------------------
      end subroutine postpro
!-----------------------------------------------------------------------
