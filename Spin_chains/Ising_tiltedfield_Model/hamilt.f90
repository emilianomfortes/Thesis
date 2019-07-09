MODULE hamilt

IMPLICIT NONE

CONTAINS

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE ising_chain(GG,alpha_x,alpha_y,alpha_z,HH,sites,BC,neig)

    ! alpha's: are couplings
    ! GG: onsite energy
    ! HH: (output) is the Ising hamiltonian
    ! sites: #of sites
    ! BC: boundary conditions 'periodic' or 'open'
    ! neig: 1 neighbour interaction, 2 next-neighbour interaction, etc.

    IMPLICIT NONE

    INTEGER :: i,n,j,s,info,ldwork,sites,dim,i1,neig
    REAL(8) :: GG,alpha_z,alpha_y,alpha_x,m1,HH(:,:)
    CHARACTER(10) :: BC                             !boundary condition
    COMPLEX(8) :: ii

    REAL(8), ALLOCATABLE, SAVE :: Hint_Z(:,:),Hint_X(:,:),Hint_Y(:,:),He(:,:),Hamil(:,:)
    REAL(8), ALLOCATABLE, SAVE :: coupl_x(:),coupl_y(:),ener_sites(:),coupl_z(:)

!!$   COMPLEX(8), ALLOCATABLE, SAVE :: work(:)
!!$   REAL(8), ALLOCATABLE, SAVE :: rwork(:),eigen1(:)
!!$   REAL(8), ALLOCATABLE, SAVE :: Had(:,:)
!!$   COMPLEX(8), ALLOCATABLE, SAVE :: Had_Y(:,:)

    ii=(0.d0,1.d0)

    dim = 2**sites

    ALLOCATE(He(dim,dim),Hint_Z(dim,dim),Hint_X(dim,dim),Hint_Y(dim,dim),Hamil(dim,dim))
    ALLOCATE(coupl_x(sites),coupl_y(sites),coupl_z(sites),ener_sites(sites))
!!$   ALLOCATE(Had_Y(dim,dim),Had(dim,dim))

    He = 0.d0   !ENERGY OF EACH SITE
    Hint_Z = 0.d0
    Hint_X = 0.d0
    Hint_Y = 0.d0

    !   print*,'START HAMILT CONSTRUCTION'

    ! COUPLINT ENERGIES -------------------
    coupl_x = 1.d0/2.d0
    coupl_y = 1.d0/2.d0
    coupl_z = 1.d0/2.d0
    ener_sites = 1.d0/2.d0

!!$   CALL RANDOM(coupl_x,sites)
!!$   coupl_x = coupl_x*3.d0
!!$   CALL RANDOM(coupl_y,sites)
!!$   CALL RANDOM(coupl_z,sites)
!!$   CALL RANDOM(ener_sites,sites)
!!$   ener_sites = ener_sites*4.d0

    ! -------------------------------------

    DO i=0,2**sites-1
       DO n=0,sites-1
          He(i+1,i+1) =  He(i+1,i+1) + ener_sites(n+1)*(-1.d0)**IBITS(i,n,1) !on-site energy
          ! interaction
          IF ((n <= sites-1-neig).OR.(BC.EQ.'periodic')) then
             i1 = IEOR(i,IBSET(0,n)) !this is a NOT in the site n of i (x and y aplication)
             i1 = IEOR(i1,IBSET(0,MOD(n+neig,sites)))  !this is a NOT in the site n+1 of i
             Hint_Y(i+1,i1+1) =   Hint_Y(i+1,i1+1) +   &
                  coupl_y(n+1)*(-1.d0)*(-1.d0)**(IBITS(i,n,1)+IBITS(i,MOD(n+neig,sites),1)) ! i*(-1)^spin_i i*(-1)^spin_(i+1)
             Hint_X(i+1,i1+1) =  Hint_X(i+1,i1+1) + coupl_x(n+1)         
             Hint_Z(i+1,i+1) =  Hint_Z(i+1,i+1) + coupl_z(n+1)*(-1.d0)**(IBITS(i,n,1)+IBITS(i,MOD(n+neig,sites),1)) 
          END IF
!!$!---------------------------------------
!!$      do j = i,2**sites-1 
!!$         m = 0.d0  
!!$         do n=0,sites-1
!!$            IF (btest(iand(i,j),n)) THEN !construction of the Haddamard matrix
!!$               m1 = m1 + 1.d0
!!$            END IF            
!!$         end do
!!$         Had(i+1,j+1) = (-1.d0)**m1
!!$         Had_Y(i+1,j+1) = (ii)**m1
!!$      end do
!!$      Had(i+1,i+1) = Had(i+1,i+1)/2.d0  !this is because at the end I add the transpose
!!$      Had_Y(i+1,i+1) = Had_Y(i+1,i+1)/2.d0
!!$!--------------------------------------
       END DO
    END DO

!!$   Had = (Had + TRANSPOSE(Had))/sqrt(2.d0**sites)
!!$   Had_Y = (Had_Y + TRANSPOSE(CONJG(Had_Y)))/sqrt(2.d0**sites)   
!!$   Hamil = GG*He + alpha_z*Hint_Z + alpha_x*MATMUL(Had,MATMUL(Hint_Z,Had)) +  &
!!$       alpha_y*MATMUL(Had_Y,MATMUL(Hint_Z,TRANSPOSE(CONJG(Had_Y))))

    Hamil = GG*He + alpha_z*Hint_Z + alpha_x*Hint_X + alpha_y*Hint_Y

    HH = Hamil

    DEALLOCATE(He,Hint_Z,Hint_X,Hint_Y,Hamil,coupl_x,coupl_y,coupl_z,ener_sites)

  end subroutine ising_chain


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE ising_chain_FN(GG,alpha_x,alpha_y,alpha_z,sites,n_partic,BC,neig,HH,state,flag)
    ! This routine generates the Hamiltonian: HH of a s=1/2 chain with fixed number of particles.
    ! OUTPUTS: state: gives the states with the binary notation given by the init_state routine
    ! flag: give the flag for this states flag -> state
    ! INPUTS: BC boundary conditions 'OPEN' 'PERIODIC'
    ! neig: neigbour interaction +1,+2, etc.
    ! sites: # sites, n_part: # particles. The other variables are the couplings.

    IMPLICIT NONE
    INTEGER :: sites,n_partic,j,label_state,k,i,dim
    INTEGER :: n,s,info,i1,neig
    INTEGER :: state(:),flag(:)
    REAL(8) :: GG,alpha_z,alpha_y,alpha_x,m1,HH(:,:)
    CHARACTER(10) :: BC                             !boundary condition

    REAL(8), ALLOCATABLE, SAVE :: Hint_Z(:,:),Hint_X(:,:),Hint_Y(:,:),He(:,:),Hamil(:,:)
    REAL(8), ALLOCATABLE, SAVE :: coupl_x(:),coupl_y(:),ener_sites(:),coupl_z(:)

    !dim = 2**sites
    dim = int(exp(logfac(sites)-logfac(n_partic)-logfac(sites-n_partic)))

    ALLOCATE(He(dim,dim),Hint_Z(dim,dim),Hint_X(dim,dim),Hint_Y(dim,dim),Hamil(dim,dim))
    ALLOCATE(coupl_x(sites),coupl_y(sites),ener_sites(sites),coupl_z(sites))

    CALL init_states(sites,n_partic,state,flag)

    PRINT*,'INISTATE'

    He = 0.d0   !ENERGY OF EACH SITE
    Hint_Z = 0.d0
    Hint_X = 0.d0
    Hint_Y = 0.d0

    ! COUPLINT ENERGIES -------------------
    coupl_x = 1.d0/1.d0
    coupl_y = 1.d0/1.d0
    coupl_z = 1.d0/1.d0
    ener_sites = 1.d0/4.d0

!!$   CALL RANDOM(coupl_x,sites)
!!$   coupl_x = coupl_x*3.d0
!!$   CALL RANDOM(coupl_y,sites)
!!$   CALL RANDOM(coupl_z,sites)
!!$   CALL RANDOM(ener_sites,sites)
!!$   ener_sites = ener_sites*4.d0
    ! -------------------------------------

    DO i=1,dim
       DO n=0,sites-1
          He(i,i) =  He(i,i) + ener_sites(n+1)*(-1.d0)**IBITS(state(i),n,1) !on-site energy
          ! interaction
          IF ((n <= sites-1-neig).OR.(BC.EQ.'periodic')) THEN
             IF ((IBITS(state(i),n,1)+IBITS(state(i),MOD(n+neig,sites),1)).EQ.1) THEN
                i1 = IEOR(state(i),IBSET(0,n)) !this is a NOT in the site n of i (x and y aplication), this is state is not in the subspace!
                i1 = flag(IEOR(i1,IBSET(0,MOD(n+neig,sites))))  !this is a NOT in the site n+neig of i
                Hint_Y(i,i1) = Hint_Y(i,i1) +  &
                     coupl_y(n+1)*(-1.d0)*(-1.d0)**(IBITS(state(i),n,1)+IBITS(state(i),MOD(n+neig,sites),1)) ! i*(-1)^spin_i i*(-1)^spin_(i+1)
                Hint_X(i,i1) =  Hint_X(i,i1) + coupl_x(n+1)        
                Hint_Z(i,i) = Hint_Z(i,i) + coupl_z(n+1)*(-1.d0)**(IBITS(state(i),n,1)+IBITS(state(i),MOD(n+neig,sites),1)) 
             END IF
          END IF
       END DO
    END DO

    Hamil = GG*He + alpha_z*Hint_Z + alpha_x*Hint_X + alpha_y*Hint_Y

    HH = Hamil

    DEALLOCATE(He,Hint_Z,Hint_X,Hint_Y,Hamil,coupl_x,coupl_y,coupl_z,ener_sites)


  END SUBROUTINE ising_chain_FN

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  SUBROUTINE init_states(sites,n_partic,state,flag)

    ! Creates a labeling for the states
    ! states: a vector where the elements are the states in binary notation with '# sites' bits
    ! flags: creates flags for the states, so if we have the binary decomposition it gives the location in states.
    ! The # of states is : sites!/(n_part! (sites-n_part)!)
    ! dim = int(exp(logfac(sites)-logfac(n_partic)-logfac(sites-n_partic)))

    IMPLICIT NONE

    INTEGER :: sites,n_partic,j,label_state,k,i,dim
    INTEGER :: state(:),flag(:)
    !INTEGER, ALLOCATABLE, SAVE :: state(:),flag(:)

    dim = int(exp(logfac(sites)-logfac(n_partic)-logfac(sites-n_partic)))
    !dim = 2**sites

    !  ALLOCATE(state(dim),flag(2**sites))

    flag = -1
    label_state = 0
    DO i=0,2**sites-1
       k = 0
       j = 0
       DO WHILE (k<=n_partic.AND.j<=sites-1)
          k = k + IBITS(i,j,1)
          j = j + 1
       END DO
       IF (k.eq.n_partic) THEN
          label_state = label_state + 1
          state(label_state) = i
          flag(i) = label_state  !i=0 never appears, we do not use cero particles
       END IF

    END DO

  END SUBROUTINE init_states

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE stat_level_spacing(bins,n_cut,energy,width,histogram,dim)

    ! IT GIVES THE HISTOGRAM OF THE LEVEL SPACING, OUTPUTS: width, histogram

    IMPLICIT NONE

    INTEGER :: i,n,j,s,m,info,ldwork,dim,i1,bins,k,n_cut
    REAL(8) :: GG,alpha_z,alpha_y,alpha_x,m1,energy(:),delta,width(:),histogram(:)
    CHARACTER(10) :: BC                             !boundary condition
    COMPLEX(8) :: ii

    REAL(8), ALLOCATABLE, SAVE :: space(:)
    INTEGER, ALLOCATABLE, SAVE :: flag(:),flag1(:)

    ALLOCATE(space(dim),flag(dim),flag1(dim))

    space = CSHIFT(energy,1) - energy
    space = CSHIFT(space,n_cut/2)        ! we put the low and high energies near, and cut n_cut of them

!!$   do i=1,dim
!!$     print*, space(i)
!!$   enddo
!!$   stop

    space(dim-n_cut:dim) = 0.d0
    space = space/(sum(space(1:dim-n_cut))/(dim-real(n_cut)))   !normalized level spacing


    delta = maxval(space)/bins

    width = (/((i-1)*delta,i=1,bins)/)
    flag = (/(i,i=1,dim)/)
    histogram = 0.d0


    print*,delta*bins,minval(space)
    j = 1
    s = dim - n_cut
    DO WHILE (j<=bins.AND.s>0)
       k = 0
       flag1 = 0
       DO i = 1,s          
          IF((delta*(j-1) <= space(flag(i))).AND.(space(flag(i)) < delta*j)) THEN
             histogram(j) = histogram(j) + 1.D0
          ELSE
             k = k + 1
             flag1(k) = flag(i)
          END IF
       END DO
       flag = flag1
       s = k
       j = j+1
    END DO

    histogram = histogram/(sum(histogram)*delta) !normalized histogram

    DEALLOCATE(space,flag,flag1)



  END SUBROUTINE stat_level_spacing

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE density_states(bins,energy,width,histogram,dim)

    ! it gives the histogram of the density of states width is the energy

    IMPLICIT NONE

    INTEGER :: i,n,j,s,m,dim,bins,k
    REAL(8) :: energy(:),delta,width(:),histogram(:)

    delta = (energy(dim)-energy(1))/bins
    print*,delta

    width = (/(energy(1)+(i-1)*delta,i=1,bins)/)

    histogram = 0.d0   
    j = 1
    DO k = 1,bins
       DO WHILE (energy(j) < width(k).AND.j < dim)
          histogram(k) = histogram(k) + 1.D0
          j=j+1
       END DO
    END DO
    !   print*,sum(histogram),dim,delta

    histogram = histogram /delta  !states over delta energy


  END SUBROUTINE density_states

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE strength_function(bins,n_state,n_average,HH,HH0,energy0,width,histogram,dim)

    ! it gives the histogram of the density of states width is the energy

    IMPLICIT NONE

    INTEGER :: i,n,j,s,m,dim,i1,bins,k,n_average,n_state,mean_value
    REAL(8) :: energy0(:),histogram(:),delta,width(:)
    REAL(8) :: HH(:,:),HH0(:,:)

    REAL(8), ALLOCATABLE, SAVE :: sf_function(:,:),sf_bins(:,:),density(:)

    ALLOCATE(sf_function(dim,n_average+1),sf_bins(bins,n_average+1),density(bins))


    ! histogram is the density of states
    CALL density_states(bins,energy0,width,density,dim)

    delta = (energy0(dim)-energy0(1))/bins

    sf_function = 0.d0
    DO i=1,dim    !overlaps with all states
       DO j=1,n_average
          sf_function(i,j) = ABS(DOT_PRODUCT(HH(:,i),HH0(:,n_state+j)))**2!*histogram(i)
       END DO
    END DO

    DO j=1,n_average
       mean_value = INT(DOT_PRODUCT(sf_function(:,j),(/(i,i=1,dim)/))/SUM(sf_function(:,j)))
       print*,mean_value
       sf_function(:,j) = CSHIFT(sf_function(:,j),mean_value)
       sf_function(:,n_average+1) = sf_function(:,n_average+1) + sf_function(:,j)/REAL(n_average)
    ENDDO
    mean_value = INT(DOT_PRODUCT(sf_function(:,n_average+1),(/(i,i=1,dim)/))/(SUM(sf_function(:,n_average+1))))
    sf_function(:,n_average+1) = CSHIFT(sf_function(:,n_average+1),-mean_value)


    width = (/(energy0(1)+(i-1)*delta,i=1,bins)/)

    histogram = 0.d0   
    j=1
    DO k = 1,bins
       DO WHILE (energy0(j) < width(k))
          histogram(k) = histogram(k) + sf_function(j,n_average+1)
          j=j+1
       END DO
    END DO


    do i=1,bins
       histogram(i) = histogram(i)*density(i)
    end do

    histogram = histogram/delta/sum(histogram)


    DEALLOCATE(sf_function,sf_bins,density)

  END SUBROUTINE strength_function

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  SUBROUTINE Sdiag_IPR(HH,HH0,Sdiag,IPR,dim)

    ! it gives the Sdiag and the IPR of HH as a function of HH0

    IMPLICIT NONE

    INTEGER :: i,n,j,s,m,dim,i1,bins,k
    REAL(8) :: Sdiag(:),IPR(:)
    REAL(8) :: HH(:,:),HH0(:,:)

    REAL(8), ALLOCATABLE, SAVE :: overlaps(:,:)

    ALLOCATE(overlaps(dim,dim))
 
    overlaps = ABS(MATMUL(TRANSPOSE(HH),HH0))**2

    DO i = 1,dim
       Sdiag(i) = -DOT_PRODUCT(overlaps(i,:),log(overlaps(i,:)))
       IPR(i) =1.d0/DOT_PRODUCT(overlaps(i,:),overlaps(i,:))
    END DO

    DEALLOCATE(overlaps)

  END SUBROUTINE Sdiag_IPR

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  SUBROUTINE Over_neighb(HH,HH0,over,dim)

    ! it gives the Sdiag and the IPR of HH as a function of HH0

    IMPLICIT NONE

    INTEGER :: i,dim
    REAL(8) :: HH(:,:),HH0(:,:),over(:)

    REAL(8), ALLOCATABLE, SAVE :: overlaps(:,:)

    ALLOCATE(overlaps(dim,dim))
 
    overlaps = ABS(MATMUL(TRANSPOSE(HH),HH0))**2

    over = 0.d0
    DO i = 1,dim-1
       over(i) =DOT_PRODUCT(overlaps(i,:),overlaps(i+1,:))
    END DO

    DEALLOCATE(overlaps)

  END SUBROUTINE Over_neighb

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE RANDOM(string,nn)
    ! Generates random string of nn elements between 0 and 1

    INTEGER, DIMENSION(1) :: OLD ! THIS PROGRAM ASSUMES K = 1
    INTEGER ::  K,x,TIMEARRAY(3),i,a,b,fr,nn,J,n, clock
    REAL(8) :: HARVEST,string(:)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    n=3
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DO I=1,NN
       CALL RANDOM_NUMBER(HARVEST)
       string(i) = HARVEST
    ENDDO

  END SUBROUTINE RANDOM
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  RECURSIVE FUNCTION factorial( n ) RESULT(res)
    INTEGER, INTENT(IN) :: n
    INTEGER :: res
    IF( n==0 ) THEN
       res = 1
    ELSE
       res = n*factorial( n-1 )
    END IF
  END FUNCTION factorial
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  FUNCTION logfac(n1)

!!$   Calculates the log of a factorial

    INTEGER :: n1,i1
    REAL(8) :: logfac
    REAL(8), ALLOCATABLE, SAVE :: vec(:)
    ALLOCATE(vec(n1))

    vec = (/ (i1, i1 = 1, n1) /)

    logfac = SUM(DLOG(vec))

    IF (n1.EQ.0) logfac = 0.

    DEALLOCATE(vec)

  END FUNCTION logfac

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module hamilt
