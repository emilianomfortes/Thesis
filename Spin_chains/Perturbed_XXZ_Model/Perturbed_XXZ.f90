MODULE global
  
  COMPLEX(8), SAVE :: ii=(0.d0,1.d0)
  REAL(8), SAVE :: pi=4.d0*atan2(1.d0,1.d0)     
  
END MODULE global

PROGRAM ISING

  USE HAMILT
  USE GLOBAL

  IMPLICIT NONE

  ! PARA LAS RUTINAS LAPACK USAR COMPLEX(8) !!!!!!!!!!!!!!!!!!!

  INTEGER :: i,j,k,l,p,q,sites,dim,info,bins,n_cut,neig,n_partic,dim_parity,dim_parity0,mean_value,dim_2
  INTEGER :: n_state,n_average,T1(3),T2(3),TIME_1,TIME_2,opera
  REAL(8) :: GG,alpha_x,alpha_y,alpha_z,a(4),mu,ele,alpha2_x,alpha2_y,alpha2_z,tau
  CHARACTER(20) ::  boundary_con

  INTEGER, ALLOCATABLE, SAVE :: IWORK(:),state(:),flag(:),statepp(:),statepn(:)
  REAL(8), ALLOCATABLE, SAVE :: width(:),histo(:),ener_parity(:),ener_vec(:),energy(:),sf_function(:,:),op(:,:)
  REAL(8), ALLOCATABLE, SAVE :: HH(:,:),HH0parity(:,:),HH1(:,:),HHparity(:,:),HH0(:,:),HH01(:,:),overlaps(:)
  REAL(8), ALLOCATABLE, SAVE :: HH11(:,:),HH011(:,:),HH12(:,:),Tr(:,:),ipr(:),HHparidad(:,:), shaper(:)
  REAL(8), ALLOCATABLE, SAVE ::energy0(:),work(:), vec1(:),vec(:),ener0_vec(:),ener0_parity(:),Sdiag(:)

  OPEN(8,FILE = 'datos.in',status='unknown') !Reads chain data file
  
  read(8,*)sites,n_partic ! sites = SPIN SITES ; n_partic = NUMBER OF SPINS UP
  read(8,*)boundary_con ! OPEN/PERIODIC
  read(8,*)alpha_x,alpha_y,alpha_z  ! neighbor interaction couplings
  read(8,*)mu ! next-nearest-neighbor interaction perturbation strength
  read(8,*)tau,opera

  dim_2 = 2**sites 
  if (MOD(sites,2) .eq. 0) then 
     dim = int(exp(logfac(sites)-logfac(n_partic)-logfac(sites-n_partic))) 
  else
     dim = int(exp(logfac(sites)-logfac(n_partic)-logfac(sites-n_partic)))
  end if
  
  PRINT*,'dim= ',dim

  OPEN(5,FILE = 'hamilt.out', status='unknown')
  OPEN(7,FILE = 'energias.out', status='unknown')
  OPEN(8,FILE = 'energiasfull.out', status='unknown')
  OPEN(13,FILE = 'states.out', status = 'unknown')
  OPEN(14,file = 'ipr.out' , status = 'unknown')
  OPEN(15,FILE = 'hamilt0full.out', status='unknown')
  open(17,file = 'hamiltfull.out', status = 'unknown')
  OPEN(18,FILE = 'hamiltodd.out',status='unknown')
  OPEN(25,FILE = 'HHnodiag.out',status='unknown')

  ALLOCATE(HH(dim,dim),energy(dim),vec1(dim),width(bins),HH1(dim,dim),vec(dim),histo(bins),HH0(dim,dim),Tr(dim,dim))
  ALLOCATE(work(3*dim**2),IWORK(6*dim),state(dim),flag(2**sites),ener_vec(dim),HHparity(dim,dim),statepp(dim),statepn(dim))
  ALLOCATE(HH0parity(dim,dim),energy0(dim),HH01(dim,dim),ener0_vec(dim))
  ALLOCATE(HH011(dim_2,dim_2),HH11(dim_2,dim_2),HH12(dim_2,dim_2))

  CALL ITIME(T1)
  TIME_1=T1(1)*3600+T1(2)*60+T1(3)
  
  GG = 0.D0 
  PRINT*, mu,alpha_z
  neig = 1
  CALL ising_chain_FN(GG,alpha_x,alpha_y,alpha_z,sites,n_partic,boundary_con,neig,HH0,state,flag)
  neig = 2
  alpha2_x = alpha_x
  alpha2_y = alpha_y
  alpha2_z = alpha_z
  CALL ising_chain_FN(GG,alpha2_x,alpha2_y,alpha2_z,sites,n_partic,boundary_con,neig,HH1,state,flag)

  HH = HH0 + mu*HH1

  do i=1,dim        
    write(25,*)HH0(:,i)
  enddo

  CALL DSYEVD('V','U',dim,HH,dim,energy,work,3*dim**2,IWORK,6*dim,info)
  CALL DSYEVD('V','U',dim,HH0,dim,energy0,work,3*dim**2,IWORK,6*dim,info)

!!!!!!!!!! PARITY SYSTEM !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1,dim        
     q = state(i)
     p = state(i)
!     write(*,*)i,state(i)
     do j=0,(sites-1)/2-1
        CALL mvbits(p,j,1,q,(sites-1)-j) 
        CALL mvbits(p,(sites-1)-j,1,q,j) 
     end do
     HH1(i,:) = HH(flag(q),:)
     HH01(i,:) = HH0(flag(q),:)
!     print*, 'flag = ',flag(q) 
     write(13,*)state(i)
     write(15,*)HH0(:,i)
     write(17,*)HH(:,i)
  enddo
  write(8,*)energy

  dim_parity = 0
  dim_parity0 = 0
  HHparity = 0.d0
  HH0parity = 0.d0
  DO i=1,dim
 !    write(*,*)i,DOT_PRODUCT(HH(:,i),HH1(:,i)),dim_parity
     IF(DOT_PRODUCT(HH(:,i),HH1(:,i)) < 0) THEN   !ODD parity for HH
        dim_parity = dim_parity +1        
        ener_vec(dim_parity) = energy(i)
        HHparity(:,dim_parity) = HH(:,i)
     END IF
  END DO

  ALLOCATE(HHparidad(dim,dim_parity))

  PRINT*,dim_parity 
  do i=1,dim_parity
   HHparidad(:,i) = HHparity(:,i)
   write(18,*)HHparity(:,i)
  end do
  
  DEALLOCATE(HHparidad)

  dim_parity = 0
  dim_parity0 = 0
  HHparity = 0.d0
  HH0parity = 0.d0
!  state_parity = 0
!  flag_parity = 0.d0
  DO i=1,dim
 !    write(*,*)i,DOT_PRODUCT(HH(:,i),HH1(:,i)),dim_parity
     IF(DOT_PRODUCT(HH(:,i),HH1(:,i)) > 0) THEN   !even parity for HH
        dim_parity = dim_parity +1        
        ener_vec(dim_parity) = energy(i)
        HHparity(:,dim_parity) = HH(:,i)
!	state_parity(dim_parity) = state(i)
!!	flag_parity(dim_parity) = flag(i)
!        write(*,*)i,state(i),dim_parity
!        write(15,*)state(i)
!	write(14,*)flag(i)
     END IF
     IF(DOT_PRODUCT(HH0(:,i),HH01(:,i)) > 0) THEN   !even parity for HH0
        dim_parity0 = dim_parity0 +1        
        ener0_vec(dim_parity0) = energy0(i)
        HH0parity(:,dim_parity0) = HH0(:,i)
     END IF

  END DO

  ALLOCATE(ener_parity(dim_parity),sf_function(dim_parity,6),ener0_parity(dim_parity))

  ener_parity = ener_vec(1:dim_parity)
  ener0_parity = ener0_vec(1:dim_parity)

  ALLOCATE(HHparidad(dim,dim_parity))

  PRINT*,dim_parity 
  do i=1,dim_parity
   HHparidad(:,i) = HHparity(:,i)
   write(5,*)HHparity(:,i)
  end do

  write(7,*)ener_parity


!!!!!!!!!! INVERSE PARTICIPATION RATIO !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*,shape(HHparidad)
  ALLOCATE(ipr(dim_parity),shaper(2))
  shaper = shape(HHparidad(1,:))
  k = shaper(1)
  shaper = shape(HHparidad(:,1))
  l = shaper(1)
  ipr = 0d0
  do i=1,k
    do j=1,l
      ipr(i) = ipr(i) + (abs(HHparidad(j,i)))**4
    enddo
    write(14,*)ipr(i)
  enddo

  CALL ITIME(T2)
  TIME_2=T2(1)*3600+T2(2)*60+T2(3)

  PRINT*,'ELAPSED TIME:',TIME_2-TIME_1,'SECONDS'

DEALLOCATE(HH,energy,work,IWORK,vec1,width,histo,vec,ener_parity,ener0_parity,ener_vec,HH1,ipr,HHparidad,shaper)

CLOSE(5)

END PROGRAM
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

