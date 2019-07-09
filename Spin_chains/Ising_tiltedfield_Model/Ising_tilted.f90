MODULE global
  
  COMPLEX(8), SAVE :: ii=(0.d0,1.d0)
  REAL(8), SAVE :: pi=4.d0*atan2(1.d0,1.d0)     
  
END MODULE global


PROGRAM ISING

  USE HAMILT
  USE GLOBAL

  IMPLICIT NONE

  ! PARA LAS RUTINAS LAPACK USAR COMPLEX(8) !!!!!!!!!!!!!!!!!!!

  INTEGER :: i,t,j,k,l,p,q,sites,dim,info,bins,n_cut,neig,n_partic,dim_parity,dim_parity0,mean_value,jj,label_state,r
  INTEGER :: n_state,n_average,T1(3),T2(3),TIME_1,TIME_2,opera
  REAL(8) :: GG,alpha_x,alpha_y,alpha_z,alpha2_x,alpha2_y,alpha2_z,a(4),mu,ele,g,lau
  REAL(8) :: alpha3_x,alpha3_y,alpha3_z,tita,pie
  CHARACTER(20) ::  boundary_con

  INTEGER, ALLOCATABLE, SAVE :: IWORK(:),state(:),flag(:),statepp(:),statepn(:),coupl(:)
  REAL(8), ALLOCATABLE, SAVE :: width(:),histo(:),ener_parity(:),ener_vec(:),energy(:),sf_function(:,:),op(:,:),ener_parity_odd(:)
  REAL(8), ALLOCATABLE, SAVE :: HH(:,:),HH0parity(:,:),HH1(:,:),HHparity(:,:),HH0(:,:),HH01(:,:),overlaps(:),HH2(:,:),HH3(:,:)
  REAL(8), ALLOCATABLE, SAVE :: energy0(:),work(:), vec1(:),vec(:),ener0_vec(:),ener0_parity(:),Sdiag(:),ener0_parity_odd(:)
  REAL(8), ALLOCATABLE, SAVE :: HHZ(:,:), shaper(:),ipr(:),HHparidad(:,:)

  OPEN(8,FILE = 'datos.in',status='unknown')

  read(8,*)sites
  read(8,*)boundary_con
  read(8,*)alpha_x,alpha_y,alpha_z
  read(8,*)mu,tita
  pie = 2.D0*ATAN(1.D0)
  dim = 2**sites 
  !dim = int(exp(logfac(sites)-logfac(n_partic)-logfac(sites-n_partic)))
  PRINT*,'dim= ',dim
  PRINT*,'theta= ', tita*pie

!  OPEN(5,FILE = 'hamilt.out',status='unknown')
  OPEN(7,FILE = 'energiasfull.out',status='unknown')
!  OPEN(8,FILE = 'h0.out',status='unknown')
!  OPEN(11,FILE = 'hh.out',status='unknown')
!  OPEN(12,FILE = 'hh_parity.out',status='unknown')
!  open(15,file = 'statesparity_even.out', status = 'unknown')
!  open(15,file = 'statesparity_odd', status = 'unknown')
  OPEN(51,FILE = 'HHparity_even.out',status='unknown')
!  OPEN(52,FILE = 'HHparity_odd.out',status='unknown')
  OPEN(71,FILE = 'ener_parity_even.out',status='unknown')
  OPEN(14,FILE = 'ipr.out',status='unknown')
  OPEN(15,FILE = 'hamilt0full.out', status='unknown')
  open(17,file = 'hamiltfull.out', status = 'unknown')
!  OPEN(72,FILE = 'ener_parity_odd.out',status='unknown')

  ALLOCATE(HH(dim,dim),energy(dim),vec1(dim),width(bins),HH1(dim,dim),vec(dim),histo(bins),HH0(dim,dim),HH2(dim,dim),HH3(dim,dim))
  ALLOCATE(work(3*dim**2),IWORK(6*dim),state(dim),flag(2**sites),ener_vec(dim),HHparity(dim,dim),statepp(dim),statepn(dim))
  ALLOCATE(HH0parity(dim,dim),energy0(dim),HH01(dim,dim),ener0_vec(dim),HHZ(dim,dim))

  CALL ITIME(T1)
  TIME_1=T1(1)*3600+T1(2)*60+T1(3)
  
   GG = 0.D0!1.5d0 !0.

  PRINT*, alpha_x,alpha_y,alpha_z
  neig = 1
  CALL ising_chain(GG,alpha_x,alpha_y,alpha_z,HH0,sites,boundary_con,neig)
  
!  neig = 2
!  CALL ising_chain(GG,alpha2_x,alpha2_y,alpha2_z,HH2,sites,boundary_con,neig)

!  neig = 3
!  CALL ising_chain(GG,alpha3_x,alpha3_y,alpha3_z,HH3,sites,boundary_con,neig)     

! campo trasverso en x
   HH1(:,:)=0d0 
   do i=1,dim
    do j=0,sites-1
!	print*, I,j,btest(I-1,j),ibset(i-1,j),ibclr(i-1,j)
    if(btest(i-1,j).eqv..true.)jj=ibclr(i-1,j)+1
    if(btest(i-1,j).eqv..false.)jj=ibset(i-1,j)+1
    HH1(i,jj)=HH1(i,jj)+1d0
!    HH1(jj,i)=HH1(i,jj)
    end do
   end do
!  do i=1,dim
!   write(55,*)HH1(i,:)
!  end do

! campo trasverso en z
   HHZ(:,:)=0d0 
   do i=1,dim
    do j=0,sites-1
!	print*, I,j,btest(I-1,j),ibset(i-1,j),ibclr(i-1,j)
    if(btest(i-1,j).eqv..true.)then
     jj=ibset(i-1,j)+1
     HHZ(i,jj) = HHZ(i,jj) - 1d0
    endif
    if(btest(i-1,j).eqv..false.)then
     jj=ibclr(i-1,j)+1
     HHZ(i,jj) = HHZ(i,jj) + 1d0
!    HH1(jj,i)=HH1(i,jj)
    endif
    end do
   end do
!  do i=1,dim
!   write(55,*)HH1(i,:)
!  end do


!  print*, "pi =", pie
  HH = HH0 + mu*(SIN(0.5*tita*pi)*HH1 + COS(0.5*tita*pi)*HHZ)

!  do i=1,dim
!   write(8,*)HH0(i,:)
!  end do
!  do i=1,dim
!   write(11,*)HH(i,:)
!  end do

  CALL DSYEVD('V','U',dim,HH,dim,energy,work,3*dim**2,IWORK,6*dim,info)
  CALL DSYEVD('V','U',dim,HH0,dim,energy0,work,3*dim**2,IWORK,6*dim,info)

  PRINT*,'dim= ',dim

  do i=1,dim
   write(15,*)HH0(:,i)
   write(17,*)HH(:,i)
  end do

  write(7,*)energy

!  goto 999

!! LOS CONTADORES FUNCIONAN AL REVES! 0 = ultimo spin, N = spin 0 y asi!!!


! PARITY SYSTEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
        if (MOD(sites,2).eq.0) then 
            do i=0,dim-1        
              q = i
              p = i
!              write(*,*)i,state(i)
              do j=0,sites/2-1
                CALL mvbits(p,j,1,q,(sites-1)-j)
                CALL mvbits(p,(sites-1)-j,1,q,j)
              end do
            HH1(i+1,:) = HH(q+1,:)
            HH01(i+1,:) = HH0(q+1,:)
            enddo
            else
            do i=0,dim-1        
              q = i
              p = i
!              write(*,*)i,state(i)
              do j=0,(sites-1)/2-1
                CALL mvbits(p,j,1,q,(sites-1)-j)
                CALL mvbits(p,(sites-1)-j,1,q,j)
              end do
            HH1(i+1,:) = HH(q+1,:)
            HH01(i+1,:) = HH0(q+1,:)
            enddo
        end if

!  do i=1,dim
!   write(12,*)HH1(:,i)
!  end do

!even parity

  dim_parity = 0
  dim_parity0 = 0
  HHparity = 0.d0
  HH0parity = 0.d0

  DO i=1,dim
      !write(*,*)i,DOT_PRODUCT(HH(:,i),HH1(:,i)),dim_parity
     IF(DOT_PRODUCT(HH(:,i),HH1(:,i)) > 0) THEN   !even parity for HH
        dim_parity = dim_parity +1        
        ener_vec(dim_parity) = energy(i)
        HHparity(:,dim_parity) = HH(:,i)
!        write(15,*)state(i)
!        write(*,*)i,state(i),dim_parity
     END IF
     IF(DOT_PRODUCT(HH0(:,i),HH01(:,i)) > 0) THEN   !even parity for HH0
        dim_parity0 = dim_parity0 +1        
        ener0_vec(dim_parity0) = energy0(i)
        HH0parity(:,dim_parity0) = HH0(:,i)
     END IF
  END DO

  PRINT*,'even parity dim = ',dim_parity 
  do i=1,dim_parity
   write(51,*)HHparity(:,i)
  end do

  ALLOCATE(ener_parity(dim_parity),sf_function(dim_parity,6),ener0_parity(dim_parity))

  ener_parity = ener_vec(1:dim_parity)
  ener0_parity = ener0_vec(1:dim_parity)

  write(71,*)ener_parity

!odd parity

!  dim_parity = 0
!  dim_parity0 = 0
!  HHparity = 0.d0
!  HH0parity = 0.d0

!  DO i=1,dim
!      write(*,*)i,DOT_PRODUCT(HH(:,i),HH1(:,i)),dim_parity
!     IF(DOT_PRODUCT(HH(:,i),HH1(:,i)) < 0) THEN   !even parity for HH
!        dim_parity = dim_parity +1        
!        ener_vec(dim_parity) = energy(i)
!        HHparity(:,dim_parity) = HH(:,i)
!        write(16,*)state(i)
!        write(*,*)i,state(i),dim_parity
!     END IF
!     IF(DOT_PRODUCT(HH0(:,i),HH01(:,i)) < 0) THEN   !even parity for HH0
!        dim_parity0 = dim_parity0 +1        
!        ener0_vec(dim_parity0) = energy0(i)
!        HH0parity(:,dim_parity0) = HH0(:,i)
!     END IF
!  END DO

!  PRINT*,'odd parity dim = ',dim_parity 
!  do i=1,dim_parity
!   write(52,*)HHparity(:,i)
!  end do

!  ALLOCATE(ener_parity_odd(dim_parity),ener0_parity_odd(dim_parity))

!  ener_parity_odd = ener_vec(1:dim_parity)
!  ener0_parity_odd = ener0_vec(1:dim_parity)

 ! write(72,*)ener_parity_odd
 

! END PARITY SYSTEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Operator S_z_i in parity basis did by Diego

  ALLOCATE(op(dim_parity,dim_parity))

!  write(*,*)opera
!  do i=1,dim_parity
!   do j=1,dim_parity

!    ele=0d0
!    do k=1,dim
!     if(ibits(state(k),opera-1,1).eq.0)ele=ele-HHparity(k,i)*HHparity(k,j)
!!     if(ibits(state(k),opera-1,1).eq.1)ele=ele+HHparity(k,i)*HHparity(k,j)
 !   end do

!    op(i,j)=ele
!   end do
!  end do

!  opera=opera+1

!  do i=1,dim_parity
!   do j=1,dim_parity

!    ele=0d0
!    do k=1,dim
!!     ele=ele+HHparity(k,i)*HHparity(k,j)
!     if(ibits(state(k),opera,1).eq.0)ele=ele-HHparity(k,i)*HHparity(k,j)
!     if(ibits(state(k),opera,1).eq.1)ele=ele+HHparity(k,i)*HHparity(k,j)
!     write(*,*)i,j,k,ele,state(k),ibits(state(k),opera,1)
!    end do

!    op(i,j)=ele+op(i,j)
!!    write(*,*)i,j,op(i,j)
!   end do
!  end do


!  do i=1,dim_parity
!   write(9,*)op(i,:)
!  end do
  ALLOCATE(HHparidad(dim,dim_parity))

  PRINT*,dim_parity 
  do i=1,dim_parity
   HHparidad(:,i) = HHparity(:,i)
!   write(5,*)HHparity(:,i)
  end do

  print*,shape(HHparidad)
  ALLOCATE(ipr(dim_parity),shaper(2))
  shaper = shape(HHparidad(1,:))
  k = shaper(1)
  shaper = shape(HHparidad(:,1))
  l = shaper(1)
  ipr = 0d0
  do i=1,k
    do j=1,l
      ipr(i) = ipr(i) + (abs(HHparity(j,i)))**4
    enddo
    write(14,*)ipr(i)
  enddo


  PRINT*,'dim= ',dim

  CALL ITIME(T2)
  TIME_2=T2(1)*3600+T2(2)*60+T2(3)

  PRINT*,'ELAPSED TIME:',TIME_2-TIME_1,'SECONDS'

!DEALLOCATE(HH,energy,work,IWORK,vec1,width,histo,vec,ener_parity,ener0_parity,ener_vec,HH1,op)
DEALLOCATE(HH,energy,work,IWORK,vec1,width,histo,ener_vec,HH1,HH2,HH3,HHZ,HHparidad,shaper,ipr)


CLOSE(5)

END PROGRAM
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

