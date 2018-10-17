
  !------------------------------------------------------------------------------
  !Program to solve a Pressure Poisson equation in 1D using Finite volume       !
  !method, uses central differenicing procedure                                 !
  !Programmed on 06 April 2015                                                  !
  !------------------------------------------------------------------------------

 !---------------------------------------------------------
 !GLOBAL VARIABLES
 !==========================================================
 module mesh_data
 implicit none
 save

  integer::np,nc
  integer::ncc,npc
  integer::flag,flag1
  double precision,allocatable,dimension(:)::x,xc,volume
  double precision,allocatable,dimension(:)::x_coarse,xc_coarse,volume_coarse
 
 end module mesh_data
!==========================================================
 module solution_data
 implicit none
 save
  
  double precision,allocatable,dimension(:)::pold,pnew
  double precision,allocatable,dimension(:)::source,residual

 end module solution_data
!==========================================================
 module choice
 implicit none
 save

  integer::solver_type,max_no_of_sweeps
  integer::no_of_sweeps_ps,no_of_sweeps_c
  integer::mg_count,no_of_mg_cycles,sigma
  double precision::domain_length
  double precision::relative_residue
  double precision::starttime,endtime
  double precision::starttime1,endtime1
  double precision::starttime3,endtime3
  double precision::time_temp
  double precision::norm,norm_iter1
  double precision::norm_mg,norm_mg1
  double precision,allocatable,dimension(:) ::norm_denominator
  double precision,allocatable,dimension(:) ::norm_numerator

 end module choice
 !==========================================================
 module linear_system
 implicit none
 save

  double precision,allocatable,dimension(:,:)::Amat,Ac_mat
  double precision,allocatable,dimension(:)  ::Bvec,Solvec,Solvec_old
  double precision,allocatable,dimension(:)  ::Solvec_coarse
  double precision,allocatable,dimension(:)  ::Solvec_old_coarse

 end module linear_system
 !==========================================================
 module amg_data
 implicit none
 save

  double precision,allocatable,dimension(:,:):: restrict_mat
  double precision,allocatable,dimension(:,:):: prolong_mat
  double precision,allocatable,dimension(:)  ::error_coarse

 end module amg_data
!==========================================================


!------------------------------------------------------------------------------
!		Main program
 program ppe_solver_1d
 use mesh_data
 use solution_data
 use choice
 use linear_system
 use amg_data
 implicit none

 !----------------------
 !Local Variables

  integer::i,ios

  ! flag1 for making old term zero at first cycle
  flag1=0    

  ! temporary variable to calculate time for coarse calculation
  time_temp=0.0

  call read_choice

  call gridgen

  call source_computation
  
  allocate(norm_numerator(nc),norm_denominator(nc),stat=ios)
   if(ios /= 0) then
    print*,"Allocation failure for norm data"
    stop
   end if

  ! Initializing old solver value to zero to ensure proper convergence
  do i = 1,nc
   pold(i) = 0.0
  end do

  if(solver_type == 1) then !Only PJ, No AMG
   flag1=1           
   call construct_linear_system(nc)

   call cpu_time(starttime)
   
   call sgs(nc,max_no_of_sweeps)

   do i = 1,nc
    pnew(i) = Solvec(i)
   end do

   deallocate(Amat,Bvec,Solvec,stat=ios)
   if(ios /= 0) then
    print*,"DeAllocation failure for linear system of size:",nc
    print*,"Note that the grid size (np) is:",np
    stop
   end if

  !====================================================
  ! The other if condition in the if loop takes solver
  ! type to 2, from the given user input
  !----------------------------------------------------
  elseif(solver_type == 2) then 
    
   allocate(restrict_mat(ncc,nc),error_coarse(ncc),prolong_mat(nc,ncc),stat=ios)
   if(ios /= 0) then
    print*,"Allocation failure for amg data"
    stop
   end if

   call gridgen_coarse

   call construct_linear_system(nc)

   call restriction 

   call construct_linear_system_coarse(ncc)

   call cpu_time(starttime)

  !===============================================================
  !  Start of Multigrid loop after construction of the fine grid
  !  and other operators.
  !=============================================================== 
   mg_loop: do mg_count = 1,no_of_mg_cycles
     flag=0
    if(mg_count==1) then
     flag1=1               !old will be initialized to zero
    elseif(mg_count/=1) then
     flag1=0 
    endif

    print*,"Multicycle:", mg_count
   
    !open(32,file='presmoothsoln.dat')
    if(flag1==1) then
      call sgs(nc,no_of_sweeps_ps)     !Presmoothing
   !   do i=1,nc
   !     write(32,*) xc(i), Solvec(i)
   !   enddo
    end if
   ! close(32)
    
    !====================================================
    !  residual calculation from the  obtained pre smooth
    !  solution in original equation                     
    !====================================================
    !call cpu_time(starttime1)
    call residual_calculation
    !call cpu_time(endtime1)

    !time_temp= time_temp + (endtime1-starttime1)

    call error_coarse_calculation 

    call sgs_coarse(ncc,no_of_sweeps_c)
    
    call cpu_time(starttime1)
    call prolongation
    call cpu_time(endtime1)
    time_temp= time_temp + (endtime1-starttime1)

    flag1=0   
 
    call sgs(nc,no_of_sweeps_ps)               !postsmoothing

    deallocate(solvec_coarse,stat=ios)
    if(ios /= 0) then
     print*,"DeAllocation failure for linear system of size:",nc
     print*,"Note that the grid size (np) is:",np
     stop
    end if

    do i = 1,nc
     pnew(i) = Solvec(i)
    end do

   if(flag==1) then 
    exit
   else
    continue
   endif
   !=================================================
    norm_mg=0.0

    !do i=1,nc
    !  norm_numerator(i)  = abs(Pold(i)-Pnew(i))
    !  norm_denominator(i)= abs(Bvec(i))
    !end do

    do i=1,nc
     norm_mg= norm_mg + (pnew(i)-pold(i))**2  !maxval(norm_numerator)/maxval(norm_denominator)
    end do

    norm_mg = sqrt(norm_mg/nc)

    if(mg_count==1) norm_mg1 = norm_mg
 
    norm_mg = norm_mg/norm_mg1
    !=======================================
    write(*,*) "Norm MG:", norm_mg
    if(norm_mg<relative_residue) then
     print*,"Multicycle success !!!"  
     exit
    else
    continue
     do i=1,nc
       pold(i)=pnew(i)
     end do
    endif
    
    if(mg_count==no_of_mg_cycles) print*,"Input more mg_cycles"
   end do mg_loop
   !call cpu_time(endtime1)
   !write(*,*) "AMG Cpu time is:", endtime1-starttime
  else

   print*,"Wrong choice of solver type:: Re-check"

  end if

  call cpu_time(endtime)

  write(*,*) "Total  Cpu time is:", endtime-starttime
  write(*,*) "Extra time is :", time_temp
  call write_output

 end program ppe_solver_1d

!------------------------------------------------------------------------------
 subroutine read_choice
 use mesh_data
 use solution_data
 use choice
 implicit none

 !Local Variables
 !---------------
  integer::ios

  open(unit=1,file='user_choice.inp',iostat=ios)
  if(ios /= 0) then
   print*,"File user_choice.inp can not be found"
   stop
  end if

  read(1,*)solver_type
  read(1,*)domain_length,np
  read(1,*)max_no_of_sweeps,no_of_sweeps_ps,no_of_sweeps_c
  read(1,*)relative_residue
  read(1,*)no_of_mg_cycles
  read(1,*)sigma
  print*," "

  if(solver_type==1) then
    print*, "The User Choice is: Only SGS"
  elseif(solver_type==2) then
    print*, "The User Choice is: Algebraic Multigrid with SGS"
  endif
  close(1)

  nc  = np - 1

  ncc = nc / 2

  npc = ncc + 1 

 !Allocation for mesh data
 !------------------------

  allocate(x(np),xc(nc),volume(nc),stat=ios)
  if(ios /= 0) then
   print*,"Allocation failure for mesh data"
   stop
  end if

!Allocation for Coarse mesh data
 !------------------------

  allocate(x_coarse(npc),xc_coarse(ncc),volume_coarse(ncc),stat=ios)
  if(ios /= 0) then
   print*,"Allocation failure for mesh data"
   stop
  end if

 !Allocation for solution data
 !----------------------------

  allocate(pnew(nc),pold(nc),source(nc),residual(nc),stat=ios)
  if(ios /= 0) then
   print*,"Allocation failure for solution data"
   stop
  end if

 end subroutine read_choice

!------------------------------------------------------
 subroutine gridgen
 use mesh_data
 use choice
 implicit none

 !Local Variables
 !---------------
  integer::i
  double precision::dx

  dx = domain_length / dble(nc)

  do i = 1,np
   x(i) = dx * dble(i-1)
  end do

  do i = 1,nc
 
   xc(i) = 0.5d0 * (x(i+1) + x(i))
 
   volume(i) = x(i+1) - x(i)
 
  end do

!  open(unit=2,file='point.dat')
!  do i = 1,np
!   write(2,*)i,x(i),0.0d0
!  end do
!
!  close(2)

!  open(unit=3,file='cent_vol.dat')
!
!  do i = 1,nc
!   write(3,*)i,xc(i),0.0d0,volume(i)
!  end do
!
!  close(3)

 end subroutine gridgen

!------------------------------------------------------------------------------
 subroutine source_computation
 use mesh_data
 use choice
 use solution_data
 implicit none

 !Local Variables
 !---------------
  integer::i
  double precision::pi

  pi = 4.0d0 * datan(1.0d0)

  do i = 1,nc

  !  source(i) = 0.0
  !source(i) = -4.0d0 * pi**2 * dsin(2.0d0 * pi * xc(i))
  source(i) = (-4.0d0 * pi**2 * dsin(2.0d0 * pi * xc(i)))/(sigma+ 4.d0*pi**2)
  end do

 ! source(nc) = 100.0
  
  open(unit=4,file='source.dat')

  !do i = 1,nc
  ! write(4,*) i, xc(i), source(i)
  !end do

  close(4)

 end subroutine source_computation

!------------------------------------------------------------------------------
subroutine construct_linear_system(lsize)
 use mesh_data
 use solution_data
 use linear_system
 use choice
 implicit none

 integer,intent(in)::lsize

 !Local Variables
 !---------------
  integer::i,ssize,ii,ios,j
  double precision::del_plus,del_minus,pi

  pi = 4.0d0 * datan(1.0d0)

  ssize = lsize-2

  allocate(Amat(ssize,3),Bvec(ssize),Solvec(lsize),stat=ios)
  if(ios /= 0) then
   print*,"Allocation failure for linear system of size:",lsize
   print*,"line 338"
   print*,"Note that the grid size (np) is:",np
   stop
  end if

 allocate(Solvec_old(lsize),stat=ios)
  if(ios /= 0) then
   print*,"Allocation failure for linear system of size:",lsize
   print*,"Note that the grid size (np) is:",ssize
   stop
  end if

  do i = 1,ssize

   ii = i + 1

   del_plus = x(ii+1) - x(ii)
   del_minus = x(ii) - x(ii-1)

   Amat(i,1) = 1.0d0 / del_minus 
   Amat(i,2) = -( (1.0d0 / del_minus) + (1.0d0 / del_plus) ) 
   Amat(i,3) = 1.0d0 / del_plus 

   Bvec(i) = source(ii) * volume(ii)

  end do

  !open(unit=5,file='Amat_Bvec.dat')
  !do i = 1,ssize
  ! write(5,*)(Amat(i,j),j=1,3),Bvec(i)
  !end do
  !close(5)

  do i=1,nc
   Solvec(i) = 0.0d0
  end do

  !Solvec(1) = dsin(2.0d0*pi*xc(1))
  !Solvec(nc) = dsin(2.0d0*pi*xc(nc))

  Solvec(1) = dsin(2.0d0*pi*xc(1)) /(sigma+ 4.0d0*pi**2)
  Solvec(nc)= dsin(2.0d0*pi*xc(nc))/(sigma+ 4.0d0*pi**2)


 end subroutine construct_linear_system

!------------------------------------------------------------------------------
subroutine sgs(lsize,nsweeps)
 use mesh_data
 use solution_data
 use choice
 use linear_system
 use choice
 implicit none

 integer,intent(in)::lsize,nsweeps

 !Local Variables
 !---------------
  integer::i,ssize,istart,iend,incr,ii,count,ios
  double precision::rhs,pi
  flag=0

  pi = 4.0d0 * datan(1.0d0)

  ssize = lsize - 2
 !--------------------------------------
 !		Value used for checking norm
 if(flag1==1) then
  do i= 2,ssize+1
   Solvec_old(i) = 0.0d0
  end do
 elseif(flag1==0) then
  do i= 2,ssize+1
   Solvec_old(i) = Solvec(i) 
  end do
 endif
  Solvec_old(1) = Solvec(1)
  Solvec_old(lsize) = Solvec(lsize)
!---------------------------------------
!		Assigning Forward or reverse sweep
  sweeps_loop: do count = 1,(2 * nsweeps)

   if(mod(count,2) /= 0) then !Forward sweep

    istart = 2
    iend = ssize + 1 
    incr = 1

   elseif(mod(count,2) == 0) then !Reverse sweep

    istart = ssize + 1
    iend = 2
    incr = -1

   end if

!-------------------------------------
!			SGS calculation  
   cell_loop: do i = istart,iend,incr

    ii = i - 1

    rhs = (Solvec_old(i-1) * Amat(ii,1)) + (Solvec_old(i+1) * Amat(ii,3)) !PJ uses old value
    !rhs = (Solvec(i-1) * Amat(ii,1)) + (Solvec(i+1) * Amat(ii,3))        !SGS uses latest value

    Solvec(i) = (Bvec(ii) - rhs) / Amat(ii,2)

   end do cell_loop

 !---------------------------------------------
 !		Checking for convergence after reverse sweep
   !=============================================
   !   L infinity norm
    norm=0.0

    do i=istart,iend,incr
      norm_numerator(i)  = abs(Solvec_old(i)-Solvec(i))
      norm_denominator(i)= abs(Bvec(i))
    end do

    norm= maxval(norm_numerator)/maxval(norm_denominator)
 
 if(solver_type==1 .and. mod(count,2) == 0) then
    write(*,*) "Sweeps", count, "Norm:", norm
    if(norm<=relative_residue) then
     print*,"PJ success !!!"  
     exit
    else
     continue
    endif
    !============================================
 elseif(solver_type==2 .and. mod(count,2) == 0) then
    write(*,*) "Sweeps on fine grid", count, "Norm:", norm
    if(norm<=relative_residue) then
     print*,"AMG success !!!"  
     flag=1
     exit
    else
      continue
    endif
  !============================================
  ! L2 Norm
  ! norm=0.0
  ! do i=istart,iend,incr 
  !   norm= norm+(Solvec(i)-Solvec_old(i))**2
  ! end do

  ! norm=sqrt(norm/abs(istart-iend))

	!if(count==2) norm_iter1 = norm

  !  norm = norm/norm_iter1
  !  !open(7,file="Sgs_norm.dat")
  !   write(*,*) "sweeps:",count, "		Norm:",norm
  !   !write(7,*) count,norm
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  if(norm<=relative_residue) then
  !   write(*,*) "Converged " !"iterations:",iter, "		Norm:",norm
  !   exit
  !  else
  !   !If(count==2*nsweeps) print*,"Given sweeps completed"
  !   continue
  !  endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-= 
  !elseif(solver_type==2 .and. mod(count,2) == 0) then
  !  write(*,*) "Sweeps on fine grid", count, "Norm:", norm
  !  if(norm<=relative_residue) then
  !   print*,"AMG success !!!"  
  !   flag=1
  !   exit
  !  else
  !    continue
  !  endif
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=

endif
  
  do i=2,ssize+1
    Solvec_old(i) = Solvec(i)
  end do

  end do sweeps_loop

   
 end subroutine sgs

!------------------------------------------------------------------------------
 subroutine write_output
 use mesh_data
 use solution_data
 use choice
 implicit none

 !Local Variables
 !---------------
  integer::i
  double precision::pi

  pi = 4.0d0 * datan(1.0d0)

  open(unit=6,file='solution.dat')

  do i = 1,nc

   write(6,*) xc(i),pnew(i), dsin(2.0d0 * pi * xc(i))/(sigma+ 4.0d0*pi**2)

  end do

  close(6)

 end subroutine write_output

!------------------------------------------------------------------------------
subroutine gridgen_coarse
 use mesh_data
 use choice
 use linear_system
 implicit none

 !Local Variables
 !---------------
  integer::i,ios
  double precision::dx

  dx = domain_length / dble(ncc)

  do i = 1,npc
   x_coarse(i) = dx * dble(i-1)
  end do

  do i = 1,ncc
 
   xc_coarse(i) = 0.5d0 * (x_coarse(i+1) + x_coarse(i))
 
   volume_coarse(i) = x_coarse(i+1) - x_coarse(i)
 
  end do

  !open(unit=11,file='coarse_point.dat')

  !do i = 1,npc
  ! write(11,*)i,x_coarse(i),0.0d0
  !end do

  !close(11)

  !open(unit=12,file='coarse_cent_vol.dat')

  !do i = 1,ncc
  ! write(12,*)i,xc_coarse(i),0.0d0,volume_coarse(i)
  !end do

  !close(12)

  allocate(Solvec_old_coarse(ncc),stat=ios)
  if(ios /= 0) then
   print*,"Allocation failure for coarse mesh data"
   stop
  end if


 end subroutine gridgen_coarse
!-----------------------------------------------

!---------------------------------------------------------
subroutine residual_calculation
 use mesh_data
 use solution_data
 use linear_system
 implicit none

 integer i
 do i=2,nc-1
  residual(i) =  Bvec(i-1) - (Amat(i-1,1)* Solvec(i-1) + Amat(i-1,2)* Solvec(i) + Amat(i-1,3)* Solvec(i+1))  !- Bvec(i-1)
 end do

 residual(1) = 0.0d0        !	Exact values at boundaries
 residual(nc)= 0.0d0        !

 open(14,file="residual.dat")
!
 do i=1,nc
 write(14,*)residual(i)
 end do 
!
 close(14)

end subroutine residual_calculation

!---------------------------------------------------------
subroutine restriction 
 use amg_data
 use mesh_data
 use solution_data
 Implicit none

 Integer i,j,k1,k2

 do i=1,ncc
  do j=1,nc
   restrict_mat(i,j) = 0.0d0
  end do
 end do

 do i=1,ncc
  k2 = 2*i
  k1 = 2*i-1
  do j=1,nc
   If(j==k1) then
    restrict_mat(i,j)= 0.5
   elseif(j==k2) then
    restrict_mat(i,j)= 0.5
   else
    restrict_mat(i,j)= 0.0d0
   end if
  end do
 end do

!Open(15,file="restrict_mat.dat")
!
!do i=1,ncc
!write(15,*) (restrict_mat(i,j),j=1,nc)
!end do
!close(15)

Close(15)


!--------------------------------
!	Prolongation operator
prolong_mat=0.0

 do i=1,nc
  do j=1, ncc
   k1=2*j-1
   k2=2*j
   if(i==k1 .or. i==k2) then
    prolong_mat(i,j)= 1.0
   else
    prolong_mat(i,j)= 0.0d0
   endif
  end do
 end do

!open(17,file="prolong_mat.dat")
!
!do i=1,nc
!write(17,*) (prolong_mat(i,j),j=1,ncc)
!enddo
! close(17)

end subroutine restriction


!-----------------------------------------------------
!  Calculate error on coarse 
subroutine error_coarse_calculation
 use mesh_data
 use linear_system
 use amg_data
 use solution_data
 implicit none

Integer i,j

 do i=1,ncc
  error_coarse(i) = 0.0
 end do


 do i=1,ncc
  
  do j=1,nc
 
   error_coarse(i)= error_coarse(i)+ (restrict_mat(i,j) * residual(j))

  enddo

 end do

 open(16,file="error_coarse.dat")
 do i=1,ncc
  write(16,*) error_coarse(i)
 end do

 Close(16)
end subroutine  error_coarse_calculation
!-----------------------------------------------------
!	Coarse matrix construction
subroutine construct_linear_system_coarse(lsize)
 use mesh_data
 use linear_system
 use amg_data
 implicit none

 integer i,j,k,k1,k2,ios,ll,lsize
 double precision atemp
 double precision, dimension(:)  ::AA(nc)
 double precision, dimension(:,:)::M_mat(ncc,nc)
 double precision, dimension(:,:)::A_temp(nc,nc),Ac_temp(ncc,ncc)

 !ssize= lsize/2
 allocate(Ac_mat(ncc,3),stat=ios)
 if(ios /= 0) then
   print*,"Allocation failure for coarse mesh data"
   stop
 end if
  
 Ac_mat = 0.0d0
 M_mat  = 0.0d0
 A_temp = 0.0d0

 A_temp(1,1) = Amat(1,2)
 A_temp(1,2) = Amat(1,3)
 A_temp(nc,nc-1) = Amat(nc-2,1)
 A_temp(nc,nc) = Amat(nc-2,2)

 do i = 2,nc-1

  do j = 1,nc

   if(j == i-1) then
 
    A_temp(i,j) = Amat(i-1,1)
 
   elseif(j == i+1) then
 
    A_temp(i,j) = Amat(i-1,3)
 
  elseif(j == i) then

    A_temp(i,j) = Amat(i-1,2)

   end if

  end do

 end do


! Open(25,file="A_temp.dat")
! do i=1,nc
!  write(25,*) (A_temp(i,j),j=1,nc)
! end do
! close(25)

 do i=1,ncc
  do j=1,nc
   do k=1,nc
    M_mat(i,j)= M_mat(i,j) + restrict_mat(i,k) * A_temp(k,j)
   end do
  end do
 end do

! open(22,file="M_mat.dat")
! do i=1,ncc
!  write(22,*) (M_mat(i,j),j=1,nc)
! end do
! close(22)

 Ac_temp = 0.0d0

 do i = 1,ncc

  do j = 1,ncc

   do k = 1,nc

    Ac_temp(i,j) = Ac_temp(i,j) + M_mat(i,k) * prolong_mat(k,j)

   end do

  end do

 end do

 
! open(20,file="Ac_temp.dat")
! do i = 1,ncc
!  write(20,*)(Ac_temp(i,j),j=1,ncc)
! end do
! close(20)

 Ac_mat = 0.0d0

 do i = 2,ncc-1

  do j = 1,ncc

   if(j == i-1) then
    Ac_mat(i,1) = Ac_temp(i,j) 
   elseif(j == i+1) then
    Ac_mat(i,3) = Ac_temp(i,j) 
   elseif(j == i) then
    Ac_mat(i,2) = Ac_temp(i,j) 
   end if

  end do

end do


 Ac_mat(1,1) = Ac_temp(1,2)
 Ac_mat(1,2) = Ac_temp(1,1)
 Ac_mat(1,3) = Ac_temp(1,2)
 Ac_mat(ncc,1) = Ac_temp(ncc,ncc-1)
 Ac_mat(ncc,2) = Ac_temp(ncc,ncc)
 Ac_mat(ncc,3) = Ac_temp(ncc,ncc-1)

! open(26,file="Ac_mat.dat")
! do i = 1,ncc
!  write(26,*)(Ac_mat(i,j),j=1,3)
! end do
! close(26)

end subroutine construct_linear_system_coarse

!-------------------------------------------------
subroutine sgs_coarse(lsize,nsweeps)
 use mesh_data
 use solution_data
 use choice
 use linear_system
 use amg_data
 implicit none

 !Local Variables
 !---------------
  integer,intent(in)::lsize,nsweeps
  integer::i,ssize,istart,iend,incr,ii,count,ios
  double precision::rhs,pi

  

  allocate(Solvec_coarse(lsize),stat=ios)
  if(ios /= 0) then
   print*,"Allocation failure for coarse mesh data 1"
   stop
  end if
  
  pi = 4.0d0 * datan(1.0d0)

!--------------------------------------
!		Value used for checking norm
 if(flag1==1) then
  print*,"1st cycle Solvec_old_coarse=0.0"
  do i= 1,lsize
   Solvec_old_coarse(i) = 0.0
  end do
 end if

 !print*," "
 ! if(solver_type==2) print*,"On coarse:"
 !----------------------------------------


 !---------------------------------------
 !		Assigning Forward or reverse sweep
  sweeps_loop: do count = 1,(2 * nsweeps)
 !print*,"SGS Coarse count", count
   if(mod(count,2) /= 0) then !Forward sweep

    istart = 1
    iend = lsize 
    incr = 1

   elseif(mod(count,2) == 0) then !Reverse sweep

    istart = lsize
    iend = 1
    incr = -1

   endif

   !-------------------------------------
   !	Coarse	SGS calculation  
   cell_loop: do i = istart,iend,incr
 
    if(i == istart) then

     if(incr == 1) then
      
      rhs = Solvec_old_coarse(i+1) * Ac_mat(i,3) 

     elseif(incr == -1) then

      rhs = Solvec_old_coarse(i-1) * Ac_mat(i,1)

     end if

    elseif(i == iend) then

     if(incr == 1) then
      
      rhs = Solvec_old_coarse(i-1) * Ac_mat(i,1)

     elseif(incr == -1) then

      rhs = Solvec_old_coarse(i+1) * Ac_mat(i,3) 

     end if

    else

     rhs = (Solvec_old_coarse(i-1) * Ac_mat(i,1)) + (Solvec_old_coarse(i+1) * Ac_mat(i,3)) 

    end if

    Solvec_coarse(i) = (error_coarse(i) - rhs) / Ac_mat(i,2)

   end do cell_loop

 !-----------------------------------------------
 ! It is solver type 2, but just giving condition
 !  and also for the reverse sweep only
  if(solver_type==2 .and. mod(count,2) == 0) then

   norm=0.0
 !   do i=istart,iend,incr
 !     norm= norm+(Solvec_coarse(i)-Solvec_old_coarse(i))**2
 !   end do

 !  norm=sqrt(norm/abs(istart-iend))

 !  if(count==2) norm_iter1 = norm

 !  norm = norm/norm_iter1

 !=============================================
   !   L infinity norm

   do i=istart,iend,incr
     norm_numerator(i)  = abs(Solvec_old_coarse(i)-Solvec_old(i))
     norm_denominator(i)= abs(error_coarse(i))
   end do

   norm= maxval(norm_numerator)  !/maxval(norm_denominator)
 
   
  ! write(*,*) "Sweeps", count, "Norm:", norm
   if(norm<=relative_residue) then
   write(*,*) "Converged on coarse" ,count !, "		Norm:",norm
   print*, ""
    exit
   else
   ! If(count==2*nsweeps) print*,"Given Sweeps completed on Coarse, Not Converged"
    continue
   endif

 endif

  do i=1,lsize
    Solvec_old_coarse(i) = Solvec_coarse(i)
  end do

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=

 end do sweeps_loop

   open(22,file="Solvec_coarse.dat")
   do i=1,lsize
    write(22,*) xc(i), Solvec_coarse(i)
   end do
   close(22)

end subroutine sgs_coarse

!-------------------------------------------------
subroutine prolongation
 use amg_data
 use mesh_data
 use solution_data
 use linear_system
 implicit none

 integer i,j,k1

 do i=2,nc-1
  do j=1,ncc
    Solvec(i) = Solvec(i) + prolong_mat(i,j) * solvec_coarse(j) 
  end do
 end do

 open(31,file="prolongsoln.dat")
 do i= 1,nc
  write(31,*) xc(i), Solvec(i)
 end do
 Close(31)
end subroutine prolongation
!--------------------------------------------------------
 !=================================================
!subroutine convergence
! use choice

! integer i
   ! norm_mg=0.0
    !do i=1,nc
    !  norm_numerator(i)  = abs(Pold(i)-Pnew(i))
    !  norm_denominator(i)= abs(Bvec(i))
    !end do

    !do i=1,nc
    !norm_mg= maxval(norm_numerator)/maxval(norm_denominator)
    !end do

    !norm_mg = sqrt(norm_mg/nc)

    !if(mg_count==1) norm_mg1 = norm_mg
 
    !norm_mg = norm_mg/norm_mg1
    !=======================================
    !write(*,*) "Norm MG:", norm_mg
    !if(norm_mg<relative_residue) then
    ! print*,"Multicycle success !!!"  
    ! exit
    !else
    !continue
    ! do i=1,nc
    !   pold(i)=pnew(i)
    ! end do
    !endif
!end subroutine convergence    
!*******************************************************************
!				              	-xXxXxX-
