module set_precision
  use iso_fortran_env, only : real64
  implicit none
  private
  public :: dp
  integer, parameter :: dp  = real64
end module set_precision

module set_constants
  use set_precision, only : dp
  implicit none
  private
  public :: zero, one, two, three, four
  public :: half, third, fourth, large, pi
  real(dp), parameter :: zero      = 0.0_dp
  real(dp), parameter :: one       = 1.0_dp
  real(dp), parameter :: two       = 2.0_dp
  real(dp), parameter :: three     = 3.0_dp
  real(dp), parameter :: four      = 4.0_dp
  real(dp), parameter :: third     = one / three
  real(dp), parameter :: fourth    = 0.25_dp
  real(dp), parameter :: half      = 0.50_dp
  real(dp), parameter :: large  = huge(one)
  real(dp), parameter :: pi     = acos(-one)
end module set_constants

module project_inputs
  implicit none
  private
  public :: verbose_level
  integer :: verbose_level = 0
end module project_inputs

module message

  implicit none

  private

  public :: error_message, warning_message
  public :: WARN_ALWAYS, WARN_SOMETIMES, WARN_RARELY

  integer, parameter :: WARN_ALWAYS    = 0
  integer, parameter :: WARN_SOMETIMES = 1
  integer, parameter :: WARN_RARELY    = 2

contains

!================================ error_message ==============================80
!>
!! Description: Writes an error message to the screen.
!!
!! Inputs:      routine_name: Routine in which error is occuring
!!              message:      Error message to print to the screen
!!
!! Outputs:     err:          Flag indicating an error
!<
!=============================================================================80
  function error_message( routine_name, message ) result( err )

    use ISO_FORTRAN_ENV, only : error_unit

    implicit none

    character(*), intent(in) :: routine_name
    character(*), intent(in) :: message
    logical                  :: err

    continue

    err = .true.

    write(error_unit,*)
    write(error_unit,*) ' ERROR: In ' // trim(routine_name)
    write(error_unit,*) '   ', trim(message)
    write(error_unit,*) ' Stopping ...'
    call abort
    stop

  end function error_message

!=============================== warning_message =============================80
!>
!! Description: Writes a warning message to the screen.
!!
!! Inputs:      warn_level:   Important level for warning output
!!              routine_name: Routine in which warning is occuring
!!              message:      Warning message to print to the screen
!!
!! Outputs:     warn:         Flag indicating an warning
!<
!=============================================================================80
  function warning_message( warn_level, routine_name, message ) result( warn )

    use ISO_FORTRAN_ENV, only : error_unit
    use project_inputs,  only : verbose_level

    implicit none

    integer,      intent(in) :: warn_level
    character(*), intent(in) :: routine_name
    character(*), intent(in) :: message
    logical                  :: warn

    continue

    ! Setup
    warn = .true.

    ! Print Warning Message
    if ( warn_level <= verbose_level ) then
      write(error_unit,*)
      write(error_unit,*) ' WARNING: In ' // trim(routine_name)
      write(error_unit,*) '   ', trim(message)
    end if

  end function warning_message

end module message

module timer_derived_type

  use set_precision, only : dp
  use set_constants, only : zero

  implicit none

  private

  public :: basic_timer_t

  type :: basic_timer_t
    private
    real(dp)         :: time_start   = zero
    real(dp), public :: time_elapsed = zero
  contains
    private
    procedure, public, pass :: tic => timer_tick
    procedure, public, pass :: toc => timer_tock
  end type basic_timer_t

contains


function get_time()
  integer(kind=8) :: ticks, ticks_per_sec, max_ticks
  real(dp) :: get_time

  call system_clock( count      = ticks,                                     &
                     count_rate = ticks_per_sec,                             &
                     count_max  = max_ticks )

  if ( ticks_per_sec == 0 ) then
    get_time = zero
  else
    get_time = real(ticks,dp) / real(ticks_per_sec,dp)
  end if
end function get_time

  subroutine timer_tick( this )
    class(basic_timer_t), intent(inout) :: this
    this%time_elapsed = zero
    this%time_start   = get_time()
  end subroutine timer_tick

  function timer_tock( this )
    class(basic_timer_t), intent(in) :: this
    real(dp)                         :: timer_tock
    timer_tock = get_time() - this%time_start
  end function timer_tock

end module timer_derived_type

module index_conversion
  implicit none
  private
  public :: global2local, global2local_bnd, global2local_ghost
  public :: local2global, local2global_bnd, local2global_ghost
  public :: in_bound, cell_face_nbors
  public :: get_face_idx_from_id

  interface cell_face_nbors
    module procedure cell_face_nbors_lin
    module procedure cell_face_nbors_sub
  end interface cell_face_nbors
contains
  
  pure function global2local(iG,nSub) result(iSub)
    integer,               intent(in) :: iG
    integer, dimension(:), intent(in) :: nSub
    integer, dimension(size(nSub)) :: iSub
    integer :: i, nDims, p, iGtmp, iTmp
    nDims = size(nSub)
    if (nDims==1) then
      iSub(1) = iG
      return
    end if
    p = product(nSub)
    iGtmp = iG
    do i = nDims,1,-1
      p = p/nSub(i)
      iTmp = mod(iGtmp-1,p) + 1
      iSub(i) = (iGtmp-iTmp)/p + 1
      iGtmp = iTmp
    end do
  end function global2local

  pure function local2global(iSub,nSub) result(iG)
    integer, dimension(:), intent(in) :: iSub, nSub
    integer :: iG
    integer :: nDims, p, i
    nDims = size(iSub)
    p = 1
    iG = 1
    do i = 1,nDims
        iG = iG + ( iSub(i) - 1 )*p
        p = p*nSub(i)
    end do
  end function local2global

  pure function global2local_ghost(iG,nSub,nGhost) result(iSub)
    integer,               intent(in) :: iG
    integer, dimension(:), intent(in) :: nSub, nGhost
    integer, dimension(size(nSub)) :: iSub, nSub2
    nSub2 = nSub + 2*nGhost
    iSub = global2local(iG,nSub2)
    iSub = iSub - nGhost
  end function global2local_ghost

  pure function local2global_ghost(iSub,nSub,nGhost) result(iG)
    integer, dimension(:), intent(in) :: iSub, nSub, nGhost
    integer, dimension(size(nSub)) :: iSub2, nSub2
    integer :: iG
    iSub2 = iSub + nGhost
    nSub2 = nSub + 2*nGhost
    iG = local2global(iSub2,nSub2)
  end function local2global_ghost

  pure function global2local_bnd(iG,lo,hi) result(iSub)
    integer,               intent(in) :: iG
    integer, dimension(:), intent(in) :: lo, hi
    integer, dimension(size(lo)) :: iSub, nSub
    nSub = hi - lo + 1
    iSub = global2local(iG,nSub)
    iSub = iSub + lo - 1
  end function global2local_bnd

  pure function local2global_bnd(iSub,lo,hi) result(iG)
    integer, dimension(:), intent(in) :: iSub, lo, hi
    integer, dimension(size(iSub)) :: idx, nSub
    integer :: iG
    idx  = iSub - lo + 1
    nSub = hi - lo + 1
    iG   = local2global(idx,nSub)
  end function local2global_bnd

  pure function in_bound( dim, idx, bnd_min, bnd_max )
    integer,                 intent(in) :: dim
    integer, dimension(dim), intent(in) :: idx, bnd_min, bnd_max
    logical                             :: in_bound
    in_bound =     all(idx>=bnd_min).and.all(idx<=bnd_max)                       &
              .or. all(idx<=bnd_min).and.all(idx>=bnd_max)
  end function in_bound

  pure subroutine cell_face_nbors_sub( dim, idx, bnd_min, bnd_max, nbor_cell_idx, nbor_face_id, n_int )
    integer,                       intent(in) :: dim
    integer, dimension(dim),       intent(in) :: idx, bnd_min, bnd_max
    integer, dimension(dim,2*dim), intent(out) :: nbor_cell_idx
    integer, dimension(2*dim),     intent(out) :: nbor_face_id
    integer,                       intent(out) :: n_int
    integer, dimension(dim,2*dim) :: nbor_cell_idx_tmp
    integer, dimension(2*dim) :: nbor_face_id_tmp
    integer, dimension(dim) :: idx_tmp
    integer :: s, j, n_ext, cnt
    cnt   = 0
    n_int = 0
    n_ext = 0
    do j = 1,dim
      do s = -1,1,2
        cnt = cnt + 1
        idx_tmp = idx
        idx_tmp(j) = idx_tmp(j) + s
        if ( in_bound(dim,idx_tmp,bnd_min,bnd_max) ) then
            n_int = n_int + 1
            nbor_cell_idx(:,n_int) = idx_tmp
            nbor_face_id(n_int) = cnt
        else
          n_ext = n_ext + 1
          nbor_cell_idx_tmp(:,n_ext) = idx_tmp
          nbor_face_id_tmp(n_ext) = cnt
        end if
      end do
    end do
    do j = 1,n_ext
      nbor_cell_idx(:,n_int+j) = nbor_cell_idx_tmp(:,j)
      nbor_face_id(n_int+j) = nbor_face_id_tmp(j)
    end do
  end subroutine cell_face_nbors_sub

  pure subroutine cell_face_nbors_lin( dim, lin_idx, bnd_min, bnd_max, &
                                       nbor_cell_idx, nbor_face_id, n_int )
    integer,                       intent(in) :: dim, lin_idx
    integer, dimension(dim),       intent(in) :: bnd_min, bnd_max
    integer, dimension(2*dim), intent(out) :: nbor_cell_idx
    integer, dimension(2*dim), intent(out) :: nbor_face_id
    integer,                       intent(out) :: n_int
    integer, dimension(dim,2*dim) :: nbor_idx
    integer, dimension(dim) :: idx
    integer :: s, j, n_ext, cnt
    idx = global2local_bnd(lin_idx,bnd_min,bnd_max)
    call cell_face_nbors_sub( dim, idx, bnd_min, bnd_max, nbor_idx, nbor_face_id, n_int )
    do j = 1,2*dim
      nbor_cell_idx(j) = local2global_bnd(nbor_idx(:,j),bnd_min,bnd_max)
    end do
  end subroutine cell_face_nbors_lin

  

  pure elemental subroutine get_face_info_from_id(face_id,dir,offset)
    integer, intent(in)  :: face_id
    integer, intent(out) :: dir, offset
    dir    = (face_id-1)/2 + 1
    offset = mod(face_id+1,2)
  end subroutine get_face_info_from_id

  pure subroutine get_face_idx_from_id(idx,face_id,dir,face_idx)
    integer, dimension(:),         intent(in) :: idx
    integer,                       intent(in)  :: face_id
    integer,                       intent(out) :: dir
    integer, dimension(size(idx)), intent(out) :: face_idx

    integer, dimension(size(idx)) :: face_offset
    integer :: offset
    call get_face_info_from_id(face_id,dir,offset)
    face_offset = 0
    face_offset(dir) = offset
    face_idx = idx + face_offset
  end subroutine get_face_idx_from_id

end module index_conversion

module combinatorics
  implicit none
  private
  public :: nchoosek
  public :: get_exponents
contains

  pure function nchoosek( n, k ) result( c )
    integer, intent(in) :: n, k
    integer             :: c
    integer :: i
    c = 0
    if (k>n) return

    c = 1
    do i = 1, min(n-k,k)
      c = c * ( n - (i-1) )
      c = c / i
    end do
  end function nchoosek

  pure subroutine get_exponents(n_dim,degree,n_terms,exponents,idx)
    use index_conversion, only : global2local
    integer, intent(in) :: n_dim, degree, n_terms
    integer, dimension(n_dim,n_terms), intent(out) :: exponents
    integer, dimension(degree+1),      intent(out) :: idx
    integer :: curr_total_degree, j, cnt, N_full_terms
    integer, dimension(n_dim) :: tmp_exp, nsub
    cnt = 0
    do curr_total_degree = 0,degree
      idx(curr_total_degree+1) = cnt + 1
      N_full_terms = (curr_total_degree+1) ** n_dim
      do j = 0,N_full_terms
        nSub = curr_total_degree + 1
        tmp_exp = global2local(j+1,nsub)-1
        if ( sum(tmp_exp) == curr_total_degree ) then
          cnt = cnt + 1
          exponents(:,cnt) = tmp_exp
        end if
      end do
    end do
  end subroutine get_exponents

end module combinatorics

module math
  use set_precision, only : dp
  implicit none
  private
  public :: cross_product, vector_norm
  public :: LUdecomp, LUsolve, mat_inv
  public :: maximal_diameter, maximal_extents
  public :: LegendrePolynomialAndDerivative, LegendreGaussNodesAndWeights
contains

    elemental subroutine LegendrePolynomialAndDerivative(N,x,LN,dLN)
    use set_constants, only : zero, one, two
    integer, intent(in) :: N
    real(dp), intent(in) :: x
    real(dp), intent(out) :: LN, dLN
    real(dp) :: LNm2, LNm1, dLNm2, dLNm1
    integer :: j
    if (N == 0) then
      LN = one
      dLN = zero
    elseif (N == 1) then
      LN = x
      dLN = one
    else
      LNm2 = one
      LNm1 = x
      dLNm2 = zero
      dLNm1 = one
      do j = 2,N
        LN = real(2*j-1,dp)/real(j,dp) * x * LNm1 &
          - real(j-1,dp)/real(j,dp) * LNm2
        dLN = dLNm2 + real(2*j-1,dp) * LNm1
        LNm2 = LNm1
        LNm1 = LN
        dLNm2 = dLNm1
        dLNm1 = dLN
      end do
    end if
  end subroutine LegendrePolynomialAndDerivative

  pure subroutine LegendreGaussNodesAndWeights(N,x,w)
    use set_constants, only : zero, one, two, four, third, pi
    integer,                  intent(in)  :: N
    real(dp), dimension(N+1), intent(out) :: x, w
    integer :: j, k
    real(dp) :: eps4, delta, LNp1, dLNp1
    integer, parameter :: quad_n_iter = 10
    eps4 = four*epsilon(one)
    x = zero
    w = zero

    if (N == 0) then
      x(1) = zero
      w(1) = two
    elseif (N == 1) then
      x(1) = -sqrt(third)
      w(1) = one
      x(2) = -x(1)
      w(2) = w(1)
    else
      do j = 0,(N+1)/2 - 1
        x(j+1) = -cos( ( real(2*j+1,dp)/real(2*N+2,dp) )*pi )
        do k = 1,quad_n_iter
          call LegendrePolynomialAndDerivative(N+1,x(j+1),LNp1,dLNp1)
          delta = -LNp1/dLNp1
          x(j+1) = x(j+1) + delta
          if ( abs(delta) <= eps4*abs(x(j+1)) ) then
            exit
          end if
        end do
        call LegendrePolynomialAndDerivative(N+1,x(j+1),LNp1,dLNp1)
        x(N+1-j) = -x(j+1)
        w(j+1) = two/( (one-x(j+1)**2)*dLNp1**2)
        w(N+1-j) = w(j+1)
      end do
      if (mod(N,2) == 0) then
        call LegendrePolynomialAndDerivative(N+1,zero,LNp1,dLNp1)
        x(N/2+1) = zero
        w(N/2+1) = two/dLNp1**2
      end if
    end if
  end subroutine LegendreGaussNodesAndWeights

  pure function cross_product( vec1, vec2 )
    real(dp), dimension(3), intent(in) :: vec1, vec2
    real(dp), dimension(3)             :: cross_product
    cross_product(1) =  ( vec1(2)*vec2(3) - vec1(3)*vec2(2) )
    cross_product(2) = -( vec1(1)*vec2(3) - vec1(3)*vec2(1) )
    cross_product(3) =  ( vec1(1)*vec2(2) - vec1(2)*vec2(1) )
  end function cross_product

  pure function vector_norm( vector )
    use set_precision, only : dp
    use set_constants, only : zero
    real(dp), dimension(:), intent(in) :: vector
    real(dp)                           :: vector_norm
    integer :: i
    vector_norm = zero
    do i = 1, size(vector)
      vector_norm = vector_norm + vector(i)**2
    end do
    vector_norm = sqrt( vector_norm )
  end function vector_norm

!================================== LUdecomp =================================80
!>
!! LU decomposition, with pivoting, for a regular NxN dense array
!<
!=============================================================================80
  subroutine LUdecomp( LU, P, A, m )
    use set_precision, only : dp
    use set_constants, only : zero, one
    real(dp), dimension(m,m), intent(out) :: LU,P
    real(dp), dimension(m,m), intent(in)  :: A
    integer,                  intent(in)  :: m
    real(dp), dimension(m) :: ctemp1, LUtemp
    integer  :: col, row, maxi, ipr
    real(dp) :: factor
    LU = A
    P = zero
    do col = 1,m
      P(col,col) = one
    end do
    do col = 1,m-1
!row pivot
      maxi=maxloc(abs(LU(col:m,col)),1)
      ipr=maxi+col-1
      if (ipr.ne.col) then
        ctemp1 = LU(ipr,:)
        LU(ipr,:) = LU(col,:)
        LU(col,:) = ctemp1
        ctemp1 = P(ipr,:)
        P(ipr,:) = P(col,:)
        P(col,:) = ctemp1
      end if
      if ( abs(LU(col,col)) > zero ) then
        do row = col+1,m
          factor = LU(row,col)/LU(col,col)
          LUtemp(col+1:m) = LU(row,col+1:m) - factor*LU(col,col+1:m)
          LU(row,col+1:m) = LUtemp(col+1:m)
          LU(row,col) = factor
        end do
      end if
    end do
  end subroutine LUdecomp

!================================== LUsolve ==================================80
!>
!! LU solve for a regular dense array
!<
!=============================================================================80
  subroutine LUsolve( x, LU, P, bin, m )
    use set_precision, only : dp
    real(dp), dimension(m),   intent(out) :: x
    real(dp), dimension(m,m), intent(in)  :: LU, P
    real(dp), dimension(m),   intent(in)  :: bin
    integer,                  intent(in)  :: m
    integer :: i, row
    real(dp), dimension(m) :: b, d
! Permute b matrix
    b = matmul(P,bin)
! Forward substitution
    d(1) = b(1)
    do row = 2,m
      d(row) = b(row) - sum( LU(row,1:row-1)*d(1:row-1) )
    end do
! Backward substitution
    x(m) = d(m)/LU(m,m)
    do i = 1,m-1
      row = m-i
      x(row) = ( d(row) - sum( LU(row,row+1:m)*x(row+1:m) ) ) / LU(row,row)
    end do
  end subroutine LUsolve

!================================== mat_inv ==================================80
!>
!! Calculates the inverse of a general nxn matrix using LU decomposition
!<
!=============================================================================80
  subroutine mat_inv( mat, inv, n )
    use set_precision, only : dp
    use set_constants, only : zero, one
    integer,                  intent(in)  :: n
    real(dp), dimension(n,n), intent(in)  :: mat
    real(dp), dimension(n,n), intent(out) :: inv
    integer                  :: i
    real(dp), dimension(n)   :: b
    real(dp), dimension(n,n) :: lu, p
    call ludecomp( lu, p, mat, n )
    inv = zero
    do i = 1,n
      b = zero
      b(i) = one
      call lusolve( inv(:,i), lu, p, b, n )
    end do
  end subroutine mat_inv

  pure function maximal_diameter(dim,n_points,points) result(d)
    use set_constants, only : zero
    integer, intent(in) :: dim, n_points
    real(dp), dimension(dim,n_points), intent(in) :: points
    real(dp) :: d
    real(dp), dimension(dim) :: delta
    integer :: i, j
    d = zero
    do j = 1,n_points-1
      do i = j+1,n_points
        delta = points(:,j) - points(:,i)
        d = max(d,dot_product(delta,delta))
      end do
    end do
    d = sqrt(d)
  end function maximal_diameter

  pure function maximal_extents(dim,n_points,points) result(d)
    use set_constants, only : half
    integer, intent(in) :: dim, n_points
    real(dp), dimension(dim,n_points), intent(in) :: points
    real(dp), dimension(dim) :: d
    d = half*(maxval(points,dim=2) - minval(points,dim=2))
  end function maximal_extents

end module math





module vector_derived_type
  use set_precision, only : dp
  use set_constants, only : zero
  implicit none
  private
  public :: face_vec
  public :: face_vec_ptr_3D
  type face_vec
    integer :: n
    real(dp), allocatable, dimension(:,:) :: v
  contains
    private
    procedure, public, pass :: create  => allocate_face_vec
    procedure, public, pass :: destroy => deallocate_face_vec
  end type face_vec

  type face_vec_ptr_3D
    type(face_vec), dimension(:,:,:), pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_face_vec_ptr_3D
  end type face_vec_ptr_3D

contains

  subroutine allocate_face_vec( this, n )
    class(face_vec), intent(inout) :: this
    integer,       intent(in)      :: n
    continue
    this%n = n
    allocate( this%v(3,n) )
    this%v = zero
  end subroutine allocate_face_vec

  pure elemental subroutine deallocate_face_vec( this )
    class(face_vec), intent(inout) :: this
    continue
    this%n = 0
    if( allocated( this%v  ) ) deallocate( this%v )
  end subroutine deallocate_face_vec

  pure elemental subroutine destroy_face_vec_ptr_3D( this )
    class(face_vec_ptr_3D), intent(inout) :: this
    this%p => null()
  end subroutine destroy_face_vec_ptr_3D
end module vector_derived_type

module pointers
  use set_precision, only : dp
  implicit none
  private
  public :: array_ptr_3D_real, array_ptr_4D_real

  type array_ptr_3D_real
    real(dp), dimension(:,:,:),     pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_real_3D
  end type array_ptr_3D_real

  type array_ptr_4D_real
    real(dp), dimension(:,:,:,:),   pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_real_4D
  end type array_ptr_4D_real

contains

  pure elemental subroutine destroy_real_3D( this )
    class(array_ptr_3D_real), intent(inout) :: this
    this%p => null()
  end subroutine destroy_real_3D

  pure elemental subroutine destroy_real_4D( this )
    class(array_ptr_4D_real), intent(inout) :: this
    this%p => null()
  end subroutine destroy_real_4D
end module pointers

module linspace_helper
  use set_precision, only : dp
  private
  public :: unit_cartesian_mesh_cat
  public :: linspace, meshgrid2, meshgrid3
contains
  pure function unit_cartesian_mesh_cat(nx,ny,nz) result(xyz)
    integer, intent(in) :: nx, ny, nz
    real(dp), dimension(3,nx,ny,nz) :: xyz
    real(dp), dimension(nx,ny,nz) :: tmp_x, tmp_y, tmp_z
    integer :: i, j, k

    call unit_cartesian_mesh(nx,ny,nz,tmp_x,tmp_y,tmp_z)

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          xyz(1,i,j,k) = tmp_x(i,j,k)
          xyz(2,i,j,k) = tmp_y(i,j,k)
          xyz(3,i,j,k) = tmp_z(i,j,k)
        end do
      end do
    end do
  end function unit_cartesian_mesh_cat

  pure subroutine unit_cartesian_mesh(nx,ny,nz,x,y,z)
    use set_constants, only : zero, one
    integer, intent(in) :: nx, ny, nz
    real(dp), dimension(nx,ny,nz), intent(out) :: x, y, z

    call meshgrid3( linspace(nx,zero,one), &
                    linspace(ny,zero,one), &
                    linspace(nz,zero,one), x,y,z)
  end subroutine unit_cartesian_mesh

  pure function linspace(N,x1,x2) result(array)
    integer,  intent(in)   :: N
    real(dp), intent(in)   :: x1, x2
    real(dp), dimension(N) :: array
    real(dp) :: range_den
    integer :: i
    if (N==0) return
    if (N==1) then
      array(1) = x1
      return
    end if
    range_den = (x2-x1)/real(N-1,dp)
    do i = 1,N
      array(i) = x1 + range_den*real(i-1,dp)
    end do
  end function linspace

  pure subroutine meshgrid2(x1,x2,x1_array,x2_array)
    real(dp), dimension(:),   intent(in)  :: x1, x2
    real(dp), dimension(:,:), intent(out) :: x1_array, x2_array
    integer :: N1, N2
    N1 = size(x1)
    N2 = size(x2)
    x1_array = spread(x1,2,N2)
    x2_array = spread(x2,1,N1)
  end subroutine meshgrid2

  pure subroutine meshgrid3(x1,x2,x3,x1_array,x2_array,x3_array)
    real(dp), dimension(:),     intent(in)  :: x1, x2, x3
    real(dp), dimension(:,:,:), intent(out) :: x1_array, x2_array, x3_array
    real(dp), dimension(size(x1),size(x2)) :: x1_tmp
    real(dp), dimension(size(x2),size(x3)) :: x2_tmp
    real(dp), dimension(size(x2),size(x3),size(x1)) :: x2_tmp2
    real(dp), dimension(size(x3),size(x1)) :: x3_tmp
    real(dp), dimension(size(x3),size(x1),size(x2)) :: x3_tmp2
    integer, parameter, dimension(3) :: o2 = [2,3,1], o3 = [3,1,2]
    integer :: N1, N2, N3
    N1 = size(x1)
    N2 = size(x2)
    N3 = size(x3)

    x1_tmp   = spread(x1,2,N2)
    x2_tmp   = spread(x2,2,N3)
    x3_tmp   = spread(x3,2,N1)
    x1_array = spread(x1_tmp,3,N3)
    x2_tmp2  = spread(x2_tmp,3,N1)
    x3_tmp2  = spread(x3_tmp,3,N2)
    x2_array = reshape(x2_tmp2,shape(x2_array),order=o2)
    x3_array = reshape(x3_tmp2,shape(x3_array),order=o3)
  end subroutine meshgrid3

end module linspace_helper

module lagrange_interpolation
  use set_precision, only : dp
  use set_constants, only : zero, one, two, half
  implicit none
  private
  ! public :: xb, wb, Dmat
  public :: generate_1D_barycentric_info, destroy_1D_barycentric_info
  public :: lagbary, lagbary_wderiv
  public :: lagbary_2D, lagbary_2D_wgrad
  public :: lagbary_3D, lagbary_3D_wgrad
  public :: calc_grid_metrics, jacobian_determinant
  public :: normal_vectors
  interface jacobian_determinant
    module procedure jacobian_determinant_2D
    module procedure jacobian_determinant_3D
  end interface jacobian_determinant

  integer :: Nmax = 10
  real(dp), dimension(:,:,:,:), allocatable :: Dmat
  real(dp), dimension(:,:),   allocatable :: xb, wb
contains

  elemental logical function almost_equal(a,b)
    real(dp), intent(in) :: a, b
    logical :: test1, test2, test3
    test1 = ( (a==zero) .or. (b==zero) )
    test2 = ( abs(a-b) <= two*epsilon(one) )
    test3 = ( ( abs(a-b) <= epsilon(abs(a)) ) .and. &
            ( abs(a-b) <= epsilon(abs(b)) ) )
    almost_equal = ( ( test1 .and. test2 ) .or. ( (.not. test1) .and. test3 ) )
  end function almost_equal

  function barycentric_weights(x) result(w)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x))       :: w
    integer :: j, k, N
    N = size(x)
    w = one
    do j = 2,N
      do k = 1,j-1
        w(k) = w(k) * ( x(k) - x(j) )
        w(j) = w(j) * ( x(j) - x(k) )
      end do
    end do
    w = one/w
  end function barycentric_weights

  pure function polynomial_derivative_matrix(x,w) result(D)
    real(dp), dimension(:), intent(in) :: x, w
    real(dp), dimension(size(x),size(x)) :: D
    integer :: i, j, N
    D = zero
    N = size(x)
    do i = 1,N
      do j = 1,N
        if (j/=i) then
          D(i,j) = w(j)/w(i) * one / ( x(i) - x(j) )
          D(i,i) = D(i,i) - D(i,j)
        end if
      end do
    end do
  end function polynomial_derivative_matrix

  pure function mth_order_polynomial_derivative_matrix(x,w,M) result(D)
    real(dp), dimension(:), intent(in) :: x, w
    integer,                intent(in) :: M
    real(dp), dimension(size(x),size(x),M) :: D
    integer :: i, j, k, N
    D = zero
    N = size(x)
    D(:,:,1) = polynomial_derivative_matrix(x,w)
    do k = 2,M
      do i = 1,N
        D(i,i,k) = zero
        do j = 1,N
          if (j/=i) then
            D(i,j,k) = ( real(k,dp) / (x(i) - x(j)) )                          &
                     * ( w(j)/w(i)*D(i,i,k-1) - D(i,j,k-1) )
            D(i,i,k) = D(i,i,k) - D(i,j,k)
          end if
        end do
      end do
    end do
  end function mth_order_polynomial_derivative_matrix

  subroutine generate_1D_barycentric_info(N)
    use linspace_helper, only : linspace
    integer, intent(in), optional :: N
    integer :: j
    real(dp) :: x1 = -one, x2 = one
    if ( present(N) ) then
      if ( N > Nmax ) then
        Nmax = N
        continue
      end if
    end if
    if ( allocated(Dmat) ) deallocate(Dmat)
    if ( allocated(xb) )   deallocate(xb)
    if ( allocated(wb) )   deallocate(wb)
    allocate( Dmat(Nmax,Nmax,Nmax,2), xb(Nmax,Nmax), wb(Nmax,Nmax) )
    Dmat = zero; xb = zero; wb = zero
    wb(1,1) = one
    do j = 2,Nmax
      ! call linspace(x1,x2,xb(1:j,j))
      xb(1:j,j) = linspace(j,x1,x2)
      wb(1:j,j) = barycentric_weights( xb(1:j,j) )
      Dmat(1:j,1:j,j,:) = mth_order_polynomial_derivative_matrix( xb(1:j,j), wb(1:j,j), 2 )
    end do
  end subroutine generate_1D_barycentric_info

  subroutine destroy_1D_barycentric_info
    if ( allocated(Dmat) ) deallocate(Dmat)
    if ( allocated(xb) ) deallocate(xb)
    if ( allocated(wb) ) deallocate(wb)
  end subroutine destroy_1D_barycentric_info

  pure subroutine lagbary(x,dir,fval,Npts,val)
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val
    real(dp) :: A, F
    real(dp) :: x1, t1
    integer :: j, N
    A = zero
    F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val = fval(j)
        return
      end if
      t1 = wb(j,N)/x1
      A = A + t1 * fval(j)
      F = F + t1
    end do
    val = A/F
  end subroutine lagbary

  pure subroutine lagbary_wderiv(x,dir,fval,Npts,val,dval)
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val, dval
    real(dp) :: A, B, C, F
    real(dp) :: x1, t1, t2, FF, AC
    integer :: j, N
    A = zero
    B = zero
    C = zero
    F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val = fval(j)
        dval = dot_product( Dmat(j,1:N,N,1), fval )
        return
      end if
      t1 = wb(j,N)/x1
      A = A + t1 * fval(j)
      F = F + t1
      t2 = t1/x1
      B = B + t2 * fval(j)
      C = C + t2
    end do
    val = A/F
    FF = F*F
    AC = A*C
    dval = (B * F - AC)/FF
  end subroutine lagbary_wderiv

  pure subroutine lagbary_wderiv2(x,dir,fval,Npts,val,dval,d2val)
    real(dp),               intent(in)  :: x
    integer,                intent(in)  :: dir
    integer,  dimension(:), intent(in)  :: Npts
    real(dp), dimension(:), intent(in)  :: fval
    real(dp),               intent(out) :: val, dval, d2val
    real(dp) :: A, B, C, D, E, F
    real(dp) :: x1, t1, t2, t3, FF, AC
    integer :: j, N
    A = zero
    B = zero
    C = zero
    D = zero
    E = zero
    F = zero
    N = Npts(dir)
    do j = 1,N
      x1 = xb(j,N) - x
      if ( almost_equal(x1,zero) ) then
        val   = fval(j)
        dval  = dot_product( Dmat(j,1:N,N,1), fval )
        d2val = dot_product( Dmat(j,1:N,N,2), fval )
        return
      end if
      t1 = wb(j,N)/x1
      A = A + t1 * fval(j)
      F = F + t1
      t2 = t1/x1
      B = B + t2 * fval(j)
      C = C + t2
      t3 = t2/x1
      D = D + t3 * fval(j)
      E = E + t3
    end do
    val = A/F
    FF = F*F
    AC = A*C
    dval = (B * F - AC)/FF
    d2val = ( two * D ) / F - ( two * E * A ) / FF - ( two * B * C ) / FF      &
          + ( two * C * AC ) / ( FF * F )
  end subroutine lagbary_wderiv2

  pure subroutine lagbary_2D(x,fval,Npts,val)
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(size(fval,2)) :: tmp
    integer :: j
    do j = 1,Npts(2)
      call lagbary( x(1), 1, fval(:,j), Npts, tmp(j) )
    end do
    call lagbary( x(2), 2, tmp, Npts, val )
  end subroutine lagbary_2D

  pure subroutine lagbary_2D_wgrad(x,fval,Npts,val,grad)
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(2),   intent(out) :: grad
    real(dp), dimension(size(fval,2)) :: tmp, gtmp
    integer :: j
    do j = 1,Npts(2)
      call lagbary_wderiv( x(1), 1, fval(:,j), Npts, tmp(j), gtmp(j) )
    end do
    call lagbary_wderiv( x(2), 2,  tmp, Npts, val, grad(2) )
    call lagbary(        x(2), 2, gtmp, Npts,      grad(1) )
  end subroutine lagbary_2D_wgrad

  pure subroutine lagbary_2D_whess(x,fval,Npts,val,grad,hess)
    real(dp), dimension(2),   intent(in)  :: x
    real(dp), dimension(:,:), intent(in)  :: fval
    integer,  dimension(2),   intent(in)  :: Npts
    real(dp),                 intent(out) :: val
    real(dp), dimension(2),   intent(out) :: grad
    real(dp), dimension(3),   intent(out) :: hess
    real(dp), dimension(size(fval,2)) :: tmp, gtmp, htmp
    integer :: j
    do j = 1,Npts(2)
      call lagbary_wderiv2( x(1), 1, fval(:,j), Npts, tmp(j), gtmp(j), htmp(j) )
    end do
    call lagbary_wderiv2( x(2), 2,  tmp, Npts, val, grad(2), hess(3) )
    call lagbary_wderiv(  x(2), 2, gtmp, Npts,      grad(1), hess(2) )
    call lagbary(         x(2), 2, htmp, Npts,               hess(1) )
  end subroutine lagbary_2D_whess

  pure subroutine lagbary_3D(x,fval,Npts,val)
    real(dp), dimension(3),     intent(in)  :: x
    real(dp), dimension(:,:,:), intent(in)  :: fval
    integer,  dimension(3),     intent(in)  :: Npts
    real(dp),                   intent(out) :: val
    real(dp), dimension(size(fval,2),size(fval,3)) :: tmp
    real(dp), dimension(size(fval,3)) :: tmp2
    integer :: k, j
    do k = 1,Npts(3)
      do j = 1,Npts(2)
        call lagbary( x(1), 1, fval(:,j,k), Npts, tmp(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call lagbary( x(2), 2, tmp(:,k), Npts, tmp2(k) )
    end do
    call lagbary( x(3), 3, tmp2, Npts, val )
  end subroutine lagbary_3D

  pure subroutine lagbary_3D_wgrad(x,fval,Npts,val,grad)
    real(dp), dimension(3),     intent(in)  :: x
    real(dp), dimension(:,:,:), intent(in)  :: fval
    integer,  dimension(3),     intent(in)  :: Npts
    real(dp),                   intent(out) :: val
    real(dp), dimension(3),     intent(out) :: grad
    real(dp), dimension(size(fval,2),size(fval,3)) :: tmp, gtmp0
    real(dp), dimension(size(fval,3)) :: tmp2, gtmp1, gtmp2
    integer :: k, j
    do k = 1,Npts(3)
      do j = 1,Npts(2)
        call lagbary_wderiv( x(1), 1, fval(:,j,k), Npts, tmp(j,k), gtmp0(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call lagbary_wderiv( x(2), 2,   tmp(:,k), Npts, tmp2(k), gtmp2(k) )
      call lagbary(        x(2), 2, gtmp0(:,k), Npts, gtmp1(k) )
    end do
    call lagbary_wderiv( x(3), 3,  tmp2, Npts, val, grad(3) )
    call lagbary(        x(3), 3, gtmp2, Npts,      grad(2) )
    call lagbary(        x(3), 3, gtmp1, Npts,      grad(1) )
  end subroutine lagbary_3D_wgrad

  pure subroutine lagbary_3D_whess(x,fval,Npts,val,grad,hess)
    real(dp), dimension(3),     intent(in)  :: x
    real(dp), dimension(:,:,:), intent(in)  :: fval
    integer,  dimension(3),     intent(in)  :: Npts
    real(dp),                   intent(out) :: val
    real(dp), dimension(3),     intent(out) :: grad
    real(dp), dimension(6),     intent(out) :: hess
    real(dp), dimension(size(fval,2),size(fval,3)) :: tmp, gtmp, htmp
    real(dp), dimension(size(fval,3)) :: tmp1, gtmp1, gtmp2, htmp1, htmp2, htmp3
    integer :: k, j
    do k = 1,Npts(3)
      do j = 1,Npts(2)
        call lagbary_wderiv2( x(1), 1, fval(:,j,k), Npts, tmp(j,k), gtmp(j,k), htmp(j,k) )
      end do
    end do
    do k = 1,Npts(3)
      call lagbary_wderiv2( x(2), 2,  tmp(:,k), Npts, tmp1(k), gtmp2(k), htmp3(k) )
      call lagbary_wderiv(  x(2), 2, gtmp(:,k), Npts,          gtmp1(k), htmp2(k) )
      call lagbary(         x(2), 2, htmp(:,k), Npts,                    htmp1(k) )
    end do
    call lagbary_wderiv2( x(3), 3,  tmp1, Npts, val, grad(3), hess(6) )
    call lagbary_wderiv(  x(3), 3, gtmp2, Npts,      grad(2), hess(5) )
    call lagbary(         x(3), 3, htmp3, Npts,               hess(4) )
    call lagbary_wderiv(  x(3), 3, gtmp1, Npts,      grad(1), hess(3) )
    call lagbary(         x(3), 3, htmp2, Npts,               hess(2) )
    call lagbary(         x(3), 3, htmp1, Npts,               hess(1) )
  end subroutine lagbary_3D_whess

  pure function covariant_base_vectors_3D(point,X1,X2,X3) result(a)
    real(dp), dimension(3),     intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: a
    real(dp) :: junk
    integer, dimension(3) :: Npts
    Npts = shape(X1)
    call lagbary_3D_wgrad(point,X1,Npts,junk,a(1,:))
    call lagbary_3D_wgrad(point,X2,Npts,junk,a(2,:))
    call lagbary_3D_wgrad(point,X3,Npts,junk,a(3,:))
  end function covariant_base_vectors_3D

  function contravariant_base_vectors_3D(point,X1,X2,X3) result(a)
    real(dp), dimension(3),     intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: a
    real(dp) :: J
    a = calc_grid_metrics(point,X1,X2,X3)
    J = jacobian_determinant_3D(point,X1,X2,X3)
    a = a / J
  end function contravariant_base_vectors_3D

  function calc_grid_metrics(point,X1,X2,X3) result(Ja)
    use pointers, only : array_ptr_3D_real
    real(dp), dimension(3), intent(in) :: point
    real(dp), dimension(:,:,:), intent(in), target :: X1, X2, X3
    real(dp), dimension(3,3) :: Ja
    real(dp), dimension(size(X1,1),size(X1,2),size(X1,3)) :: X_l, X_m, tmp
    type(array_ptr_3D_real), dimension(3) :: X
    integer, dimension(3) :: Npts
    integer :: i
    real(dp) :: junk
    real(dp), dimension(3) :: dX_l, dX_m, dd1, dd2, dd3
    Ja = zero
    X(1)%p => X1
    X(2)%p => X2
    X(3)%p => X3
    Npts = shape(X1)
    do i = 1,3
      X_l = X(mod(4+i,3)+1)%p
      X_m = X(mod(3+i,3)+1)%p
      call lagbary_3D_wgrad( point, X_l, Npts, junk, dX_l )
      call lagbary_3D_wgrad( point, X_m, Npts, junk, dX_m )
      tmp = X_l*dX_m(1) - X_m*dX_l(1)
      call lagbary_3D_wgrad( point, tmp, Npts, junk, dd1 )
      tmp = X_l*dX_m(2) - X_m*dX_l(2)
      call lagbary_3D_wgrad( point, tmp, Npts, junk, dd2 )
      tmp = X_l*dX_m(3) - X_m*dX_l(3)
      call lagbary_3D_wgrad( point, tmp, Npts, junk, dd3 )
      Ja(i,1) = -half*( dd3(2) - dd2(3) );
      Ja(i,2) = -half*( dd1(3) - dd3(1) );
      Ja(i,3) = -half*( dd2(1) - dd1(2) );
    end do
    nullify( X(1)%p, X(2)%p, X(3)%p ) ! may not be necessary
  end function calc_grid_metrics


  pure function normal_vectors(point,X1,X2,X3) result(Nvec)
    use math, only : cross_product, vector_norm
    implicit none
    real(dp), dimension(3), intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp), dimension(3,3) :: Nvec
    integer, dimension(3) :: Npts
    real(dp), dimension(3,3) :: A
    real(dp) :: junk
    Npts = shape(X1)
    call lagbary_3D_wgrad(point,X1,Npts,junk,A(:,1))
    call lagbary_3D_wgrad(point,X2,Npts,junk,A(:,2))
    call lagbary_3D_wgrad(point,X3,Npts,junk,A(:,3))

    Nvec(:,1) = cross_product(A(:,2),A(:,3))
    Nvec(:,2) = cross_product(A(:,1),A(:,3))
    Nvec(:,3) = cross_product(A(:,1),A(:,2))

    Nvec(:,1) = Nvec(:,1)/vector_norm(Nvec(:,1))
    Nvec(:,2) = Nvec(:,2)/vector_norm(Nvec(:,2))
    Nvec(:,3) = Nvec(:,3)/vector_norm(Nvec(:,3))
  end function normal_vectors

  pure function jacobian_determinant_2D(point,X1,X2,X3) result(Jac)
    implicit none
    real(dp), dimension(2), intent(in) :: point
    real(dp), dimension(:,:), intent(in) :: X1, X2, X3
    real(dp) :: Jac, junk
    integer, dimension(2) :: Npts
    real(dp), dimension(2,3) :: A
    Npts = shape(X1)
    call lagbary_2D_wgrad(point,X1,Npts,junk,A(:,1))
    call lagbary_2D_wgrad(point,X2,Npts,junk,A(:,2))
    call lagbary_2D_wgrad(point,X3,Npts,junk,A(:,3))
    Jac = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  end function jacobian_determinant_2D

  pure function jacobian_determinant_3D(point,X1,X2,X3) result(Jac)
    implicit none
    real(dp), dimension(3), intent(in) :: point
    real(dp), dimension(:,:,:), intent(in) :: X1, X2, X3
    real(dp) :: Jac, junk
    integer, dimension(3) :: Npts
    real(dp), dimension(3,3) :: A
    Npts = shape(X1)
    call lagbary_3D_wgrad(point,X1,Npts,junk,A(:,1))
    call lagbary_3D_wgrad(point,X2,Npts,junk,A(:,2))
    call lagbary_3D_wgrad(point,X3,Npts,junk,A(:,3))
    Jac = A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(3,2)*A(1,3) + A(3,1)*A(1,2)*A(2,3) &
        - A(3,1)*A(2,2)*A(1,3) - A(2,1)*A(1,2)*A(3,3) - A(1,1)*A(3,2)*A(2,3)
  end function jacobian_determinant_3D

end module lagrange_interpolation

module quadrature_derived_type

  use set_precision,       only : dp
  use set_constants,       only : zero
  implicit none
  private
  public :: quad_t
  public :: quad_ptr, quad_ptr_3D
  public :: create_quad_ref_1D, create_quad_ref_2D, create_quad_ref_3D
  public :: map_quad_ref_to_physical_1D
  public :: map_quad_ref_to_physical_2D
  public :: map_quad_ref_to_physical_3D
  type quad_t
    integer :: n_quad = 0
    real(dp), allocatable, dimension(:,:) :: quad_pts
    real(dp), allocatable, dimension(:)   :: quad_wts
  contains
    private
    procedure, public, pass :: create  => allocate_quad
    procedure, public, pass :: destroy => deallocate_quad
    generic,   public :: integrate => integrate_scalar, integrate_vector
    procedure :: integrate_scalar
    procedure :: integrate_vector
  end type quad_t

  type quad_ptr
    type(quad_t), pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_quad_ptr
  end type quad_ptr

  type quad_ptr_3D
    type(quad_t), dimension(:,:,:), pointer :: p => null()
  contains
    private
    procedure, public, pass :: destroy => destroy_quad_ptr_3D
  end type quad_ptr_3D

contains

  pure elemental subroutine destroy_quad_ptr_3D( this )
    class(quad_ptr_3D), intent(inout) :: this
    this%p => null()
  end subroutine destroy_quad_ptr_3D

  pure elemental subroutine destroy_quad_ptr( this )
    class(quad_ptr), intent(inout) :: this
    this%p => null()
  end subroutine destroy_quad_ptr

  pure elemental subroutine allocate_quad( this, n_quad )
    use set_constants, only : zero
    class(quad_t), intent(inout) :: this
    integer,       intent(in)    :: n_quad
    this%n_quad = n_quad
    allocate( this%quad_pts(3,n_quad) )
    this%quad_pts = zero
    allocate( this%quad_wts(n_quad) )
    this%quad_wts = zero
  end subroutine allocate_quad

  pure elemental subroutine deallocate_quad( this )
    class(quad_t), intent(inout) :: this
    this%n_quad = 0
    if( allocated( this%quad_wts  ) ) deallocate( this%quad_wts  )
    if( allocated( this%quad_pts  ) ) deallocate( this%quad_pts  )
  end subroutine deallocate_quad

  pure function integrate_scalar( this, f ) result( integral )
    use set_precision, only : dp
    class(quad_t),                    intent(in) :: this
    real(dp), dimension(this%n_quad), intent(in) :: f
    real(dp)                                     :: integral
    integral = dot_product(f,this%quad_wts)
  end function integrate_scalar

  pure function integrate_vector( this, neq, f ) result( integral )
    use set_precision, only : dp
    class(quad_t),                        intent(in) :: this
    integer,                              intent(in) :: neq
    real(dp), dimension(neq,this%n_quad), intent(in) :: f
    real(dp), dimension(neq)                         :: integral
    integer :: n
    do n = 1, neq
      integral(n) = dot_product(f(n,:),this%quad_wts)
    end do
  end function integrate_vector

  pure function gauss_1D_size( polynomial_order ) result( N_quad )
    use set_constants, only : half
    integer, intent(in) :: polynomial_order
    integer             :: N_quad
    N_quad = ceiling( half*(polynomial_order + 1) )
  end function gauss_1D_size

  pure subroutine gauss_1D( n_quad, pts_1D, wts_1D )
    use math, only : LegendreGaussNodesAndWeights
    integer,                       intent(in)  :: n_quad
    real(dp), dimension( n_quad ), intent(out) :: pts_1D
    real(dp), dimension( n_quad ), intent(out) :: wts_1D
    call LegendreGaussNodesAndWeights(n_quad-1, pts_1D, wts_1D)
  end subroutine gauss_1D

  pure subroutine create_quad_ref_1D( quad_order, quad_ref )
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: xtmp
    integer :: n_quad
    n_quad = gauss_1D_size( quad_order )
    call quad_ref%destroy()
    call quad_ref%create( n_quad )
    call gauss_1D( n_quad, xtmp, quad_ref%quad_wts )
    quad_ref%quad_pts(1,:) = xtmp
  end subroutine create_quad_ref_1D

  pure subroutine create_quad_ref_2D( quad_order, quad_ref )
    use set_constants, only : zero
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    integer :: n_quad
    integer :: i, j, cnt
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: pts_1D
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: wts_1D
    n_quad = gauss_1D_size( quad_order )
    call gauss_1D(n_quad, pts_1D, wts_1D)
    call quad_ref%destroy()
    call quad_ref%create( n_quad**2 )
    cnt = 0
    do j = 1, n_quad
      do i = 1, n_quad
        cnt = cnt + 1
        quad_ref%quad_pts(:,cnt) = [ pts_1D(i), pts_1D(j), zero ]
        quad_ref%quad_wts(cnt) = wts_1D(i)*wts_1D(j)
      end do
    end do
  end subroutine create_quad_ref_2D

  pure subroutine create_quad_ref_3D( quad_order, quad_ref )
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    integer :: n_quad
    integer :: i, j, k, cnt
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: pts_1D
    real(dp), dimension( gauss_1D_size( quad_order ) ) :: wts_1D
    n_quad = gauss_1D_size( quad_order )
    call gauss_1D(n_quad, pts_1D, wts_1D)
    call quad_ref%destroy()
    call quad_ref%create( n_quad**3 )
    cnt = 0
    do k = 1, n_quad
      do j = 1, n_quad
        do i = 1, n_quad
          cnt = cnt + 1
          quad_ref%quad_pts(:,cnt) = [ pts_1D(i), pts_1D(j), pts_1D(k) ]
          quad_ref%quad_wts(cnt) = wts_1D(i)*wts_1D(j)*wts_1D(k)
        end do
      end do
    end do
  end subroutine create_quad_ref_3D

  pure subroutine map_quad_ref_to_physical_1D( x1nodes, x2nodes, x3nodes,      &
                                               quad_ref, quad_physical )
    use set_constants,           only : zero, half
    use math,                    only : vector_norm
    use lagrange_interpolation,  only : lagbary_wderiv

    real(dp), dimension(:), intent(in)    :: x1nodes, x2nodes, x3nodes
    type(quad_t),           intent(in)  :: quad_ref
    type(quad_t),           intent(out) :: quad_physical
    integer  :: n, i
    integer, dimension(1) :: Npts
    real(dp) :: det_jac
    real(dp), dimension(3) :: node_diff

    real(dp), dimension(3) :: tangent, normal, binormal
    real(dp) :: kappa

    continue
    Npts = shape(x1nodes)

    if( quad_ref%n_quad /= quad_physical%n_quad ) then
      call quad_physical%destroy()
      call quad_physical%create(quad_ref%n_quad)
    end if

    quad_physical%quad_pts = zero
    do n = 1,quad_ref%n_quad
      call lagbary_wderiv( quad_ref%quad_pts(1,n), 1, x1nodes, Npts, &
                           quad_physical%quad_pts(1,n), node_diff(1) )
      call lagbary_wderiv( quad_ref%quad_pts(1,n), 1, x2nodes, Npts, &
                           quad_physical%quad_pts(2,n), node_diff(2) )
      call lagbary_wderiv( quad_ref%quad_pts(1,n), 1, x3nodes, Npts, &
                           quad_physical%quad_pts(3,n), node_diff(3) )
      det_jac = vector_norm(node_diff)
      quad_physical%quad_wts(n) = det_jac * quad_ref%quad_wts(n)

    end do
  end subroutine map_quad_ref_to_physical_1D

  pure subroutine map_quad_ref_to_physical_2D( x1nodes, x2nodes, x3nodes,    &
                                                  quad_ref, quad_physical )
    use set_constants,           only : zero
    use lagrange_interpolation,  only : lagbary_2D, jacobian_determinant

    real(dp),     dimension(:,:), intent(in)  :: x1nodes, x2nodes, x3nodes
    type(quad_t),                 intent(in)  :: quad_ref
    type(quad_t),                 intent(out) :: quad_physical
    integer  :: n
    real(dp) :: det_jac
    integer,      dimension(2) :: Npts
    continue

    Npts = shape(x1nodes)
    if( quad_ref%n_quad /= quad_physical%n_quad ) then
      call quad_physical%destroy()
      call quad_physical%create(quad_ref%n_quad)
    end if

    quad_physical%quad_pts = zero
    do n = 1,quad_ref%n_quad
      call lagbary_2D( quad_ref%quad_pts(1:2,n), x1nodes, Npts, &
                  quad_physical%quad_pts(  1,n) )
      call lagbary_2D( quad_ref%quad_pts(1:2,n), x2nodes, Npts, &
                  quad_physical%quad_pts(  2,n) )
      call lagbary_2D( quad_ref%quad_pts(1:2,n), x3nodes, Npts, &
                  quad_physical%quad_pts(  3,n) )
      det_jac = jacobian_determinant( quad_ref%quad_pts(1:2,n), &
                  x1nodes, x2nodes, x3nodes )
      quad_physical%quad_wts(n) = abs(det_jac) * quad_ref%quad_wts(n)
    end do
  end subroutine map_quad_ref_to_physical_2D

  pure subroutine map_quad_ref_to_physical_3D( x1nodes, x2nodes, x3nodes,    &
                                                  quad_ref, quad_physical )
    use set_constants,           only : zero
    use lagrange_interpolation,  only : lagbary_3D, jacobian_determinant
    real(dp), dimension(:,:,:), intent(in)  :: x1nodes, x2nodes, x3nodes
    type(quad_t),               intent(in)  :: quad_ref
    type(quad_t),               intent(out) :: quad_physical
    integer  :: n
    real(dp) :: det_jac
    integer,  dimension(3) :: Npts
    continue

    Npts = shape(x1nodes)
    if( quad_ref%n_quad /= quad_physical%n_quad ) then
      call quad_physical%destroy()
      call quad_physical%create(quad_ref%n_quad)
    end if

    quad_physical%quad_pts = zero
    do n = 1,quad_ref%n_quad
      call lagbary_3D( quad_ref%quad_pts(:,n), x1nodes, Npts, &
                      quad_physical%quad_pts(  1,n) )
      call lagbary_3D( quad_ref%quad_pts(:,n), x2nodes, Npts, &
                      quad_physical%quad_pts(  2,n) )
      call lagbary_3D( quad_ref%quad_pts(:,n), x3nodes, Npts, &
                      quad_physical%quad_pts(  3,n) )
      det_jac = jacobian_determinant( quad_ref%quad_pts(:,n), &
                                      x1nodes, x2nodes, x3nodes )
      quad_physical%quad_wts(n) = det_jac * quad_ref%quad_wts(n)
    end do
  end subroutine map_quad_ref_to_physical_3D
end module quadrature_derived_type

module grid_derived_type
  use set_precision,           only : dp
  use quadrature_derived_type, only : quad_t, quad_ptr_3D
  use vector_derived_type,     only : face_vec, face_vec_ptr_3D
  use pointers,                only : array_ptr_3D_real, array_ptr_4D_real
  implicit none
  private
  public :: derived_grid_vars
  public :: grid_block
  public :: grid_type
  public :: deallocate_grid
  public :: allocate_grid_block, deallocate_grid_block
  public :: allocate_derived_grid, deallocate_derived_grid

  public :: pack_cell_node_coords
  public :: get_face_quad_ptrs

  type derived_grid_vars
    real(dp),       allocatable, dimension(:,:,:,:) :: cell_c
    real(dp),       allocatable, dimension(:,:,:)   :: volume
    real(dp),       allocatable, dimension(:,:,:)   :: xi_area
    real(dp),       allocatable, dimension(:,:,:)   :: eta_area
    real(dp),       allocatable, dimension(:,:,:)   :: zeta_area
    type(quad_t),   allocatable, dimension(:,:,:)   :: quad
    type(quad_t),   allocatable, dimension(:,:,:)   :: xi_face_quad
    type(quad_t),   allocatable, dimension(:,:,:)   :: eta_face_quad
    type(quad_t),   allocatable, dimension(:,:,:)   :: zeta_face_quad
    type(face_vec), allocatable, dimension(:,:,:)   :: xi_nv
    type(face_vec), allocatable, dimension(:,:,:)   :: eta_nv
    type(face_vec), allocatable, dimension(:,:,:)   :: zeta_nv
    type(quad_ptr_3D),           dimension(3)       :: face_quads
    type(face_vec_ptr_3D),       dimension(3)       :: normals
    integer, dimension(:), pointer :: n_cells, n_ghost
    integer,               pointer :: n_dim
  contains
    private
    procedure, public, pass :: setup   =>   allocate_derived_grid
    procedure, public, pass :: destroy => deallocate_derived_grid
  end type derived_grid_vars

  type :: grid_block
    integer, dimension(3) :: n_nodes
    integer, dimension(3) :: n_cells
    integer, dimension(3) :: n_ghost
    integer  :: n_dim
    integer  :: total_cells
    real(dp) :: total_volume
    real(dp), allocatable, dimension(:,:,:,:) ::  node_coords
    type(derived_grid_vars) :: grid_vars
  contains
    private
    procedure, public, pass :: setup   =>   allocate_grid_block
    procedure, public, pass :: destroy => deallocate_grid_block
  end type grid_block

  type grid_type
    integer  :: n_blocks
    integer  :: total_int_cells
    real(dp) :: total_int_volume
    type(grid_block), allocatable, dimension(:) :: gblock
  contains
    private
    procedure, public, pass   :: setup => init_grid_type
    procedure, public, pass :: destroy => deallocate_grid
  end type grid_type

contains

  pure function pack_cell_node_coords(idx,bnd_min,bnd_max,coords_in) result(coords_out)
  integer, dimension(3),                            intent(in)  :: idx, bnd_min, bnd_max
    real(dp), dimension( 3, bnd_min(1):bnd_max(1), &
                            bnd_min(2):bnd_max(2), &
                            bnd_min(3):bnd_max(3) ), intent(in)  :: coords_in
    real(dp), dimension(3,8)                                     :: coords_out
    integer :: i,j,k,cnt
    cnt = 0
    do k = idx(3),idx(3)+1
      do j = idx(2),idx(2)+1
        do i = idx(1),idx(1)+1
          cnt = cnt + 1
          coords_out(:,cnt) = coords_in(:,i,j,k)
        end do
      end do
    end do
  end function pack_cell_node_coords

  pure function cell_node_coords(idx,stride,bnd_min,bnd_max,coords_in) result(coords_out)
    integer, dimension(3),                            intent(in)  :: idx, stride, bnd_min, bnd_max
    real(dp), dimension( 3, bnd_min(1):bnd_max(1), &
                            bnd_min(2):bnd_max(2), &
                            bnd_min(3):bnd_max(3) ), intent(in)  :: coords_in
    real(dp), dimension(stride(1)+1,stride(2)+1,stride(3)+1,3)   :: coords_out
    integer :: i,j,k,ii,jj,kk
    kk = 0
    do k = idx(3),idx(3)+stride(3)
      kk = kk + 1
      jj = 0
      do j = idx(2),idx(2)+stride(2)
        jj = jj + 1
        ii = 0
        do i = idx(1),idx(1)+stride(1)
          ii = ii + 1
          coords_out(ii,jj,kk,1) = coords_in(1,i,j,k)
          coords_out(ii,jj,kk,2) = coords_in(2,i,j,k)
          coords_out(ii,jj,kk,3) = coords_in(3,i,j,k)
        end do
      end do
    end do
  end function cell_node_coords
  
  subroutine get_face_quad_ptrs(gblock,cell_idx,face_ids,fquad)
    use quadrature_derived_type, only : quad_ptr
    use index_conversion,        only : get_face_idx_from_id
    type(grid_block), target,                  intent(in)  :: gblock
    integer,        dimension(3),              intent(in)  :: cell_idx
    integer,        dimension(:),              intent(in)  :: face_ids
    type(quad_ptr), dimension(size(face_ids)), intent(out) :: fquad
    integer :: i, dir
    integer, dimension(3) :: face_idx
    call fquad%destroy()
    do i = 1,size(face_ids)
      call get_face_idx_from_id(cell_idx,face_ids(i),dir,face_idx)
      fquad(i)%p => gblock%grid_vars%face_quads(3)%p(face_idx(1),face_idx(2),face_idx(3))
    end do
  end subroutine get_face_quad_ptrs

  pure subroutine init_grid_type( this, n_blocks )
    use set_constants, only : zero
    class(grid_type), intent(inout) :: this
    integer, intent(in) :: n_blocks

    this%n_blocks = n_blocks
    this%total_int_cells = 0
    this%total_int_volume = zero
    allocate( this%gblock(n_blocks) )
  end subroutine init_grid_type

  pure subroutine allocate_grid_block( this, n_dim, n_nodes, n_ghost )
    use set_constants, only : zero
    integer,               intent(in)  :: n_dim
    integer, dimension(3), intent(in)  :: n_nodes, n_ghost
    class(grid_block),     intent(inout) :: this
    integer, dimension(3) :: lo, hi
    this%n_dim   = n_dim
    this%n_nodes = n_nodes
    this%n_ghost = n_ghost
    this%n_cells = 1; this%n_cells(1:n_dim) = n_nodes(1:n_dim) - 1
    lo = 1; lo(1:n_dim) = 1 - n_ghost(1:n_dim)
    hi = 1; hi(1:n_dim) = n_nodes(1:n_dim) + n_ghost(1:n_dim)
    allocate( this%node_coords( n_dim, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3) ) )
    this%node_coords = zero
    this%total_volume = zero
    this%total_cells  = product(this%n_cells)
  end subroutine allocate_grid_block

  subroutine allocate_derived_grid( this, gblock )
    use set_constants,  only : zero
    class(derived_grid_vars), target, intent(inout)   :: this
    class(grid_block),       target, intent(inout) :: gblock
    integer, dimension(3) :: lo, hi
    this%n_cells => gblock%n_cells
    this%n_ghost => gblock%n_ghost
    this%n_dim   => gblock%n_dim
    lo = 1
    lo(1:this%n_dim) = 1 - this%n_ghost(1:this%n_dim)
    hi = 1
    hi(1:this%n_dim) = this%n_cells(1:this%n_dim) + this%n_ghost(1:this%n_dim)
    allocate( this%volume(    lo(1):hi(1), lo(2):hi(2), lo(3):hi(3) ) )
    allocate( this%cell_c( 3, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3) ) )
    allocate( this%quad(      lo(1):hi(1), lo(2):hi(2), lo(3):hi(3) ) )
    this%volume = zero
    this%cell_c = zero

    lo = 1
    hi = 1
    hi(1:this%n_dim) = this%n_cells(1:this%n_dim)
    allocate( this%xi_area(   lo(1):hi(1)+1, lo(2):hi(2),   lo(3):hi(3)   ) )
    allocate( this%eta_area(  lo(1):hi(1),   lo(2):hi(2)+1, lo(3):hi(3)   ) )
    allocate( this%zeta_area( lo(1):hi(1),   lo(2):hi(2),   lo(3):hi(3)+1 ) )
    allocate( this%xi_face_quad(   lo(1):hi(1)+1, lo(2):hi(2),   lo(3):hi(3)   ) )
    allocate( this%eta_face_quad(  lo(1):hi(1),   lo(2):hi(2)+1, lo(3):hi(3)   ) )
    allocate( this%zeta_face_quad( lo(1):hi(1),   lo(2):hi(2),   lo(3):hi(3)+1 ) )
    allocate( this%xi_nv(   lo(1):hi(1)+1, lo(2):hi(2),   lo(3):hi(3)   ) )
    allocate( this%eta_nv(  lo(1):hi(1),   lo(2):hi(2)+1, lo(3):hi(3)   ) )
    allocate( this%zeta_nv( lo(1):hi(1),   lo(2):hi(2),   lo(3):hi(3)+1 ) )
    this%xi_area   = zero
    this%eta_area  = zero
    this%zeta_area = zero

    call compute_quadrature_points( gblock, gblock, [1,1,1], 2 )

    this%face_quads(1)%p => this%xi_face_quad
    this%face_quads(2)%p => this%eta_face_quad
    this%face_quads(3)%p => this%zeta_face_quad
    this%normals(1)%p    => this%xi_nv
    this%normals(2)%p    => this%eta_nv
    this%normals(3)%p    => this%zeta_nv
  end subroutine allocate_derived_grid

  pure subroutine map_cell_quad_points( dim, n_skip, quad_order, coords, mask, ref_quads, quad )
    use quadrature_derived_type, only : map_quad_ref_to_physical_1D,           &
                                        map_quad_ref_to_physical_2D,           &
                                        map_quad_ref_to_physical_3D
    integer,                                                        intent(in)    :: quad_order, dim
    integer,      dimension(3),                                     intent(in)    :: n_skip
    real(dp),     dimension(n_skip(1)+1,n_skip(1)+1,n_skip(1)+1,3), intent(in)    :: coords
    logical,      dimension(n_skip(1)+1,n_skip(1)+1,n_skip(1)+1),   intent(in)    :: mask
    type(quad_t), dimension(3),                                     intent(in)    :: ref_quads
    type(quad_t),                                                   intent(inout) :: quad
    real(dp),     dimension(product(n_skip+1)) :: Xtmp, Ytmp, Ztmp
    integer :: n_mask, n_quad

    n_mask = count(mask)

    n_quad = ref_quads(1)%n_quad

    Xtmp(1:n_mask) = pack(coords(:,:,:,1),mask)
    Ytmp(1:n_mask) = pack(coords(:,:,:,2),mask)
    Ztmp(1:n_mask) = pack(coords(:,:,:,3),mask)

    select case(dim)
    case(1)
      call map_quad_ref_to_physical_1D( Xtmp(1:n_mask), &
                                        Ytmp(1:n_mask), &
                                        Ztmp(1:n_mask), &
                                        ref_quads(1), quad )
    case(2)
      call map_quad_ref_to_physical_2D( reshape( Xtmp(1:n_mask), [n_quad, n_quad ] ), &
                                        reshape( Ytmp(1:n_mask), [n_quad, n_quad ] ), &
                                        reshape( Ztmp(1:n_mask), [n_quad, n_quad ] ), &
                                        ref_quads(2), quad )
    case(3)
      call map_quad_ref_to_physical_3D( reshape( Xtmp(1:n_mask), [n_quad, n_quad, n_quad ] ), &
                                        reshape( Ytmp(1:n_mask), [n_quad, n_quad, n_quad ] ), &
                                        reshape( Ztmp(1:n_mask), [n_quad, n_quad, n_quad ] ), &
                                        ref_quads(3), quad )
    end select
  end subroutine map_cell_quad_points


  subroutine compute_quadrature_points( gblock1, gblock, n_skip, quad_order )
    use quadrature_derived_type, only : create_quad_ref_1D,                    &
                                        create_quad_ref_2D,                    &
                                        create_quad_ref_3D
    class(grid_block),     intent(in)    :: gblock1
    class(grid_block),     intent(inout) :: gblock
    integer, dimension(3), intent(in)    :: n_skip
    integer,               intent(in)    :: quad_order
    type(quad_t), dimension(3) :: ref_quads
    logical,  dimension(n_skip(1)+1,n_skip(1)+1,n_skip(1)+1)   :: mask1, mask2, mask
    real(dp), dimension(n_skip(1)+1,n_skip(1)+1,n_skip(1)+1,3) :: coords_tmp
    integer :: i, j, k
    integer, dimension(3) :: idx, bnd_min, bnd_max, sz

    sz = n_skip + 1
    
    call create_quad_ref_1D( quad_order, ref_quads(1) )
    call create_quad_ref_2D( quad_order, ref_quads(2) )
    call create_quad_ref_3D( quad_order, ref_quads(3) )

    bnd_min = [ lbound(gblock1%node_coords,2), lbound(gblock1%node_coords,3), lbound(gblock1%node_coords,4) ]
    bnd_max = [ ubound(gblock1%node_coords,2), ubound(gblock1%node_coords,3), ubound(gblock1%node_coords,4) ]
    ! first the volume quads
    mask = .false.
    select case( gblock%n_dim )
    case(1)
      mask(:,1,1) = .true.
    case(2)
      mask(:,:,1) = .true.
    case(3)
      mask(:,:,:) = .true.
    end select
    do k = 1,gblock%n_cells(3)
      do j = 1,gblock%n_cells(2)
        do i = 1,gblock%n_cells(1)
          idx = [i,j,k]
          coords_tmp = cell_node_coords( idx, n_skip, bnd_min, bnd_max, gblock1%node_coords )
          call map_cell_quad_points( gblock%n_dim, n_skip, quad_order, coords_tmp, mask, ref_quads, gblock%grid_vars%quad(i,j,k) )
          gblock%grid_vars%volume = sum( gblock%grid_vars%quad(i,j,k)%quad_wts )
        end do
      end do
    end do

    if ( gblock%n_dim == 1) return
    ! now the face quads

    ! xi-faces
    mask1 = .false.
    mask2 = .false.
    select case( gblock%n_dim )
    case(2)
      mask1(1,:,1) = .true.
      mask2(sz(1),:,1) = .true.
    case(3)
      mask1(1,:,:) = .true.
      mask2(sz(1),:,:) = .true.
    end select
    do k = 1,gblock%n_cells(3)
      do j = 1,gblock%n_cells(2)
        do i = 1,1
          idx = [i,j,k]
          coords_tmp = cell_node_coords( idx, n_skip, bnd_min, bnd_max, gblock1%node_coords )
          call map_cell_quad_points( gblock%n_dim-1, n_skip, quad_order, coords_tmp, mask1, ref_quads, gblock%grid_vars%xi_face_quad(i,j,k) )
          gblock%grid_vars%xi_area = sum( gblock%grid_vars%xi_face_quad(i,j,k)%quad_wts )
        end do
        do i = 1,gblock%n_cells(1)
          idx = [i,j,k]
          coords_tmp = cell_node_coords( idx, n_skip, bnd_min, bnd_max, gblock1%node_coords )
          call map_cell_quad_points( gblock%n_dim-1, n_skip, quad_order, coords_tmp, mask2, ref_quads, gblock%grid_vars%xi_face_quad(i,j,k) )
          gblock%grid_vars%xi_area = sum( gblock%grid_vars%xi_face_quad(i,j,k)%quad_wts )
        end do
      end do
    end do

    ! eta-faces
    mask1 = .false.
    mask2 = .false.
    select case( gblock%n_dim )
    case(2)
      mask1(:,1,1) = .true.
      mask2(:,sz(2),1) = .true.
    case(3)
      mask1(:,1,:) = .true.
      mask2(:,sz(2),:) = .true.
    end select
    do k = 1,gblock%n_cells(3)
      do j = 1,1
        do i = 1,gblock%n_cells(1)
          idx = [i,j,k]
          coords_tmp = cell_node_coords( idx, n_skip, bnd_min, bnd_max, gblock1%node_coords )
          call map_cell_quad_points( gblock%n_dim-1, n_skip, quad_order, coords_tmp, mask1, ref_quads, gblock%grid_vars%eta_face_quad(i,j,k) )
          gblock%grid_vars%eta_area = sum( gblock%grid_vars%eta_face_quad(i,j,k)%quad_wts )
        end do
      end do
      do j = 1,gblock%n_cells(2)
        do i = 1,gblock%n_cells(1)
          idx = [i,j,k]
          coords_tmp = cell_node_coords( idx, n_skip, bnd_min, bnd_max, gblock1%node_coords )
          call map_cell_quad_points( gblock%n_dim-1, n_skip, quad_order, coords_tmp, mask2, ref_quads, gblock%grid_vars%eta_face_quad(i,j,k) )
          gblock%grid_vars%eta_area = sum( gblock%grid_vars%eta_face_quad(i,j,k)%quad_wts )
        end do
      end do
    end do

    if ( gblock%n_dim == 2) return

    ! zeta-faces
    mask1 = .false.; mask1(:,:,1) = .true.
    mask2 = .false.; mask2(:,:,sz(3)) = .true.
    do k = 1,1
      do j = 1,gblock%n_cells(2)
        do i = 1,gblock%n_cells(1)
          idx = [i,j,k]
          coords_tmp = cell_node_coords( idx, n_skip, bnd_min, bnd_max, gblock1%node_coords )
          call map_cell_quad_points( gblock%n_dim-1, n_skip, quad_order, coords_tmp, mask1, ref_quads, gblock%grid_vars%zeta_face_quad(i,j,k) )
          gblock%grid_vars%zeta_area = sum( gblock%grid_vars%zeta_face_quad(i,j,k)%quad_wts )
        end do
      end do
    end do
    do k = 1,gblock%n_cells(3)
      do j = 1,gblock%n_cells(2)
        do i = 1,gblock%n_cells(1)
          idx = [i,j,k]
          coords_tmp = cell_node_coords( idx, n_skip, bnd_min, bnd_max, gblock1%node_coords )
          call map_cell_quad_points( gblock%n_dim-1, n_skip, quad_order, coords_tmp, mask2, ref_quads, gblock%grid_vars%zeta_face_quad(i,j,k) )
          gblock%grid_vars%zeta_area = sum( gblock%grid_vars%zeta_face_quad(i,j,k)%quad_wts )
        end do
      end do
    end do

    call ref_quads%destroy()

  end subroutine compute_quadrature_points

  pure elemental subroutine deallocate_grid(this)
    use set_constants, only : zero
    class(grid_type), intent(inout) :: this
    if ( allocated(this%gblock) ) then
      call this%gblock%destroy()
      deallocate( this%gblock)
    end if
    this%n_blocks = 0
    this%total_int_cells = 0
    this%total_int_volume = zero
  end subroutine deallocate_grid

  pure elemental subroutine deallocate_grid_block( this )
    use set_constants, only : zero
    class(grid_block), intent(inout) :: this
    call this%grid_vars%destroy()
    if (allocated( this%node_coords ) ) deallocate( this%node_coords )
    this%total_volume = zero
    this%total_cells  = 0
    this%n_dim        = 0
    this%n_nodes      = 0
    this%n_cells      = 0
  end subroutine deallocate_grid_block

  pure elemental subroutine deallocate_derived_grid( this )
    class(derived_grid_vars), intent(inout) :: this
    call this%face_quads(1)%destroy()
    call this%face_quads(2)%destroy()
    call this%face_quads(3)%destroy()
    call this%normals(1)%destroy()
    call this%normals(2)%destroy()
    call this%normals(3)%destroy()
    if ( allocated(this%cell_c)    ) deallocate( this%cell_c    )
    if ( allocated(this%volume)    ) deallocate( this%volume    )
    if ( allocated(this%xi_area)   ) deallocate( this%xi_area   )
    if ( allocated(this%eta_area)  ) deallocate( this%eta_area  )
    if ( allocated(this%zeta_area) ) deallocate( this%zeta_area )
    if ( allocated(this%quad) ) then
      call this%quad%destroy()
      deallocate( this%quad )
    end if
    if ( allocated(this%xi_face_quad) ) then
      call this%xi_face_quad%destroy()
      deallocate( this%xi_face_quad )
    end if
    if ( allocated(this%eta_face_quad) ) then
      call this%eta_face_quad%destroy()
      deallocate( this%eta_face_quad )
    end if
    if ( allocated(this%zeta_face_quad) ) then
      call this%zeta_face_quad%destroy()
      deallocate( this%zeta_face_quad )
    end if
    if ( allocated(this%xi_nv) ) then
      call this%xi_nv%destroy()
      deallocate(this%xi_nv)
    end if
    if ( allocated(this%eta_nv) ) then
      call this%eta_nv%destroy()
      deallocate(this%eta_nv)
    end if
    if ( allocated(this%zeta_nv) ) then
      call this%zeta_nv%destroy()
      deallocate(this%zeta_nv)
    end if
  end subroutine deallocate_derived_grid

end module grid_derived_type

module monomial_basis_derived_type
  implicit none
  private
  public :: monomial_basis_t
  type :: monomial_basis_t
    private
    integer, public :: total_degree
    integer, public :: n_dim
    integer, public :: n_terms
    integer, public, allocatable, dimension(:)   :: idx
    integer, public, allocatable, dimension(:,:) :: exponents
  contains
    private
    procedure, public, pass :: eval  => evaluate_monomial
    procedure, public, pass :: deval => evaluate_monomial_derivative
    procedure, public, pass :: destroy => destroy_monomial_basis_t
  end type monomial_basis_t

  interface monomial_basis_t
    procedure constructor
  end interface monomial_basis_t

contains

  pure function constructor( total_degree, n_dim ) result(this)
    use combinatorics, only : nchoosek, get_exponents
    integer, intent(in) :: total_degree, n_dim
    type(monomial_basis_t) :: this

    call this%destroy()

    this%total_degree  = total_degree
    this%n_dim   = n_dim
    this%n_terms = nchoosek( n_dim + total_degree, total_degree )
    allocate( this%exponents( this%n_dim, this%n_terms ) )
    allocate( this%idx(this%total_degree+1) )
    call get_exponents( this%n_dim, this%total_degree, this%n_terms, this%exponents, this%idx )
  end function constructor
  
  pure subroutine destroy_monomial_basis_t(this)
    class(monomial_basis_t), intent(inout) :: this
    if ( allocated(this%exponents) ) deallocate( this%exponents )
    if ( allocated(this%idx) )       deallocate( this%idx )
  end subroutine destroy_monomial_basis_t

  pure subroutine evaluate_monomial(this,term,x,val,coef)
    use set_precision, only : dp
    use set_constants, only : one
    class(monomial_basis_t),         intent(in)  :: this
    integer,                         intent(in)  :: term
    real(dp), dimension(this%n_dim), intent(in)  :: x
    real(dp),                        intent(out) :: val
    integer,                         intent(out) :: coef
    integer :: d, i
    val  = one
    coef = 1
    do d = 1,this%n_dim
      do i = this%exponents(d,term),1,-1
        val  = val * x(d)
        coef = coef * i
      end do
    end do
  end subroutine evaluate_monomial

  pure subroutine evaluate_monomial_derivative(this,term,x,order,dval,dcoef,coef)
    use set_precision, only : dp
    use set_constants, only : zero, one
    class(monomial_basis_t),         intent(in)  :: this
    integer,                         intent(in)  :: term
    real(dp), dimension(this%n_dim), intent(in)  :: x
    integer,  dimension(this%n_dim), intent(in)  :: order
    real(dp),                        intent(out) :: dval
    integer,                         intent(out) :: dcoef, coef
    integer :: d, i
    
    dcoef = 1
    coef  = 1
    dval  = zero
    if ( any( this%exponents(:,term)-order < 0 ) ) return

    dval  = one
    do d = 1,this%n_dim
      do i = this%exponents(d,term),this%exponents(d,term)-order(d)+1,-1
        dcoef = dcoef * i
      end do
      do i = this%exponents(d,term)-order(d),1,-1
        dval  = dval * x(d)
        coef = coef * i
      end do
    end do
  end subroutine evaluate_monomial_derivative

end module monomial_basis_derived_type

module zero_mean_basis_derived_type
  use set_precision, only : dp
  use quadrature_derived_type, only : quad_t
  use monomial_basis_derived_type, only : monomial_basis_t
  implicit none
  private
  public :: zero_mean_basis_t
  type :: zero_mean_basis_t
    private
    type(monomial_basis_t), pointer :: mono_basis => null()
    real(dp), public, allocatable, dimension(:)   :: moments
    real(dp), public, allocatable, dimension(:)   :: x_ref
    real(dp), public, allocatable, dimension(:)   :: h_ref
    integer, public, pointer :: total_degree
    integer, public, pointer :: n_dim
    integer, public, pointer :: n_terms
    integer, public, pointer, dimension(:)   :: idx
    integer, public, pointer, dimension(:,:)   :: exponents
  contains
    private
    procedure, pass :: compute_grid_moment
    procedure, pass :: transform
    
    procedure, nopass :: length_scale  => get_length_scale_vector
    procedure, public, pass :: eval  => evaluate_basis
    procedure, public, pass :: deval => evaluate_basis_derivative
    procedure, public, pass :: scaled_basis_derivative
    procedure, public, pass :: scaled_basis_derivatives
    procedure, public, pass :: destroy => destroy_zero_mean_basis_t

  end type zero_mean_basis_t

  interface zero_mean_basis_t
    procedure constructor
  end interface

contains

function constructor( mono_basis, quad, h_ref ) result(this)
  use monomial_basis_derived_type, only : monomial_basis_t
  type(monomial_basis_t), target, intent(in) :: mono_basis
  type(quad_t), intent(in)                   :: quad
  real(dp),     intent(in), dimension(:) :: h_ref
  type(zero_mean_basis_t) :: this
  integer :: n

  call this%destroy()
  this%mono_basis => mono_basis
  this%total_degree => mono_basis%total_degree
  this%n_dim        => mono_basis%n_dim
  this%n_terms      => mono_basis%n_terms
  this%idx          => mono_basis%idx
  this%exponents    => mono_basis%exponents

    allocate( this%moments( this%n_terms ) )
  allocate( this%x_ref( this%n_dim ) )
  allocate( this%h_ref( this%n_dim ) )
  this%h_ref = h_ref
  this%x_ref = quad%integrate( this%n_dim, quad%quad_pts(1:this%n_dim,:) ) / sum( quad%quad_wts )
  do n = 1,this%n_terms
    this%moments(n) = this%compute_grid_moment(n,quad)
  end do
end function constructor

pure subroutine destroy_zero_mean_basis_t(this)
  class(zero_mean_basis_t), intent(inout) :: this
  if ( allocated(this%moments) ) deallocate( this%moments )
  if ( allocated(this%x_ref)   ) deallocate( this%x_ref   )
  if ( allocated(this%h_ref)   ) deallocate( this%h_ref   )
end subroutine destroy_zero_mean_basis_t

pure function transform(this,x) result(x_bar)
  class(zero_mean_basis_t), intent(in) :: this
  real(dp), dimension(this%n_dim), intent(in) :: x
  real(dp), dimension(this%n_dim) :: x_bar
  x_bar = (x-this%x_ref)/this%h_ref
end function transform

pure elemental function compute_grid_moment(this,term,quad) result(moment)
  class(zero_mean_basis_t), intent(in) :: this
  integer,                  intent(in) :: term
  type(quad_t),             intent(in) :: quad
  real(dp), dimension(quad%n_quad) :: tmp
  real(dp), dimension(size(quad%quad_pts,1)) :: xtmp
  real(dp) :: moment
  integer :: q, coef
    do q = 1,quad%n_quad
      xtmp = this%transform( quad%quad_pts(1:this%n_dim,q) )
      call this%mono_basis%eval( term, xtmp, tmp(q), coef )
    end do
    moment = quad%integrate(tmp) / sum( quad%quad_wts )
end function compute_grid_moment

pure function evaluate_basis(this,term,point) result(B)
  use set_constants, only : one
  class(zero_mean_basis_t), intent(in) :: this
  integer,                  intent(in) :: term
  real(dp), dimension(this%n_dim), intent(in) :: point
  real(dp) :: B
  integer :: coef
  B = one
  if (term == 1) return
  call this%mono_basis%eval(term,point,B,coef)
  B = B - this%moments(term)
end function evaluate_basis

pure function evaluate_basis_derivative(this,term,point,order) result(dB)
  use set_constants, only : zero, one
  class(zero_mean_basis_t), intent(in) :: this
  integer,                  intent(in) :: term
  integer,  dimension(this%n_dim), intent(in) :: order
  real(dp), dimension(this%n_dim), intent(in) :: point
  real(dp) :: dB
  integer :: dcoef,coef

  if (all(order==0)) then
    dB =  this%eval(term,point)
    return
  end if

  dB = zero
  if (term==1) return

  call this%mono_basis%deval( term, point,order,dB,dcoef,coef)
  dB = dB * dcoef / product( this%h_ref ** order )
end function evaluate_basis_derivative

pure function get_length_scale_vector( order, scale ) result(L)
  use set_constants, only : one
  integer,  dimension(:),           intent(in) :: order
  real(dp), dimension(size(order)), intent(in) :: scale
  real(dp)                                     :: L
  integer :: n, d, den

  L   = one
  den = 1
  do d = 1,size(order)
    do n = order(d),1,-1
      L   = L * scale(d)
      den = den * n
    end do
  end do
  L = L / real(den,dp)

end function get_length_scale_vector

pure function scaled_basis_derivative( this, term_idx, diff_idx, point, scale ) result(derivative)
  class(zero_mean_basis_t),        intent(in) :: this
  integer,                         intent(in) :: term_idx, diff_idx
  real(dp), dimension(this%n_dim), intent(in) :: point
  real(dp), dimension(this%n_dim), intent(in) :: scale
  real(dp)                                    :: derivative
  real(dp) :: L
  derivative = this%deval( term_idx, point, this%exponents(:,diff_idx) )
  L = this%length_scale( this%exponents(:,diff_idx), scale )
  derivative = derivative * L
end function scaled_basis_derivative

pure function scaled_basis_derivatives( this, term_start, term_end, point, scale ) result(derivatives)
  use set_constants, only : zero
  class(zero_mean_basis_t),        intent(in) :: this
  integer,                         intent(in) :: term_start, term_end
  real(dp), dimension(this%n_dim), intent(in) :: point
  real(dp), dimension(this%n_dim), intent(in) :: scale
  real(dp), dimension(term_end, term_end - term_start) :: derivatives
  integer :: i, j
  derivatives = zero
  ! outer loop over basis functions
  do j = term_start+1,term_end
    ! inner loop over derivatives
    do i = 1,j
      derivatives(i,j-term_start) = this%scaled_basis_derivative(j,i,point,scale);
    end do
  end do
end function scaled_basis_derivatives

end module zero_mean_basis_derived_type
module var_rec_derived_type
  use set_precision,                only : dp
  use monomial_basis_derived_type,  only : monomial_basis_t
  use zero_mean_basis_derived_type, only : zero_mean_basis_t
  use quadrature_derived_type,      only : quad_t
  implicit none
  private
  public :: var_rec_t
  type :: var_rec_t
    ! private
    integer :: n_vars
    integer :: n_interior
    integer :: n_boundary
    integer :: self_idx
    integer :: self_block
    integer, dimension(:), allocatable      :: nbor_idx, face_id
    type(zero_mean_basis_t)                 :: basis
    ! type(quad_t), pointer                   :: quad
    real(dp), dimension(:,:),   allocatable :: coefs
    real(dp), dimension(:,:),   allocatable :: A_inv
    real(dp), dimension(:,:),   allocatable :: A_inv_D
    real(dp), dimension(:,:,:), allocatable :: A_inv_B
    real(dp), dimension(:,:,:), allocatable :: A_inv_C
    real(dp), dimension(:,:),   allocatable :: A, D
    real(dp), dimension(:,:,:), allocatable :: B, C
  contains
    private
    ! procedure, public, pass :: setup => setup_reconstruction
    ! procedure, public, pass :: eval  => eval_reconstruction
    procedure, public, pass :: allocate_reconstruction_matrices
    procedure, public, pass :: get_nbor_contribution
    ! procedure, pass :: get_LHS_interior
    ! procedure, pass :: get_RHS_interior
    procedure, public, pass :: get_RHS_contribution_update
    procedure, public, pass :: get_LHS_contribution_update
    ! procedure, pass :: get_LHS_dirichlet
    ! procedure, pass :: get_RHS_dirichlet
    procedure, public, pass :: destroy => destroy_var_rec_t
  end type var_rec_t

  interface var_rec_t
    procedure constructor
  end interface var_rec_t

contains

  function constructor( self_block, self_idx, n_interior, nbor_idx, face_id, n_vars, mono_basis, quad, h_ref ) result(this)
    use set_constants, only : zero
    integer,                    intent(in) :: self_block, self_idx, n_interior
    integer, dimension(:),      intent(in) :: nbor_idx, face_id
    integer,                    intent(in) :: n_vars
    type(monomial_basis_t), target, intent(in) :: mono_basis
    type(quad_t),               intent(in) :: quad
    real(dp), dimension(mono_basis%n_dim), intent(in) :: h_ref
    type(var_rec_t)                        :: this

    call this%destroy()
    this%self_block = self_block
    this%self_idx   = self_idx
    this%n_interior = n_interior
    this%n_boundary = size(nbor_idx) - n_interior
    this%n_vars     = n_vars
    this%basis      = zero_mean_basis_t(mono_basis,quad,h_ref)
    
    allocate( this%nbor_idx( size(nbor_idx) ) )
    this%nbor_idx   = nbor_idx

    allocate( this%face_id( size(face_id) ) )
    this%face_id    = face_id
    
    allocate( this%coefs(this%basis%n_terms,this%n_vars) )
    this%coefs      = zero
  end function constructor

  pure elemental subroutine destroy_var_rec_t(this)
    class(var_rec_t), intent(inout) :: this
    call this%basis%destroy()
    if ( allocated(this%nbor_idx) ) deallocate( this%nbor_idx )
    if ( allocated(this%face_id)  ) deallocate( this%face_id )
    if ( allocated(this%coefs  )  ) deallocate( this%coefs   )
    if ( allocated(this%A_inv  )  ) deallocate( this%A_inv   )
    if ( allocated(this%A_inv_B)  ) deallocate( this%A_inv_B )
    if ( allocated(this%A_inv_C)  ) deallocate( this%A_inv_C )
    if ( allocated(this%A_inv_D)  ) deallocate( this%A_inv_D )
    if ( allocated(this%A      )  ) deallocate( this%A       )
    if ( allocated(this%B      )  ) deallocate( this%B       )
    if ( allocated(this%C      )  ) deallocate( this%C       )
    if ( allocated(this%D      )  ) deallocate( this%D       )
  end subroutine destroy_var_rec_t

  subroutine allocate_reconstruction_matrices( this, term_start, term_end )
    class(var_rec_t),                             intent(inout) :: this
    integer,                                      intent(in)    :: term_start, term_end
    if ( allocated(this%A_inv  ) ) deallocate( this%A_inv   )
    if ( allocated(this%A_inv_B) ) deallocate( this%A_inv_B )
    if ( allocated(this%A_inv_C) ) deallocate( this%A_inv_C )
    if ( allocated(this%A_inv_D) ) deallocate( this%A_inv_D )
    if ( allocated(this%A      ) ) deallocate( this%A       )
    if ( allocated(this%B      ) ) deallocate( this%B       )
    if ( allocated(this%C      ) ) deallocate( this%C       )
    if ( allocated(this%D      ) ) deallocate( this%D       )
    allocate( this%A_inv(   term_end-term_start,term_end-term_start ) )
    allocate( this%A_inv_B( term_end-term_start,term_end-term_start, this%n_interior ) )
    allocate( this%A_inv_C( term_end-term_start, term_start, this%n_interior ) )
    allocate( this%A_inv_D( term_end-term_start, term_start ) )
    allocate( this%A(       term_end-term_start,term_end-term_start ) )
    allocate( this%B( term_end-term_start,term_end-term_start, this%n_interior ) )
    allocate( this%C( term_end-term_start, term_start, this%n_interior ) )
    allocate( this%D( term_end-term_start, term_start ) )
  end subroutine allocate_reconstruction_matrices

  pure subroutine get_nbor_contribution( this, nbor, fquad, term_start, term_end, A, B, D, C )
    use set_constants, only : zero, one
    class(var_rec_t), intent(in) :: this
    class(var_rec_t), intent(in) :: nbor
    class(quad_t),    intent(in) :: fquad
    integer,          intent(in) :: term_start, term_end
    real(dp), dimension(term_end - term_start, term_end - term_start), intent(out) :: A, B
    real(dp), dimension(term_end - term_start, term_start),            intent(out) :: D, C
    real(dp), dimension(term_end,term_end-term_start) :: d_basis_i
    real(dp), dimension(term_end,max(term_end-term_start,term_start)) :: d_basis_j
    integer :: q, l, m
    real(dp), dimension(this%basis%n_dim) :: dij
    real(dp) :: xdij_mag

    A = zero
    B = zero
    C = zero
    D = zero
    dij = abs( this%basis%x_ref - nbor%basis%x_ref )
    xdij_mag = one/norm2(dij)
    do q = 1,fquad%n_quad
      ! LHS
      d_basis_i = this%basis%scaled_basis_derivatives( term_start, term_end, fquad%quad_pts(:,q), dij )
      d_basis_j(:,1:term_end-term_start) = nbor%basis%scaled_basis_derivatives( term_start, term_end, fquad%quad_pts(:,q), dij )
      do m = 1,term_end-term_start
        do l = 1,term_end-term_start
          A(l,m) = A(l,m) + fquad%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_i(:,m) )
          B(l,m) = B(l,m) + fquad%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_j(:,m) )
        end do
      end do

      ! RHS
      d_basis_j = zero
      d_basis_j(1:term_start,1:term_start) = nbor%basis%scaled_basis_derivatives( 0, term_start, fquad%quad_pts(:,q), dij )
      do m = 1,term_start
        do l = 1,term_end-term_start
          D(l,m) = D(l,m) + fquad%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_i(:,m) )
          C(l,m) = C(l,m) + fquad%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_j(:,m) )
        end do
      end do
    end do
  end subroutine get_nbor_contribution

  pure function get_RHS_contribution_update( this, nbor, nbor_id, term_start, term_end, var_idx ) result(update)
    use set_constants, only : zero, one
    class(var_rec_t),    intent(in) :: this
    class(var_rec_t),    intent(in) :: nbor
    integer,             intent(in) :: nbor_id, term_start, term_end
    integer,             intent(in) :: var_idx
    real(dp), dimension( term_end - term_start ) :: update
    update = zero
    update = update + matmul( this%A_inv_C(:,:,nbor_id), nbor%coefs(1:term_start,var_idx) )
    update = update + matmul( this%A_inv_D, this%coefs(1:term_start,var_idx) )
  end function get_RHS_contribution_update

  pure function get_LHS_contribution_update( this, nbor, nbor_id, term_start, term_end, var_idx ) result(update)
    use set_constants, only : zero, one
    class(var_rec_t),    intent(in) :: this
    class(var_rec_t),    intent(in) :: nbor
    integer,             intent(in) :: nbor_id, term_start, term_end
    integer,             intent(in) :: var_idx
    real(dp), dimension( term_end - term_start ) :: update
    update = matmul( this%A_inv_B(:,:,nbor_id), nbor%coefs(term_start+1:term_end,var_idx) )
  end function get_LHS_contribution_update

  ! pure subroutine get_RHS_interior_update( this, nbors, term_start, term_end, var_idx, b )
  !   use set_constants, only : zero, one
  !   class(var_rec_t),                             intent(in) :: this
  !   class(var_rec_t), dimension(this%n_interior), intent(in) :: nbors
  !   integer,                                      intent(in) :: term_start, term_end
  !   integer,          dimension(:),               intent(in) :: var_idx
  !   real(dp), dimension( term_end - term_start,size(var_idx) ), intent(out) :: b
  !   integer :: v, j

  !   b = zero
  !   do v = 1,size(var_idx)
  !     do j = 1,this%n_interior
  !       b(:,v) = b(:,v) + matmul( this%A_inv_C(:,:,j), nbors(j)%coefs(term_start+1:term_end,var_idx(v)) )
  !     end do
  !     b(:,v) = b(:,v) + matmul( this%A_inv_D, this%coefs(term_start+1:term_end,var_idx(v)) )
  !   end do
  ! end subroutine get_RHS_interior_update

  subroutine get_LHS_dirichlet( this, fquad, term_start, term_end, A )
    use set_constants, only : zero, one
    class(var_rec_t),                             intent(in) :: this
    class(quad_t),                                intent(in) :: fquad
    integer,                                      intent(in) :: term_start, term_end
    real(dp), dimension(term_end - term_start, term_end - term_start), intent(out) :: A
    real(dp), dimension(term_end-term_start) :: basis
    integer :: q, l, m
    real(dp), dimension(this%basis%n_dim) :: dib, face_x_ref
    real(dp) :: xdib_mag

    A = zero
    face_x_ref = fquad%integrate(this%basis%n_dim,fquad%quad_pts) / sum( fquad%quad_wts )

    dib = abs( this%basis%x_ref - face_x_ref )
    xdib_mag = one/norm2(dib)
    do q = 1,fquad%n_quad
      do m = 1,term_end-term_start
        basis(m) = this%basis%scaled_basis_derivative(m,1,fquad%quad_pts(:,q),dib)
      end do
      do m = 1,term_end-term_start
        do l = 1,term_end-term_start
          A(l,m)   = A(l,m)   + fquad%quad_wts(q) * xdib_mag * basis(l) * basis(m)
        end do
      end do
    end do

  end subroutine get_LHS_dirichlet

  pure subroutine get_RHS_dirichlet( this, fquad, term_start, term_end, var_idx, bc_eval, b )
    use set_constants, only : zero, one
    class(var_rec_t),                                         intent(in)  :: this
    class(quad_t),                                            intent(in)  :: fquad
    integer,                                                  intent(in)  :: term_start, term_end
    integer,  dimension(:),                                   intent(in)  :: var_idx
    real(dp), dimension(size(var_idx),fquad%n_quad),          intent(in)  :: bc_eval
    real(dp), dimension(term_end - term_start,size(var_idx)), intent(out) :: b
    real(dp), dimension(term_end-term_start) :: basis
    integer :: v, q, l, m
    real(dp), dimension(this%basis%n_dim) :: dib, face_x_ref
    real(dp) :: xdib_mag

    face_x_ref = fquad%integrate(this%basis%n_dim,fquad%quad_pts) / sum( fquad%quad_wts )

    dib = abs( this%basis%x_ref - face_x_ref )
    xdib_mag = one/norm2(dib)
    do q = 1,fquad%n_quad
      do l = 1,term_end-term_start
        basis(l) = this%basis%eval( l,fquad%quad_pts(:,q) )
      end do
      do v = 1,size(var_idx)
        do l = 1,term_end-term_start
          b(l,v) = b(l,v) + fquad%quad_wts(q) * xdib_mag * basis(l) * ( bc_eval(v,q) - this%coefs(1,var_idx(v)) )
        end do
      end do
    end do

  end subroutine get_RHS_dirichlet

end module var_rec_derived_type

module var_rec_holder_type
  use set_precision,               only : dp
  use monomial_basis_derived_type, only : monomial_basis_t
  use var_rec_derived_type,        only : var_rec_t
  
  implicit none
  private
  public :: reconstruction_holder
  type :: reconstruction_holder
    private
    integer :: n_cells
    integer :: degree
    integer :: term_start, term_end
    type(var_rec_t), dimension(:), allocatable :: rec
    type(monomial_basis_t) :: monomial_basis
  contains
    private
    procedure, public, pass :: create      => create_reconstruction_holder
    procedure, public, pass :: setup       => setup_reconstruction_holder
    procedure, public, pass :: setup_basic => setup_reconstruction_holder_basic
    procedure, public, pass :: destroy     => destroy_reconstruction_holder
    procedure, public, pass :: solve       => perform_iterative_reconstruction_SOR
    procedure,         pass :: SOR_iteration, SOR_cell_residual, SOR_residual_norm
    procedure,         pass :: set_cell_avgs
  end type reconstruction_holder

  abstract interface
    pure function spatial_function(n_var,x) result(val)
      use set_precision, only : dp
      integer,                intent(in) :: n_var
      real(dp), dimension(:), intent(in) :: x
      real(dp), dimension(n_var)         :: val
    end function spatial_function
  end interface
contains
  subroutine create_reconstruction_holder( this, block_num, degree, n_vars, gblock )
    use math, only : maximal_extents
    use grid_derived_type, only : grid_block, pack_cell_node_coords
    use index_conversion,  only : cell_face_nbors
    class(reconstruction_holder), intent(inout) :: this
    integer,                     intent(in)    :: block_num, degree, n_vars
    type(grid_block),            intent(in)    :: gblock
    real(dp), dimension(gblock%n_dim) :: h_ref
    real(dp), dimension(3,8) :: nodes
    integer, dimension(3) :: lo, hi, bnd_min, bnd_max
    integer, dimension(2*gblock%n_dim) :: nbor_cell_idx, nbor_face_id
    integer :: i, j, k, cnt, n_int

    call this%destroy()
    this%degree         = degree
    this%term_start     = 0
    this%term_end       = 0
    this%monomial_basis = monomial_basis_t(this%degree,gblock%n_dim)
    this%n_cells = product( gblock%n_cells )
    lo = [1,1,1]        - gblock%n_ghost
    hi = gblock%n_nodes + gblock%n_ghost
    bnd_min = [1,1,1]
    bnd_max = gblock%n_cells
    allocate( this%rec( this%n_cells ) )
    cnt = 0
    do k = 1,gblock%n_cells(3)
      do j = 1,gblock%n_cells(2)
        do i = 1,gblock%n_cells(1)
          cnt = cnt + 1
          nodes = pack_cell_node_coords( [i,j,k],lo,hi,gblock%node_coords)
          h_ref = maximal_extents( gblock%n_dim, 8, nodes(1:gblock%n_dim,:) )
          call cell_face_nbors( gblock%n_dim, cnt, bnd_min(1:gblock%n_dim), &
                                bnd_max(1:gblock%n_dim), nbor_cell_idx, nbor_face_id, n_int )
          this%rec(cnt) = var_rec_t(block_num, cnt, n_int, nbor_cell_idx, nbor_face_id, n_vars, this%monomial_basis, gblock%grid_vars%quad(i,j,k), h_ref)
        end do
      end do
    end do

  end subroutine create_reconstruction_holder

  subroutine setup_reconstruction_holder( this, term_start, term_end, gblock, eval_fun )
    use set_constants,           only : zero
    use index_conversion,        only : get_face_idx_from_id
    use math,                    only : mat_inv
    use grid_derived_type,       only : grid_block
    use var_rec_derived_type,    only : var_rec_t
    use quadrature_derived_type, only : quad_t
    class(reconstruction_holder), target, intent(inout) :: this
    integer,                              intent(in)    :: term_start, term_end
    type(grid_block),             target, intent(in)    :: gblock
    procedure(spatial_function), optional               :: eval_fun
    type(var_rec_t), pointer :: nbor  => null()
    type(quad_t),    pointer :: fquad => null()
    real(dp), dimension(term_end - term_start, term_end - term_start )                :: dA
    real(dp), dimension(term_end - term_start,            term_start )                :: dD
    ! real(dp), dimension(term_end - term_start, term_end - term_start, 2*gblock%n_dim) :: B
    ! real(dp), dimension(term_end - term_start,            term_start, 2*gblock%n_dim) :: C
    integer, dimension(3) :: face_idx
    integer :: i, j, k, n, cnt, dir
    this%term_start = term_start
    this%term_end   = term_end

    if (present(eval_fun)) call this%set_cell_avgs(gblock,eval_fun)

    cnt = 0
    do k = 1,gblock%n_cells(3)
      do j = 1,gblock%n_cells(2)
        do i = 1,gblock%n_cells(1)
          cnt = cnt + 1
          call this%rec(cnt)%allocate_reconstruction_matrices( term_start, term_end )
          this%rec(cnt)%A = zero; this%rec(cnt)%B = zero; this%rec(cnt)%D = zero; this%rec(cnt)%C = zero
          do n = 1,this%rec(cnt)%n_interior
            call get_face_idx_from_id([i,j,k],this%rec(cnt)%face_id(n),dir,face_idx)
            fquad => gblock%grid_vars%face_quads(dir)%p(face_idx(1),face_idx(2),face_idx(3))
            nbor  => this%rec( this%rec(cnt)%nbor_idx(n) )
            call this%rec(cnt)%get_nbor_contribution( nbor, fquad, term_start, term_end, dA, this%rec(cnt)%B(:,:,n), dD, this%rec(cnt)%C(:,:,n) )
            this%rec(cnt)%A = this%rec(cnt)%A + dA
            this%rec(cnt)%D = this%rec(cnt)%D + dD
          end do
          call mat_inv( this%rec(cnt)%A, this%rec(cnt)%A_inv, term_end - term_start )
          this%rec(cnt)%A_inv_D = matmul(this%rec(cnt)%A_inv,this%rec(cnt)%D)
          do n = 1,this%rec(cnt)%n_interior
            this%rec(cnt)%A_inv_B(:,:,n) = matmul( this%rec(cnt)%A_inv, this%rec(cnt)%B(:,:,n) )
            this%rec(cnt)%A_inv_C(:,:,n) = matmul( this%rec(cnt)%A_inv, this%rec(cnt)%C(:,:,n) )
          end do
        end do
      end do
    end do
  end subroutine setup_reconstruction_holder

  subroutine setup_reconstruction_holder_basic( this, gblock, eval_fun )
    use grid_derived_type, only : grid_block
    class(reconstruction_holder), intent(inout) :: this
    type(grid_block),             intent(in)    :: gblock
    procedure(spatial_function), optional       :: eval_fun
    integer :: term_start, term_end
    term_start = 2
    term_end   = this%monomial_basis%n_terms
    if (present(eval_fun)) then
      call this%setup(term_start,term_end,gblock,eval_fun=eval_fun)
    else
      call this%setup(term_start,term_end,gblock)
    end if
  end subroutine setup_reconstruction_holder_basic

  pure elemental subroutine destroy_reconstruction_holder( this )
    class(reconstruction_holder), intent(inout) :: this
    if ( allocated(this%rec) ) then
      call this%rec%destroy()
      deallocate( this%rec )
    end if
    call this%monomial_basis%destroy()
    this%n_cells    = 0
    this%degree     = 0
    this%term_start = 0
    this%term_end   = 0
  end subroutine destroy_reconstruction_holder

  subroutine set_cell_avgs(this,gblock,eval_fun)
    use set_constants,           only : zero
    use grid_derived_type,       only : grid_block
    class(reconstruction_holder), intent(inout) :: this
    type(grid_block),             intent(in)    :: gblock
    procedure(spatial_function)                 :: eval_fun
    real(dp), dimension(:,:), allocatable :: tmp_val
    real(dp), dimension(3) :: pt
    integer :: i, j, k, n, cnt, n_quad, n_var
    n_quad = maxval( gblock%grid_vars%quad%n_quad )
    n_var  = maxval( this%rec%n_vars )
    allocate( tmp_val(n_var,n_quad) )

    !  pure function integrate_vector( this, neq, f ) result( integral )
    ! use set_precision, only : dp
    ! class(quad_t),                        intent(in) :: this
    ! integer,                              intent(in) :: neq
    ! real(dp), dimension(neq,this%n_quad), intent(in) :: f
    cnt = 0
    do k = 1,gblock%n_cells(3)
      do j = 1,gblock%n_cells(2)
        do i = 1,gblock%n_cells(1)
          associate( quad   => gblock%grid_vars%quad(i,j,k), &
                     volume => gblock%grid_vars%volume(i,j,k) )
            cnt = cnt + 1
            n_var  = this%rec(cnt)%n_vars
            n_quad = quad%n_quad
            tmp_val = zero
            do n = 1,n_quad
              pt = quad%quad_pts(:,n)
              tmp_val(1:n_var,n) = eval_fun(n_var,pt)
            end do
            this%rec(cnt)%coefs(1,1:n_var) = quad%integrate(n_var,tmp_val(1:n_var,1:n_quad)) / volume
          end associate
        end do
      end do
    end do
  end subroutine set_cell_avgs

  subroutine perform_iterative_reconstruction_SOR(this,var_idx,omega,tol,n_iter,converged,residual)
    use set_constants, only : zero, one
    class(reconstruction_holder),       intent(inout) :: this
    integer,  dimension(:),             intent(in)    :: var_idx
    real(dp), optional,                 intent(in)    :: omega, tol
    integer,  optional,                 intent(in)    :: n_iter
    logical,  optional,                 intent(out)   :: converged
    real(dp), optional, dimension(:,:), allocatable, intent(out) :: residual
    real(dp) :: omega_, tol_
    real(dp), dimension(size(var_idx)) :: res_tmp, res_init
    integer :: n, n_var, n_iter_
    logical :: converged_

    n_var = size(var_idx)
    omega_ = 1.3_dp
    tol_   = 1.0e-8_dp
    n_iter_ = 100
    if ( present(omega)     ) omega_     = omega
    if ( present(tol)       ) tol_       = tol
    if ( present(n_iter)    ) n_iter_    = n_iter
    if ( present(converged) ) converged  = .false.
    converged_ = .false.
    if (present(residual)) then
      if (allocated(residual)) deallocate(residual)
      allocate(residual(n_var,n_iter_))
    end if
    res_init = this%SOR_residual_norm(var_idx)

    if (any(res_init < epsilon(one))) then
      converged_ = .true.
      if ( present(converged) ) converged = converged_
      return
    end if

    do n = 1,n_iter_
      call this%SOR_iteration(var_idx,omega_,res_tmp)
      res_tmp = res_tmp / res_init
      write(*,*) n, res_tmp
      converged_ = all( res_tmp < tol)
      if ( present(residual ) ) residual(:,n) = res_tmp
      if ( present(converged) ) converged     = converged_
      if ( converged_ ) return
    end do
  end subroutine perform_iterative_reconstruction_SOR

  subroutine SOR_iteration(this,var_idx,omega,residual)
    use set_constants, only : zero, one
    class(reconstruction_holder), intent(inout) :: this
    integer,  dimension(:),       intent(in)    :: var_idx
    real(dp),                     intent(in)    :: omega
    ! real(dp), dimension(this%term_end-this%term_start,size(var_idx)), intent(out) :: residual
    real(dp), dimension(size(var_idx)), intent(out) :: residual
    real(dp), dimension(this%term_end-this%term_start,size(var_idx)) :: update
    integer :: i, n, j, vv, v, n_var
    
    residual = zero
    n_var = size(var_idx)
    do i = 1,this%n_cells
      update = zero
      do n = 1,this%rec(i)%n_interior
        j = this%rec(i)%nbor_idx(n)
        do vv = 1,n_var
          v = var_idx(vv)
          update(:,vv) = update(:,vv) + this%rec(i)%get_LHS_contribution_update(this%rec(j),n,this%term_start,this%term_end,var_idx(vv))
          update(:,vv) = update(:,vv) + this%rec(i)%get_RHS_contribution_update(this%rec(j),n,this%term_start,this%term_end,var_idx(vv))
        end do
      end do
      do vv = 1,n_var
        v = var_idx(vv)
        this%rec(i)%coefs(this%term_start+1:this%term_end,v) = (one-omega)*this%rec(i)%coefs(this%term_start+1:this%term_end,v) + omega*update(:,vv)
        residual = residual + sum( ( this%SOR_cell_residual(i,var_idx) )**2,dim=1 )
      end do
    end do
    residual = sqrt( residual )
  end subroutine SOR_iteration

  pure function SOR_residual_norm(this,var_idx) result(residual)
    use set_constants, only : zero
    class(reconstruction_holder), intent(in)    :: this
    integer,  dimension(:),       intent(in)    :: var_idx
    real(dp), dimension(size(var_idx))          :: residual
    integer :: i, vv, v, n_var

    residual = zero
    n_var = size(var_idx)
    do i = 1,this%n_cells
      do vv = 1,n_var
        v = var_idx(vv)
        residual = residual + sum( ( this%SOR_cell_residual(i,var_idx) )**2,dim=1 )
      end do
    end do
    residual = sqrt( residual )
  end function SOR_residual_norm

  pure function SOR_cell_residual(this,cell_idx,var_idx) result(residual)
    use set_constants, only : zero, one
    class(reconstruction_holder), intent(in)    :: this
    integer,                      intent(in)    :: cell_idx
    integer,  dimension(:),       intent(in)    :: var_idx
    real(dp), dimension(this%term_end-this%term_start,size(var_idx)) :: residual
    integer :: n, j, vv, v, n_var

    residual = zero

    n_var = size(var_idx)
    do vv = 1,n_var
      v = var_idx(vv)
      residual(:,vv) = residual(:,vv) + matmul( this%rec(cell_idx)%A, this%rec(cell_idx)%coefs(this%term_start+1:this%term_end,v) )
      residual(:,vv) = residual(:,vv) - matmul( this%rec(cell_idx)%D, this%rec(cell_idx)%coefs(1:this%term_start,v) )
      do n = 1,this%rec(cell_idx)%n_interior
        j = this%rec(cell_idx)%nbor_idx(n)
        residual(:,vv) = residual(:,vv) - matmul( this%rec(cell_idx)%B(:,:,n), this%rec(j)%coefs(this%term_start+1:this%term_end,v) )
        residual(:,vv) = residual(:,vv) - matmul( this%rec(cell_idx)%C(:,:,n), this%rec(j)%coefs(1:this%term_start,v) )
      end do
    end do
  end function SOR_cell_residual

end module var_rec_holder_type

module test_problem
  use set_precision, only : dp
  implicit none
  private
  public :: setup_grid_and_reconstruction
contains

  pure function test_function_1(n_var,x) result(val)
    integer,                intent(in) :: n_var
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(n_var)         :: val

    val(1) = 999.0_dp * x(1) - 888.0_dp * x(2) + 777.0_dp * x(3) - 666.0_dp
  end function test_function_1

  pure function test_function_2(n_var,x) result(val)
  use set_constants, only : pi
    integer,                intent(in) :: n_var
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(n_var)         :: val
    val(1) = sin(pi*x(1)) * sin(pi*x(2)) * sin(pi*x(3))
  end function test_function_2

  subroutine setup_grid_and_reconstruction( degree, n_vars, n_dim, n_nodes, n_ghost, grid, rec )
    use grid_derived_type, only : grid_type
    use linspace_helper,   only : unit_cartesian_mesh_cat
    use var_rec_holder_type, only : reconstruction_holder
    use lagrange_interpolation, only : generate_1D_barycentric_info, destroy_1D_barycentric_info
    integer,                     intent(in)  :: degree, n_vars, n_dim
    integer, dimension(3),       intent(in)  :: n_nodes, n_ghost
    type(grid_type),             intent(out) :: grid
    type(reconstruction_holder), intent(out) :: rec
    logical :: converged
    ! setup_reconstruction_holder( this, block_num, degree, n_vars, gblock )
    call generate_1D_barycentric_info()
    call grid%setup(1)
    call grid%gblock(1)%setup(n_dim,n_nodes,n_ghost)
    grid%gblock(1)%node_coords = unit_cartesian_mesh_cat(n_nodes(1),n_nodes(2),n_nodes(3))
    call grid%gblock(1)%grid_vars%setup( grid%gblock(1) )
    call rec%create(1,degree,n_vars,grid%gblock(1))
    call rec%setup_basic(grid%gblock(1),eval_fun=test_function_1)
    call rec%solve([1],omega=1.3_dp,tol=1e-10_dp,n_iter=100,converged=converged)
    write(*,*) 'converged =', converged
    call destroy_1D_barycentric_info()
  end subroutine setup_grid_and_reconstruction

end module test_problem

program main
  use set_precision, only : dp
  use set_constants, only : zero, one
  use test_problem,  only : setup_grid_and_reconstruction
  use grid_derived_type, only : grid_type
  use var_rec_holder_type, only : reconstruction_holder
  use timer_derived_type, only : basic_timer_t

  implicit none

  type(grid_type) :: grid
  type(reconstruction_holder) :: rec
  type(basic_timer_t) :: timer
  integer :: degree, n_vars, n_dim
  integer, dimension(3) :: n_nodes, n_ghost
  integer :: i, j, N, N_repeat
  real(dp) :: avg

  degree  = 3
  n_vars  = 1
  n_dim   = 2
  n_nodes = [5,5,2]
  n_ghost = [0,0,0]
  call setup_grid_and_reconstruction(degree, n_vars, n_dim, n_nodes, n_ghost, grid, rec )

  write(*,*) 'Here'
  call rec%destroy()
  call grid%destroy()
end program main


! program main
!   use set_precision, only : dp
!   use set_constants, only : zero, one
!   use quadrature_derived_type, only : quad_t, create_quad_ref_1D, create_quad_ref_2D
!   use monomial_basis_derived_type, only : monomial_basis_t
!   use var_rec_derived_type, only : var_rec_t
!   use timer_derived_type, only : basic_timer_t

!   implicit none

!   type(quad_t)    :: quad
!   type(monomial_basis_t) :: mono_basis
!   type(var_rec_t) :: rec
!   type(basic_timer_t) :: timer
!   integer :: self_block, self_idx, n_vars, degree, n_dim, n_interior
!   integer, dimension(:), allocatable :: nbor_idx, face_id
!   real(dp), dimension(:), allocatable :: h_ref

!   integer :: i, j, N, N_repeat
!   real(dp) :: avg

!   self_block = 1
!   self_idx   = 1
!   n_vars     = 4
!   degree     = 4
!   n_dim      = 3
!   n_interior = 2
!   allocate(nbor_idx(2)); nbor_idx = [2,3]
!   allocate(face_id(2)); face_id   = [2,3]
!   allocate(h_ref(n_dim)); h_ref = one

!   call create_quad_ref_2D(degree+1,quad)
!   mono_basis = monomial_basis_t(degree,n_dim)

!   N_repeat = 1000
!   avg = zero
!   do j = 1,N_repeat
!     call timer%tic()
!     rec =  var_rec_t( self_block, self_idx, n_interior, nbor_idx, face_id, n_vars, mono_basis, quad, h_ref )
!     avg = avg + timer%toc()
!   end do
!   write(*,*) 'Average Elapsed time: ', avg/real(N_repeat,dp)

!   deallocate( nbor_idx, face_id, h_ref )
!   call quad%destroy()
!   call rec%destroy()

! end program main