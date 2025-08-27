module set_precision
  use iso_c_binding, only : c_double
  implicit none
  private
  public :: r8, dp
  integer, parameter :: r8  = c_double
  integer, parameter :: dp  = r8
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

end module combinatorics

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

module quadrature_derived_type

  use set_precision,       only : dp
  use set_constants,       only : zero
  implicit none
  private
  public :: quad_t
  public :: create_quad_ref_1D, create_quad_ref_2D
  type quad_t
    integer :: n_quad = 0
    real(dp), allocatable, dimension(:,:) :: quad_pts
    real(dp), allocatable, dimension(:)   :: quad_wts
  contains
    private
    procedure, public :: allocate_quad
    procedure, public :: deallocate_quad
    generic,   public :: integrate => integrate_scalar, integrate_vector
    procedure :: integrate_scalar
    procedure :: integrate_vector
  end type quad_t

contains

  subroutine allocate_quad( this, n_quad )
    use set_constants, only : zero
    class(quad_t), intent(inout) :: this
    integer,       intent(in)    :: n_quad
    this%n_quad = n_quad
    allocate( this%quad_pts(3,n_quad) )
    this%quad_pts = zero
    allocate( this%quad_wts(n_quad) )
    this%quad_wts = zero
  end subroutine allocate_quad

  subroutine deallocate_quad( this )
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

  elemental subroutine LegendrePolynomialAndDerivative(N,x,LN,dLN)
    use set_constants, only : one, two
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
    use set_constants, only : one, two, four, third, pi
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

  function gauss_1D_size( polynomial_order ) result( N_quad )
    use set_constants, only : half
    integer, intent(in) :: polynomial_order
    integer             :: N_quad
    N_quad = ceiling( half*(polynomial_order + 1) )
  end function gauss_1D_size

  subroutine create_quad_ref_1D( quad_order, quad_ref )
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    real(dp), dimension(:), allocatable :: xtmp
    integer :: n_quad
    n_quad = gauss_1D_size( quad_order )
    allocate( xtmp(n_quad) )
    call quad_ref%deallocate_quad()
    call quad_ref%allocate_quad( n_quad )
    call LegendreGaussNodesAndWeights( n_quad, xtmp, quad_ref%quad_wts )
    quad_ref%quad_pts(1,:) = xtmp
    deallocate( xtmp )
  end subroutine create_quad_ref_1D

  subroutine create_quad_ref_2D( quad_order, quad_ref )
    use set_constants,           only : zero
    integer,      intent(in)  :: quad_order
    type(quad_t), intent(out) :: quad_ref
    integer :: n_quad
    integer :: i, j, cnt
    real(dp), dimension(:), allocatable :: pts_1D
    real(dp), dimension(:), allocatable :: wts_1D
    continue
    n_quad = gauss_1D_size( quad_order )
    allocate( pts_1D(n_quad) )
    allocate( wts_1D(n_quad) )
    call LegendreGaussNodesAndWeights(n_quad-1, pts_1D, wts_1D)
    call quad_ref%deallocate_quad()
    call quad_ref%allocate_quad( n_quad**2 )
    cnt = 0
    do j = 1, n_quad
      do i = 1, n_quad
        cnt = cnt + 1
        quad_ref%quad_pts(:,cnt) = [ pts_1D(i), pts_1D(j), zero ]
        quad_ref%quad_wts(cnt) = wts_1D(i)*wts_1D(j)
      end do
    end do
    deallocate( wts_1D )
    deallocate( pts_1D )
  end subroutine create_quad_ref_2D
end module quadrature_derived_type

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
    integer :: n

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
  class(zero_mean_basis_t),        intent(in) :: this
  integer,                         intent(in) :: term_start, term_end
  real(dp), dimension(this%n_dim), intent(in) :: point
  real(dp), dimension(this%n_dim), intent(in) :: scale
  real(dp), dimension(term_end, term_end - term_start) :: derivatives
  integer :: i, j
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
    private
    integer :: n_vars
    integer :: n_interior
    integer :: n_boundary
    integer :: self_idx
    integer :: self_block
    integer, dimension(:), allocatable      :: int_idx, bnd_idx
    type(zero_mean_basis_t)                 :: basis
    ! type(quad_t), pointer                   :: quad
    real(dp), dimension(:,:),   allocatable :: coefs
    real(dp), dimension(:,:),   allocatable :: A_inv
    real(dp), dimension(:,:),   allocatable :: A_inv_D
    real(dp), dimension(:,:,:), allocatable :: A_inv_B
    real(dp), dimension(:,:,:), allocatable :: A_inv_C
  contains
    private
    procedure, public, pass :: setup => setup_reconstruction
    ! procedure, public, pass :: eval  => eval_reconstruction
    procedure, pass :: setup_reconstruction_interior_faces
    procedure, pass :: get_LHS_interior
    procedure, pass :: get_RHS_interior
    procedure, pass :: get_RHS_interior_update
    procedure, pass :: get_LHS_dirichlet
    procedure, pass :: get_RHS_dirichlet
    procedure, public, pass :: destroy => destroy_var_rec_t

  end type var_rec_t

  interface var_rec_t
    procedure constructor
  end interface var_rec_t

contains

function constructor( self_block, self_idx, int_idx, bnd_idx, n_vars, mono_basis, quad, h_ref ) result(this)
  integer,                    intent(in) :: self_block, self_idx
  integer, dimension(:),      intent(in) :: int_idx, bnd_idx
  integer,                    intent(in) :: n_vars
  type(monomial_basis_t),     intent(in) :: mono_basis
  type(quad_t),               intent(in) :: quad
  real(dp), dimension(mono_basis%n_dim), intent(in) :: h_ref
  type(var_rec_t)                        :: this

  call this%destroy()

  this%self_block = self_block
  this%self_idx   = self_idx
  this%n_interior = size(int_idx)
  allocate( this%int_idx( this%n_interior ) )
  this%int_idx    = int_idx

  this%n_boundary = size(bnd_idx)
  allocate( this%bnd_idx( this%n_boundary ) )
  this%bnd_idx    = bnd_idx
  this%n_vars = n_vars
  this%basis = zero_mean_basis_t(mono_basis,quad,h_ref)
end function constructor

pure subroutine destroy_var_rec_t(this)
  class(var_rec_t), intent(inout) :: this
  call this%basis%destroy()
  if ( allocated(this%int_idx) ) deallocate( this%int_idx )
  if ( allocated(this%bnd_idx) ) deallocate( this%bnd_idx )
  if ( allocated(this%coefs  ) ) deallocate( this%coefs   )
  if ( allocated(this%A_inv  ) ) deallocate( this%A_inv   )
  if ( allocated(this%A_inv_B) ) deallocate( this%A_inv_B )
  if ( allocated(this%A_inv_C) ) deallocate( this%A_inv_C )
  if ( allocated(this%A_inv_D) ) deallocate( this%A_inv_D )
end subroutine destroy_var_rec_t

subroutine setup_reconstruction_interior_faces(this,nbors,fquads,term_start,term_end)
  class(var_rec_t),                             intent(inout) :: this
  class(var_rec_t), dimension(this%n_interior), intent(in)    :: nbors
  class(quad_t),    dimension(this%n_interior), intent(in)    :: fquads
  integer,                                      intent(in)    :: term_start, term_end
  real(dp), dimension(term_end - term_start, term_end - term_start)                  :: A
  real(dp), dimension(term_end - term_start, term_end - term_start, this%n_interior) :: B
  real(dp), dimension(term_end - term_start, term_start)                  :: D
  real(dp), dimension(term_end - term_start, term_start, this%n_interior) :: C

  integer :: j

  call this%get_LHS_interior(nbors,fquads,term_start,term_end,A,B)
  call this%get_RHS_interior(nbors,fquads,term_start,term_end,D,C)

  if ( allocated(this%A_inv  ) ) deallocate( this%A_inv   )
  if ( allocated(this%A_inv_B) ) deallocate( this%A_inv_B )
  if ( allocated(this%A_inv_C) ) deallocate( this%A_inv_C )
  if ( allocated(this%A_inv_D) ) deallocate( this%A_inv_D )

  allocate( this%A_inv(   term_end-term_start,term_end-term_start ) )
  allocate( this%A_inv_B( term_end-term_start,term_end-term_start, this%n_interior ) )
  allocate( this%A_inv_C( term_end-term_start, term_start, this%n_interior ) )
  allocate( this%A_inv_D( term_end-term_start, term_start ) )

  ! real(dp), dimension( term_end - term_start, term_start ), intent(out) :: D
  ! real(dp), dimension( term_end - term_start, term_start, this%n_interior), intent(out) :: C

end subroutine setup_reconstruction_interior_faces

subroutine setup_reconstruction( this, nbors, fquads, term_start, term_end )
  class(var_rec_t),                             intent(inout) :: this
  class(var_rec_t), dimension(this%n_interior), intent(in)    :: nbors
  class(quad_t),    dimension(this%n_interior), intent(in)    :: fquads
  integer,                                      intent(in)    :: term_start, term_end

  call this%setup_reconstruction_interior_faces(nbors,fquads,term_start,term_end)

  ! TODO: boundary face contributions
  ! this%A_inv = inv(A)
  ! do j = 1,this%n_interior
  !   this%A_inv_B(:,:,j) = matmul( this%A_inv, B(:,:,j) )
end subroutine setup_reconstruction

pure subroutine get_LHS_interior( this, nbors, fquads, term_start, term_end, A, B )
  use set_constants, only : zero, one
  class(var_rec_t),                             intent(in) :: this
  class(var_rec_t), dimension(this%n_interior), intent(in) :: nbors
  class(quad_t),    dimension(this%n_interior), intent(in) :: fquads
  integer,                                      intent(in) :: term_start, term_end
  real(dp), dimension(term_end - term_start, term_end - term_start), intent(out) :: A
  real(dp), dimension(term_end - term_start, term_end - term_start, this%n_interior), intent(out) :: B
  real(dp), dimension(term_end,term_end-term_start) :: d_basis_i, d_basis_j
  integer :: j, q, l, m
  real(dp), dimension(this%basis%n_dim) :: dij
  real(dp) :: xdij_mag

  A = zero
  B = zero
  do j = 1,this%n_interior
    dij = abs( this%basis%x_ref - nbors(j)%basis%x_ref )
    xdij_mag = one/norm2(dij)
    do q = 1,fquads(j)%n_quad
      d_basis_i = this%basis%scaled_basis_derivatives( term_start, term_end, fquads(j)%quad_pts(:,q), dij )
      d_basis_j = nbors(j)%basis%scaled_basis_derivatives( term_start, term_end, fquads(j)%quad_pts(:,q), dij )
      do m = 1,term_end-term_start
        do l = 1,term_end-term_start
          A(l,m)   = A(l,m)   + fquads(j)%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_i(:,m) )
          B(l,m,j) = B(l,m,j) + fquads(j)%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_j(:,m) )
        end do
      end do
    end do
  end do

end subroutine get_LHS_interior

pure subroutine get_RHS_interior( this, nbors, fquads, term_start, term_end, D, C )
  use set_constants, only : zero, one
  class(var_rec_t),                             intent(in) :: this
  class(var_rec_t), dimension(this%n_interior), intent(in) :: nbors
  class(quad_t),    dimension(this%n_interior), intent(in) :: fquads
  integer,                                      intent(in) :: term_start, term_end
  real(dp), dimension( term_end - term_start, term_start ), intent(out) :: D
  real(dp), dimension( term_end - term_start, term_start, this%n_interior), intent(out) :: C
  real(dp), dimension( term_end, term_end - term_start ) :: d_basis_i
  real(dp), dimension( term_end, term_start )            :: d_basis_j
  integer :: j, q, l, m
  real(dp), dimension(this%basis%n_dim) :: dij
  real(dp) :: xdij_mag

  C = zero
  D = zero
  do j = 1,this%n_interior
    dij = abs( this%basis%x_ref - nbors(j)%basis%x_ref )
    xdij_mag = one/norm2(dij)
    do q = 1,fquads(j)%n_quad
      d_basis_i =     this%basis%scaled_basis_derivatives( term_start, term_end, fquads(j)%quad_pts(:,q), dij )
      d_basis_j = nbors(j)%basis%scaled_basis_derivatives( 0, term_start, fquads(j)%quad_pts(:,q), dij )
      do m = 1,term_start
        do l = 1,term_end-term_start
          D(l,m)   = D(l,m)   + fquads(j)%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_i(:,m) )
          C(l,m,j) = C(l,m,j) + fquads(j)%quad_wts(q) * xdij_mag * dot_product( d_basis_i(:,l), d_basis_j(:,m) )
        end do
      end do
    end do
  end do

end subroutine get_RHS_interior

pure subroutine get_RHS_interior_update( this, nbors, term_start, term_end, var_idx, b )
  use set_constants, only : zero, one
  class(var_rec_t),                             intent(in) :: this
  class(var_rec_t), dimension(this%n_interior), intent(in) :: nbors
  integer,                                      intent(in) :: term_start, term_end
  integer,          dimension(:),               intent(in) :: var_idx
  real(dp), dimension( term_end - term_start,size(var_idx) ), intent(out) :: b
  integer :: v, j

  b = zero
  do v = 1,size(var_idx)
    do j = 1,this%n_interior
      b(:,v) = b(:,v) + matmul( this%A_inv_C(:,:,j), nbors(j)%coefs(term_start+1:term_end,var_idx(v)) )
    end do
    b(:,v) = b(:,v) + matmul( this%A_inv_D, this%coefs(term_start+1:term_end,var_idx(v)) )
  end do
end subroutine get_RHS_interior_update

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

program main
  use set_precision, only : dp
  use set_constants, only : zero, one
  use quadrature_derived_type, only : quad_t, create_quad_ref_1D, create_quad_ref_2D
  use monomial_basis_derived_type, only : monomial_basis_t
  use var_rec_derived_type, only : var_rec_t
  use timer_derived_type, only : basic_timer_t

  type(quad_t)    :: quad
  type(monomial_basis_t) :: mono_basis
  type(var_rec_t) :: rec
  type(basic_timer_t) :: timer
  integer :: self_block, self_idx, n_vars, degree, n_dim
  integer, dimension(:), allocatable :: int_idx, bnd_idx
  real(dp), dimension(:), allocatable :: h_ref

  integer :: i, j, N, N_repeat
  real(dp) :: avg

  self_block = 1
  self_idx   = 1
  n_vars     = 1
  degree     = 10
  n_dim      = 3
  allocate(int_idx(2)); int_idx = [2,3]
  allocate(bnd_idx(0))
  allocate(h_ref(n_dim)); h_ref = one

  call create_quad_ref_2D(degree+1,quad)
  mono_basis = monomial_basis_t(degree,n_dim)

  N_repeat = 1000
  avg = zero
  do j = 1,N_repeat
    call timer%tic()
    rec =  var_rec_t( self_block, self_idx, int_idx, bnd_idx, n_vars, mono_basis, quad, h_ref )
    avg = avg + timer%toc()
  end do
  write(*,*) 'Average Elapsed time: ', avg/real(N_repeat,dp)

  deallocate( int_idx, bnd_idx, h_ref )
  call quad%deallocate_quad()
  call rec%destroy()

end program main