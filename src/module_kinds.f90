module module_kinds
  !================================================================!
  ! Kind parameters and constants used throughout the program.     !
  !================================================================!
  implicit none

  ! kinds of integer, real, complex and logical variables
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  integer, parameter :: lgt = kind(.true.)

  ! mathematical constants
  real(dp), parameter :: pi_d     = 3.141592653589793238462643383279502884_dp
  real(dp), parameter :: twopi_d  = 2.0_dp*pi_d
  real(dp), parameter :: fourpi_d = 4.0_dp*pi_d
  complex(dpc), parameter :: zero = (0.0_dp,0.0_dp)
  complex(dpc), parameter :: one  = (1.0_dp,0.0_dp)
  complex(dpc), parameter :: ii   = (0.0_dp,1.0_dp)

  ! gravitational constant in SI units (the value used in minos)
  real(dp), parameter :: bigg = 6.6723e-11_dp

end module module_kinds
