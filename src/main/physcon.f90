
module physcon 

 implicit none

!--Mathematical constants
 real(kind=8), parameter :: pi       =  3.1415926536d0

!--Physical constants
 real(kind=8), parameter :: c = 2.997924d10                     !Speed of light            cm/s
 real(kind=8), parameter :: gg = 6.672041d-8                    !Gravitational constant    dyn cm^2 g^-2 or cm^3 s^-2 g^-1
 real(kind=8), parameter :: mass_electron_cgs = 9.10938291d-28  !Electron mass             g
 real(kind=8), parameter :: mass_proton_cgs = 1.67262158d-24    !Proton mass               g
 real(kind=8), parameter :: kboltz = 1.38066d-16                !Boltzmann constant        erg/K
 real(kind=8), parameter :: steboltz   = 5.67051d-5             !Stefan-Boltzmann constant erg cm^-2K^-4 s^-1

!--Mass scale
 real(kind=8), parameter :: solarm = 1.9891d33                  !Mass of the Sun           g
 real(kind=8), parameter :: solarr = 6.959500d10                !Radius of the Sun         cm
 real(kind=8), parameter :: solarl = 3.9d33                     !Luminosity of the Sun     erg/s
 real(kind=8), parameter :: gram = 1.d0

!--Distance scale
 real(kind=8), parameter :: au = 1.496d13                       !Astronomical unit         cm
 real(kind=8), parameter :: ly = 9.4605d17                      !Light year                cm
 real(kind=8), parameter :: pc = 3.086d18                       !Parsec                    cm
 real(kind=8), parameter :: kpc = 3.086d21                      !Kiloparsec                cm
 real(kind=8), parameter :: Mpc = 3.086d24                      !Megaparsec                cm
 real(kind=8), parameter :: cm = 1.d0                           !Centimetre                cm

!--Time scale
 real(kind=8), parameter :: years = 3.1556926d7


end module physcon 