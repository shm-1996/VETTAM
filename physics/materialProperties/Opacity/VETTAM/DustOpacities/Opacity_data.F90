module Opacity_data
    implicit none
#include "constants.h"
    !dust opacity types in each band
    character(len=MAX_STRING_LENGTH), save :: dust_opacity_ir, dust_opacity_fuv, dust_opacity_euv
    !Opacity in different bands if opacity type is constant
    real, save :: dusttoGasRatio, dust_iropac_value,dust_fuvopac_value, dust_euvopac_value, dust_lwopac_value, dust_peopac_value
end module Opacity_data