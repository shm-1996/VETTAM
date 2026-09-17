# VETTAM Runtime Parameters and Variables

Generated from the `Config` files in this release by `preparePublic.py`.
Do not edit by hand -- change the `Config` files instead.
Defaults shown are the Config defaults; override them in your `flash.par`.

## `physics/RadTrans/RadTransMain/VETTAM/Config`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `rt_compute_Dt` | Boolean | `FALSE` | compute a time step associated with cooling |
| `rt_debug` | Boolean | `FALSE` | debug flag |
| `rt_etens` | String | `"vet"` | Eddington tensor type |
| `rt_Dt_type` | Integer | `2` | choose prescription for the timestep calculation. 1 = light crossing time, 2 = Radn modified CFL |
| `rt_dtFactor` | Real | `1.0` | cooling timestep factor |
| `rt_T_min` | Real | `0.0` | cool only in regimes where temperature > T_rt_min |
| `rt_rtol` | Real | `1.e-8` | relative error tolerance |
| `rt_abstol` | Real | `1.e-20` | absolute error tolereance |
| `rt_dtol` | Real | `1.e5` | tolerance before diverged |
| `rt_maxits` | Integer | `1000` | maximum number of iterations |
| `rt_minits` | Integer | `1` | minimum number of iterations |
| `rt_preconditioner` | String | `"asm"` | preconditioner for linear solver, e.g bjacobi, jacobi, shell, none |
| `rt_linearsolver` | String | `"gmres"` | linear solver, e.g. gmres, cg, bcgs |
| `rt_smalle` | Real | `0.0` | lower limit for erad |
| `rt_smallf` | Real | `1.e-20` | lower limit for frad |
| `rt_hydro_type` | Integer | `2` | Temperature update type: 1=LTE, 2=Picard Iteration |
| `rt_temp_type` | Integer | `1` | Temperature update type: 1=Newton Rhapson, 2=Linearised, 3: Quartic |
| `rt_picard_rtol` | Real | `1.e-3` | Picard iteration rtol |
| `rt_picard_abstol` | Real | `1.e-10` | Picard iteration rtol |
| `rt_picard_maxits` | Integer | `1000` | Picard iteration maxsteps |
| `rt_picard_theta` | Real | `1.0` | Picard theta for next guess. 1 = newly computed guess, 2- previous guess |
| `rt_nr_rtol` | Real | `1.e-4` | NR rtol |
| `rt_nr_maxits` | Integer | `10000` | NR maxits |
| `rt_nr_dtol` | Real | `1.e5` | NR Dtol |
| `rt_hlle_correction` | Boolean | `TRUE` | HLLE Correction Flag |
| `rt_sink_implicit` | Boolean | `FALSE` | Implicit Sink Heating |
| `rt_picard_correction` | Boolean | `TRUE` | Picard temp correction |
| `rt_ovcterms` | Boolean | `TRUE` | O(v/c) correction terms flag |
| `noSink_VETTAM` | Boolean | `TRUE` | Flag whether to use VETTAM when there are no sink particles |
| `rt_dt_maxvel` | Real | `3.e10` | Parameters to control frequency of VETTAM (only use rt every dt = cfl*dx/rt_dt_maxvel; default is to resolve the light speed) |
| `xl_rad_BC` | Integer | `2` | boundary condition: 1=Reflecting, 2=Outflow, 3=Fixed, 4=Periodic, 5=Vacuum, 6=Marshak |
| `xr_rad_BC` | Integer | `2` | -- |
| `yl_rad_BC` | Integer | `2` | -- |
| `yr_rad_BC` | Integer | `2` | -- |
| `zl_rad_BC` | Integer | `2` | -- |
| `zr_rad_BC` | Integer | `2` | -- |
| `rt_update_hydro` | Boolean | `TRUE` | Update hydro quantites after radiation update |
| `rt_boundary_T` | Real | `10.` | boundary condition, fixed boundary radiation temperature |
| `fradx_boundary_value` | Real | `0.0` | boundary condition, fixed boundary radiation flux |
| `frady_boundary_value` | Real | `0.0` | -- |
| `fradz_boundary_value` | Real | `0.0` | -- |
| `sink_useRayTrace` | Boolean | `FALSE` | Use ray-trace for direct sink contribution; false by default, i.e. uses moment update |
| `sigma_star` | Real | `2.0` | sigma_star - smoothing length if Gaussian input for sinks used |
| `rt_sinkheat` | Boolean | `TRUE` | flag for sink heating contribution |
| `rt_sinkmom` | Boolean | `TRUE` | flag for sink momentum contribution |
| `rt_freqbands` | Integer | `3` | Multiband stuff -- by default 3 frequency bands |
| `rt_band1` | String | `"IR"` | Frequency Band information : by default IR, FUV and EUV respectively |
| `rt_band2` | String | `"FUV"` | -- |
| `rt_band3` | String | `"EUV"` | -- |
| `rt_band4` | String | `"PE"` | -- |
| `rt_band5` | String | `"LYMAN_WERNER"` | -- |

| Variable | Kind | Description |
| --- | --- | --- |
| `erad` | variable | Variables |
| `mohx` | variable | -- |
| `mohy` | variable | -- |
| `mohz` | variable | -- |
| `taup` | variable | Opacities: taup planck opacities, and taur Rosseland opacities |
| `taur` | variable | -- |
| `erao` | variable | Old time states of radiation variables |
| `mhxo` | variable | -- |
| `mhyo` | variable | -- |
| `mhzo` | variable | -- |
| `tmpo` | variable | -- |
| `erap` | variable | Picard updated versions of radiation variables |
| `mhxp` | variable | -- |
| `mhyp` | variable | -- |
| `mhzp` | variable | -- |
| `tmpp` | variable | -- |
| `xxed` | variable | Components of the Eddington Tensor (f_Eddington: 3x3 matrix) |
| `xyed` | variable | -- |
| `xzed` | variable | -- |
| `yyed` | variable | -- |
| `yzed` | variable | -- |
| `zzed` | variable | -- |
| `jstr` | variable | Source term for sink radiation |
| `reir` | variable | Reprocessed radiation to IR from higher energy radiation |

## `physics/RadTrans/RadTransMain/VETTAM/LWBand/Config`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `useH2Dissociate` | Boolean | `TRUE` | flag to use ionising radiation module |
| `lwdiss_type` | String | `"raytrace"` | H2 dissociation type |
| `bfive` | Real | `1.0` | Doppler b parameter used in fshield equation in units of km/s (see Draine and Bertoldi 1996); default value of 1 |
| `fpump` | Real | `6.94` | The parameter fpump that quantifies the probability of UV pumping vs dissociation; Eq 45 of Baczynski et al 2015 |
| `energyPerDissociation` | Real | `6.408705848000027e-13` | The thermal energy deposited per photodissociation event; using 0.4 eV by default (Black & Dalgarno 1977) |
| `avgEnergyLW` | Real | `1.954655283640008e-11` | The average energy of a LW photon -- this would depend on the spectrum, but we use 12.2 eV by default -- typical of a stellar population (Kim+ 23) |
| `useHIshield` | Boolean | `TRUE` | Flag to switch on/off H2 shielding due to H |

| Variable | Kind | Description |
| --- | --- | --- |
| `fshm` | variable | Variable Storing the Shielding Factor of mol H2 for the cell |
| `fsha` | variable | Variable Storing the Shielding Factor of atomic H2 for the cell |
| `lwfl` | variable | Variable storing the self-shielded LW dissociation flux |

## `physics/RadTrans/RadTransMain/VETTAM/LWBand/MetalKrome/Config`

| Variable | Kind | Description |
| --- | --- | --- |
| `ulwd` | variable | Variable storing the absorbed energy in LW band by dust (for the PE heating rate, and for subsequent self-shielded H2, CO and C dissociation/ionization) |
| `uped` | variable | Variable storing the absorbed energy in the PE band by dust (for the PE heating rate) |
| `dih2` | variable | Variable storing the shielded dissociation rate of H2 |
| `ic` | variable | Variable storing the shielded ionization rate of Carbon |
| `dico` | variable | Variable storing the shielded dissociation rate of CO |

## `physics/RadTrans/RadTransMain/VETTAM/PhotoIonization/Config`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `useEUVIonize` | Boolean | `TRUE` | flag to use ionising radiation module |
| `useH2Ionize` | Boolean | `TRUE` | flag to ionise molecular hydrogen |
| `ion_OTS` | Boolean | `TRUE` | flag to use the on-the-spot approximation |
| `ion_sigmaH` | Real | `3.1e-18` | Cross section of H to LyC photons with >13.6eV averaged over SB99 luminosities (Kim+23) |
| `ion_sigmaH_13p6_15p2` | Real | `5.6e-18` | Cross section of H to LyC photons with energies 13.6<hnu<15.2 (Baczynski+ 2015) |
| `ion_sigmaH_15p2_infty` | Real | `3.1e-18` | Cross section of H to LyC photons with energies hnu>15.2 (Baczynski+ 2015) |
| `ion_sigmaH2_15p2_infty` | Real | `5.5e-18` | Cross section of H2 to LyC photons with energies > 15.2eV (Baczynski+ 2015) |
| `ion_sigmaH2` | Real | `4.6e-18` | Cross section of H2 to LyC photons with energies > 13.6eV averaged over SB99 luminosities (Kim+ 23) |
| `alpha_type` | String | `"default"` | Recombination coefficient type: either "default" for a temp-dependent fit or "constant"; latter useful for standard tests |
| `alpha_B` | Real | `2.6e-13` | Recombination rate coefficients: Case B and Case A coefficient values - used only if alpha_type is "constant" |
| `alpha_A` | Real | `4.1e-13` | -- |
| `ion_type` | Integer | `1` | -- |
| `ion_implicit` | Boolean | `TRUE` | -- |
| `C_tss` | Real | `1.0` | Substepping timestep factor for ionisation-recombination balance (for methods 1 and 2) |
| `C_therm` | Real | `100.0` | Substepping timestep factor for heating-cooling balance (if Ideal Gas EOS used) |
| `ion_conv_rtol` | Real | `1.e-2` | Tolerances of ionisation-radiation coupling (methods 3 and 4) |
| `ion_conv_atol` | Real | `1.e-4` | -- |
| `ion_conv_maxits` | Integer | `100` | -- |
| `energyPerIonH` | Real | `5.447399970800022e-12` | Energy per Ionisation of H in ergs: Default is 3.4 eV, app for Star Cluster (Kim+23, Table 3) |
| `energyPerIonH2` | Real | `7.04957643280003e-12` | Energy per Ionisation of H2 in ergs: Default is 4.4 eV, app for Star Cluster (Kim+23, Table 3) |
| `hnu` | Real | `2.883917631600012e-11` | Average energy per ionising photons (18 eV by default. Will affect dust heating by EUV) |

| Variable | Kind | Description |
| --- | --- | --- |
| `ionh` | variable | Cumulative heat added per unit volume via ionisations for this cell |
| `ionp` | variable | -- |
| `iono` | variable | -- |
| `ueuv` | variable | -- |
| `tauh` | variable | -- |
| `heuv` | variable | -- |

## `physics/RadTrans/RadTransMain/VETTAM/PhotoIonization/OnlyH/Config`

| Variable | Kind | Description |
| --- | --- | --- |
| `iony` | mass scalar | Neutral Fraction Mass Scalar |

## `physics/RadTrans/RadTransMain/VETTAM/PhotoIonization/ThermoChemistry/Config`

| Variable | Kind | Description |
| --- | --- | --- |
| `ih` | variable | -- |
| `ih2` | variable | -- |

## `physics/RadTrans/RadTransMain/VETTAM/RayTrace/3DRT/Config`

| Variable | Kind | Description |
| --- | --- | --- |
| `mean` | variable | mean intensity |
| `rpxx` | variable | second moment of radiation intensity (radiation pressure tensor) |
| `rpxy` | variable | -- |
| `rpxz` | variable | -- |
| `rpyy` | variable | -- |
| `rpyz` | variable | -- |
| `rpzz` | variable | -- |

## `physics/RadTrans/RadTransMain/VETTAM/RayTrace/3DRT/healpix/Config`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `rt_healpix_nSide` | Integer | `1` | -- |

## `physics/RadTrans/RadTransMain/VETTAM/RayTrace/Config`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `useRayTrace` | Boolean | `TRUE` | useRayTrace flag to indicate whether to use the raytracer |
| `rt_ALI` | Integer | `1` | lambda iteration 0 = simple lambda iteration 1 = ALI with diagonal lambda operator 2 = ALI with tri-diagonal lambda operator |
| `rt_nPhi` | Integer | `1` | -- |
| `rt_nTheta` | Integer | `1` | -- |
| `rt_nrOfAngleGroups` | Integer | `1` | the number of angle groups for communication |
| `rt_healpix_randomize` | Integer | `1` | rt_healpix_randomize Randomize the healpix sphere orientation every nth call of the raytracer. 0 = No randomization |
| `rt_maxNrOfBoundIter` | Integer | `0` | max number of periodic boundary conditions iterations |
| `rt_dirX` | Real | `1.0` | -- |
| `rt_dirY` | Real | `0.0` | -- |
| `rt_dirZ` | Real | `0.0` | -- |
| `rt_irradiation` | Real | `0.0` | -- |
| `rt_epsilon` | Real | `1.0` | rt_epsilon This is the thermal coupling parameter for the source function S = (1-epsilon) * J + epsilon * B |

| Variable | Kind | Description |
| --- | --- | --- |
| `sour` | variable | source function |

## `physics/materialProperties/Opacity/RadTrans/DustOpacities/Config`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `opac_floor` | Real | `1.e-5` | -- |
| `dust_opacity_ir` | String | `"Semenov"` | Flag for opacities. Values possible: "Semenov", "Fixed". |
| `dust_opacity_fuv` | String | `"Fixed"` | -- |
| `dust_opacity_euv` | String | `"Fixed"` | -- |
| `dust_euvopac_value` | Real | `1.25e3` | Opacities in the bands if Fixed option is chosen. EUV opacity assumes pure hydrogen mixture; divide by eos_singlespeciesA These values are motivated by star cluster SED-weighted averages shown in Figure 22, Kim et al. 2023 Bands: FUV, EUV, Photoelectric (PE), Lyman-Werner (LW), Infrared (IR) |
| `dust_fuvopac_value` | Real | `600.0` | -- |
| `dust_peopac_value` | Real | `514.0` | -- |
| `dust_lwopac_value` | Real | `900.0` | -- |
| `dust_iropac_value` | Real | `10.0` | -- |
| `dusttoGasRatio` | Real | `1.0` | Dust to gas ratio relative to the solar value |

| Variable | Kind | Description |
| --- | --- | --- |
| `opac` | variable | total extinction coefficient |
