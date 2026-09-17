# VETTAM

**V**ariable **E**ddington **T**ensor-closed **T**ransport on **A**daptive **M**eshes -- a
multi-frequency radiation-hydrodynamics module coupled to the the FLASH MHD code. $\texttt{VETTAM}$ solves the two-moment radiation
hydrodynamical equations closed with a non-local Variable Eddington Tensor (VET) obtained with a ray-trace solution. The implementation also offers the option to use the Flux Limited Diffusion (FLD) and Moment-1 ($M_1$) closures instead of the VET. The temporal update for VETTAM utilizes an implicit approach, solved using sparse matrix solvers available with the open-source $\texttt{PETSc}$ library. 

If you use VETTAM, please cite **Menon et al. 2022, MNRAS 512, 401**.

## Repository structure

VETTAM is a FLASH unit, not a standalone program. This repository holds the unit
directories only, laid out exactly as they sit under a FLASH `source/` tree:

```
physics/RadTrans/RadTransMain/VETTAM/          #the module
physics/materialProperties/Opacity/RadTrans/   #dust opacities
```
The subdirectories include additional functionalities which you should include if your problem requires them:
```
VETTAM/RayTrace          #the ray-tracer for the VET
VETTAM/PhotoIonization   #Photoionization (of H & H2 for now)
VETTAM/LWBand         #Lyman-werner band dissociation (H2,CO) and ionization (CI)
```

FLASH itself is not included. You supply your own copy, obtained from the Flash Center
for Computational Science (https://flash.rochester.edu/site/flashcode/).

## Requirements

- FLASH 4.0.1 or a derivative, with $\texttt{PARAMESH}$ AMR functionality (VETTAM reads the `tree`)
- $\texttt{PETSc}$, built with MPI and Fortran capability. 
- An MPI Fortran compiler

## Installing

From the root of this repository:

```bash
./install.sh /path/to/your/flash/source
```

That copies both unit directories into place. Pass `--link` to symlink instead, which is
convenient while developing against a checkout. Nothing outside these two directories is
touched, and no FLASH core file is modified.

Then point your `Makefile.h` file in the machine-specific directory within `sites/` to the install path of $\texttt{PETSc}$ :

```make
PETSC_PATH   = /path/to/petsc
FFLAGS_PETSC = -I$(PETSC_PATH)/include
LIB_PETSC    = -L$(PETSC_PATH)/lib -lpetsc
```

## Using it

For the barebone module, include the following unit from your simulation `Config`:

```
REQUIRES physics/RadTrans/RadTransMain/VETTAM
REQUIRES physics/RadTrans/RadTransMain/VETTAM/RayTrace # If using VET
REQUIRES physics/RadTrans/RadTransMain/VETTAM/Photoionization # If including photoionization
REQUIRES physics/RadTrans/RadTransMain/VETTAM/LWBand # If including LW dissociation/ionization
```

### Frequency bands

`rt_freqbands` (default 3; maximum 5 for now) bands are processed sequentially each hydro step, named by
`rt_band1`..`rt_band5` (defaults: IR, FUV, EUV, PE, LYMAN_WERNER; order them from low to
high energy). The radiation energy density (`ERAD`) and fluxes (`MOHX/Y/Z`) are reused by each band in turn, so to retain a
band's radiation moments in the output files define `BnER`, `BnFX`, `BnFY`, `BnFZ` (n = band number) in your
simulation `Config`. Using the $M_1$ closure with more than one band *requires* these and aborts at startup without them.

### Sink particle sources

To include the contribution of radiation field from point sources (such as sink particles), build with the sink unit *and* any
stellar-evolution sub-unit that sets the luminosities in different bands. The implementation expects the variables `LUMINOSITY_PART_PROP` and
`STELLAR_RADIUS_PART_PROP`to be set by the sink routines to maintain the modularity of $\texttt{VETTAM}$. 

### Optional submodules

Each is opt-in via `REQUIRES` and gated at compile time. Their prerequisites are not all
declared in their own `Config` files, so they are spelled out here:

| Submodule | Adds | Also requires |
| --- | --- | --- |
| `RayTrace/` | 3D ray-traced Eddington tensors (`RAYTRACE_3DRT`) | Nothing beyond the core. Needed for `rt_etens="vet"`, the default closure |
| `PhotoIonization/OnlyH` | EUV photoionization via the `iony` mass scalar | Nothing beyond the core |
| `PhotoIonization/ThermoChemistry` | EUV photoionization writing rates for a chemistry network | Multispecies `H`, `HP`, `ELEC` (and `H2` for H2 ionization). This is the `PhotoIonization` Config default, and it `REQUESTS` KROME |
| `LWBand/Manual` | Lyman-Werner H2 dissociation, computed in-module | `RayTrace/`, `PhotoIonization/` with Multispecies `H` + `H2`, and the sink unit with `StellarEvolution` |
| `LWBand/PrimKrome`, `LWBand/MetalKrome` | LW rates handed to a KROME network | As `LWBand/Manual`, plus KROME |

KROME is not distributed with this release. Without it, use `PhotoIonization/OnlyH`,
which needs nothing beyond the core:

```
REQUIRES physics/RadTrans/RadTransMain/VETTAM
REQUIRES physics/RadTrans/RadTransMain/VETTAM/RayTrace
REQUIRES physics/RadTrans/RadTransMain/VETTAM/PhotoIonization/OnlyH
```

Every `LWBand` flavour needs Multispecies `H_SPEC` and `H2_SPEC` at compile time. The
`#ifndef H2_SPEC` check inside `rt_lw.F90` is a runtime abort and will not spare you a
compile error, so declare the species in your simulation `Config`, or pull in a
chemistry unit that does.

## Documentation

- `docs/VETTAM_Module_Documentation.md` -- architecture and feature guide
- `docs/RUNTIME_PARAMETERS.md` -- every runtime parameter and solution variable,
  generated from the `Config` files in this release

## Contact

Shyam Harimohan Menon -- smenon@flatironinstitute.org
