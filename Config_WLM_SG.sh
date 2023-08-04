#!/bin/bash            # this line only there to enable syntax highlighting in this file

####################################################################################################
#  Enable/Disable compile-time options as needed: this is where you determine how the code will act
#  From the list below, please activate/deactivate the
#       options that apply to your run. If you modify any of these options,
#       make sure that you recompile the whole code by typing "make clean; make".
#
#  Consult the User Guide before enabling any option. Some modules are proprietary -- access to them
#    must be granted separately by the code authors (just having the code does NOT grant permission).
#    Even public modules have citations which must be included if the module is used for published work,
#    these are all given in the User Guide.
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel (volker.springel@h-its.org). The code has been modified
#   substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO (to add new modules and clean
#   up the naming conventions and changed many of them to match the new GIZMO conventions)
#
####################################################################################################

GALSF
STELLAR_FEEDBACK
PHOTO_IONIZATION
STOCHASTIC_IMF
INSTANT_SN_FEEDBACK_FOR_SOME_TIME
OUTPUT_STARFORMATION_INFO
#JEANS_LENGTH_THRESHOLD
JEANS_MASS_THRESHOLD
SF_ONLY_IN_CONV_FLOWS
SF_INSTANT_CUTOFF
OUTPUT_VELGRAD

#MAGNETIC                       # master switch for MHD, regardless of which Hydro solver is used
#MHD_B_SET_IN_PARAMS            # set initial fields (Bx,By,Bz) in parameter file

#SAMPLE_IMF
#N_STELLAR_MASS=300
#INDIVIDUAL_STARS_SPLIT

WSS_CIE_COOL 
CHEMCOOL 
CHEMISTRYNETWORK=5 
#OUTPUT_INDIVIDUAL_COOLRATES
#OUTPUT_SHIELD_FAC

G0_VARIABLE

#G0_SCALE_WITH_TOTAL_SFR
#CR_SCALE_WITH_G0   #need G0_SCALE_WITH_TOTAL_SFR

NEED_HEALPIX
TREE_RAD
TREE_RAD_H2
TREE_RAD_CO
NSIDE=1    #NPIX=12*NSIDE*NSIDE
#OUTPUTCOL
#TREE_RAD_NO_GRAVITY
#SELFGRAVITY_OFF                # turn off self-gravity (compatible with GRAVITY_ANALYTIC); setting NOGRAVITY gives identical functionality

#DGR_SCALE_WITH_Z
#BS19_XH2
#CR_HEATING_DM72

COOLING
COOLING_OPERATOR_SPLIT
METALS

NO_ISEND_IRECV_IN_DOMAIN       # MPI debugging: slower, but fixes memory errors during exchange in the domain decomposition (ANY RUN with >2e9 particles MUST SET THIS OR FAIL!)


#DUST_ONE_FLUID             #including dust with a one-fluid approach (passive scalar)
#N_DUST_ELEM=2                                                                                                                                                  
#DUST_EVOLUTION
#DUST_IN_AGB
#DUST_GROWTH_IN_ISM
#MULTI_GRAIN_SIZE=8  #number of grains sizes
#GRAIN_SIZE_MRN
#GRAIN_SIZE_MRN_NO_SHUFFLE


#METALS_IC_TURBBOX_MIX #set metallicity in ICs for turb box (right=1 left=0.1)
#LES_FILTER                 # calculate filtered quantities for an a priori test in LES

#COOLING_ONEZONE_TEST

#RANDOM_SN_INJECT
#SINGLE_SN_INJECT            # need RANDOM_SN_INJECT
#EQUAL_WEIGHT_SN_INJECT

####################################################################################################
# --------------------------------------- Boundary Conditions & Dimensions
####################################################################################################
#BOX_PERIODIC               # Use this if periodic boundaries are needed (otherwise open boundaries are assumed)
#BOX_BND_PARTICLES          # particles with ID=0 are forced in place (their accelerations are set =0): use for special boundary conditions where these particles represent fixed "walls"
#BOX_LONG_X=140             # modify box dimensions (non-square periodic box): multiply X (BOX_PERIODIC and SELFGRAVITY_OFF required)
#BOX_LONG_Y=1               # modify box dimensions (non-square periodic box): multiply Y
#BOX_LONG_Z=1               # modify box dimensions (non-square periodic box): multiply Z
#BOX_REFLECT_X              # make the x-boundary reflecting (assumes a box 0<x<1, unless BOX_PERIODIC is set)
#BOX_REFLECT_Y              # make the y-boundary reflecting (assumes a box 0<y<1, unless BOX_PERIODIC is set)
#BOX_REFLECT_Z              # make the z-boundary reflecting (assumes a box 0<z<1, unless BOX_PERIODIC is set)
#BOX_SHEARING=1             # shearing box boundaries: 1=r-z sheet (r,z,phi coordinates), 2=r-phi sheet (r,phi,z), 3=r-phi-z box, 4=as 3, with vertical gravity
#BOX_SHEARING_Q=(3./2.)     # shearing box q=-dlnOmega/dlnr; will default to 3/2 (Keplerian) if not set
#BOX_SPATIAL_DIMENSION=3    # sets number of spatial dimensions evolved (default=3D). Switch for 1D/2D test problems: if =1, code only follows the x-line (all y=z=0), if =2, only xy-plane (all z=0). requires SELFGRAVITY_OFF
####################################################################################################



####################################################################################################
# --------------------------------------- Hydro solver method
####################################################################################################
# --------------------------------------- Finite-volume Godunov methods (choose one, or SPH)
HYDRO_MESHLESS_FINITE_MASS     # solve hydro using the mesh-free Lagrangian (fixed-mass) finite-volume Godunov method
#HYDRO_MESHLESS_FINITE_VOLUME   # solve hydro using the mesh-free (quasi-Lagrangian) finite-volume Godunov method (control mesh motion with HYDRO_FIX_MESH_MOTION)
#HYDRO_REGULAR_GRID             # solve hydro equations on a regular (recti-linear) Cartesian mesh (grid) with a finite-volume Godunov method
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- Options to explicitly control the mesh motion (for use with the MFV or grid solvers): only set for non-standard behavior
#HYDRO_FIX_MESH_MOTION=0        # mesh with arbitrarily-defined mesh-generating velocities: (0=non-moving, 1=fixed-v [set in ICs] cartesian, 2=fixed-v [ICs] cylindrical, 3=fixed-v [ICs] spherical, 4=analytic function, 5=smoothed-Lagrangian, 6=glass-generating, 7=fully-Lagrangian)
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- SPH methods (enable one of these flags to use SPH):
#HYDRO_PRESSURE_SPH             # solve hydro using SPH with the 'pressure-sph' formulation ('P-SPH')
#HYDRO_DENSITY_SPH              # solve hydro using SPH with the 'density-sph' formulation (GADGET-2 & GASOLINE SPH)
# --------------------------------------- SPH artificial diffusion options (use with SPH; not relevant for Godunov/Mesh modes)
#SPH_DISABLE_CD10_ARTVISC       # for SPH only: Disable Cullen & Dehnen 2010 'inviscid sph' (viscosity suppression outside shocks); just use Balsara switch
#SPH_DISABLE_PM_CONDUCTIVITY    # for SPH only: Disable mixing entropy (J.Read's improved Price-Monaghan conductivity with Cullen-Dehnen switches)
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- Kernel Options
#KERNEL_FUNCTION=3              # Choose the kernel function (2=quadratic peak, 3=cubic spline [default], 4=quartic spline, 5=quintic spline, 6=Wendland C2, 7=Wendland C4, 8=2-part quadratic)
####################################################################################################



####################################################################################################
# --------------------------------------- Additional Fluid Physics
####################################################################################################
## ----------------------------------------------------------------------------------------------------
# --------------------------------------- Gas (or Material) Equations-of-State
#ISOTHERMAL

####################################################################################################
## ------------------------ Gravity & Cosmological Integration Options ---------------------------------
####################################################################################################
# --------------------------------------- TreePM Options (recommended for cosmological sims)
#PMGRID=1024                     # adds Particle-Mesh grid for faster (but less accurate) long-range gravitational forces: value sets resolution (e.g. a PMGRID^3 grid will overlay the box, as the 'top level' grid)
#PM_PLACEHIGHRESREGION=1+2+16   # adds a second-level (nested) PM grid before the tree: value denotes particle types (via bit-mask) to place high-res PMGRID around. Requires PMGRID.
#PM_HIRES_REGION_CLIPPING=1000  # optional additional criterion for boundaries in 'zoom-in' type simulations: clips gas particles that escape the hires region in zoom/isolated sims, specifically those whose nearest-neighbor distance exceeds this value (in code units)
#PM_HIRES_REGION_CLIPDM         # split low-res DM particles that enter high-res region (completely surrounded by high-res)
## -----------------------------------------------------------------------------------------------------
# ---------------------------------------- Adaptive Grav. Softening (including Lagrangian conservation terms!)
#ADAPTIVE_GRAVSOFT_FORGAS       # allows variable softening length for gas particles (scaled with local inter-element separation), so gravity traces same density field seen by hydro
#ADAPTIVE_GRAVSOFT_FORALL=100   # enable adaptive gravitational softening lengths for all particle types (ADAPTIVE_GRAVSOFT_FORGAS should be disabled). the softening is set to the distance
                                # enclosing a neighbor number set in the parameter file. baryons search for other baryons, dm for dm, sidm for sidm, etc. If set to numerical value, the maximum softening is this times All.ForceSoftening[for appropriate particle type]. cite Hopkins et al., arXiv:1702.06148
## -----------------------------------------------------------------------------------------------------
#GRAVITY_NOT_PERIODIC           # self-gravity is not periodic, even though the rest of the box is periodic
## -----------------------------------------------------------------------------------------------------
#GRAVITY_ANALYTIC               # specific analytic gravitational force to use instead of or with self-gravity. If set to a numerical value
                                #  > 0 (e.g. =1), then BH_CALC_DISTANCES will be enabled, and it will use the nearest BH particle as the center for analytic gravity computations
                                #  (edit "gravity/analytic_gravity.h" to actually assign the analytic gravitational forces). 'ANALYTIC_GRAVITY' gives same functionality
#GALSF_GENERATIONS=1             # the number of star particles a gas particle may spawn (defaults to 1, set otherwise)

####################################################################################################
# --------------------------------------- Multi-Threading and Parallelization options
####################################################################################################
#OPENMP=2                       # Masterswitch for explicit OpenMP implementation
#PTHREADS_NUM_THREADS=4         # custom PTHREADs implementation (don't enable with OPENMP)
MULTIPLEDOMAINS=4             # Multi-Domain option for the top-tree level (alters load-balancing)
####################################################################################################



####################################################################################################
# --------------------------------------- Input/Output options
####################################################################################################
#OUTPUT_ADDITIONAL_RUNINFO      # enables extended simulation output data (can slow down machines significantly in massively-parallel runs)
#OUTPUT_IN_DOUBLEPRECISION      # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION       # input files assumed to be in double precision (otherwise float is assumed)
#OUTPUT_POSITIONS_IN_DOUBLE     # input/output files in single, but positions in double (used in hires, hi-dynamic range sims when positions differ by < float accuracy)
#INPUT_POSITIONS_IN_DOUBLE      # as above, but specific to the ICs file
#OUTPUT_POTENTIAL               # forces code to compute+output potentials in snapshots
#OUTPUT_ACCELERATION            # output physical acceleration of each particle in snapshots
#OUTPUT_CHANGEOFENERGY          # outputs rate-of-change of internal energy of gas particles in snapshots
#OUTPUT_VORTICITY               # outputs the vorticity vector
OUTPUT_TIMESTEP                # outputs timesteps for each particle
#OUTPUT_COOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_LINEOFSIGHT				# enables on-the-fly output of Ly-alpha absorption spectra
#OUTPUT_LINEOFSIGHT_SPECTRUM    # computes power spectrum of these (requires additional code integration)
#OUTPUT_LINEOFSIGHT_PARTICLES   # computes power spectrum of these (requires additional code integration)
#OUTPUT_POWERSPEC               # compute and output power spectra (not used)
#OUTPUT_RECOMPUTE_POTENTIAL     # update potential every output even it EVALPOTENTIAL is set
#INPUT_READ_HSML                # force reading hsml from IC file (instead of re-computing them; in general this is redundant but useful if special guesses needed)
#OUTPUT_TWOPOINT_ENABLED        # allows user to calculate mass 2-point function by enabling and setting restartflag=5
#IO_DISABLE_HDF5                # disable HDF5 I/O support (for both reading/writing; use only if HDF5 not install-able)
####################################################################################################



####################################################################################################
# -------------------------------------------- De-Bugging & special (usually test-problem only) behaviors
####################################################################################################
# --------------------
# ----- General De-Bugging and Special Behaviors
#DEVELOPER_MODE                 # allows you to modify various numerical parameters (courant factor, etc) at run-time
#STOP_WHEN_BELOW_MINTIMESTEP    # forces code to quit when stepsize wants to go below MinSizeTimestep specified in the parameterfile
#DEBUG                          # enables core-dumps and FPU exceptions
# --------------------
# ----- Hydrodynamics
#FREEZE_HYDRO                   # zeros all fluxes from RP and doesn't let particles move (for testing additional physics layers)
#EOS_ENFORCE_ADIABAT=(1.0)      # if set, this forces gas to lie -exactly- along the adiabat P=EOS_ENFORCE_ADIABAT*(rho^GAMMA)
#SLOPE_LIMITER_TOLERANCE=1      # sets the slope-limiters used. higher=more aggressive (less diffusive, but less stable). 1=default. 0=conservative. use on problems where sharp density contrasts in poor particle arrangement may cause errors. 2=same as AGGRESSIVE_SLOPE_LIMITERS below
#AGGRESSIVE_SLOPE_LIMITERS      # use the original GIZMO paper (more aggressive) slope-limiters. more accurate for smooth problems, but
                                # these can introduce numerical instability in problems with poorly-resolved large noise or density contrasts (e.g. multi-phase, self-gravitating flows)
#ENERGY_ENTROPY_SWITCH_IS_ACTIVE # enable energy-entropy switch as described in GIZMO methods paper. This can greatly improve performance on some problems where the
                                # the flow is very cold and highly super-sonic. it can cause problems in multi-phase flows with strong cooling, though, and is not compatible with non-barytropic equations of state
#FORCE_ENTROPIC_EOS_BELOW=(0.01) # set (manually) the alternative energy-entropy switch which is enabled by default in MFM/MFV: if relative velocities are below this threshold, it uses the entropic EOS
#DISABLE_SPH_PARTICLE_WAKEUP    # don't let gas particles move to lower timesteps based on neighbor activity (use for debugging)
# --------------------
# ----- Additional Fluid Physics and Gravity
#MHD_ALTERNATIVE_LEAPFROG_SCHEME # use alternative leapfrog where magnetic fields are treated like potential/positions (per Federico Stasyszyn's suggestion): still testing
#SUPER_TIMESTEP_DIFFUSION       # use super-timestepping to accelerate integration of diffusion operators [for testing or if there are stability concerns]
#EVALPOTENTIAL                  # computes gravitational potential
# --------------------
# ----- Particle IDs
#TEST_FOR_IDUNIQUENESS          # explicitly check if particles have unique id numbers (only use for special behaviors)
#LONGIDS                        # use long ints for IDs (needed for super-large simulations)
#ASSIGN_NEW_IDS                 # assign IDs on startup instead of reading from ICs
#NO_CHILD_IDS_IN_ICS            # IC file does not have child IDs: do not read them (used for compatibility with snapshot restarts from old versions of the code)
# --------------------
# ----- Particle Merging/Splitting/Deletion/Boundaries
#PREVENT_PARTICLE_MERGE_SPLIT   # don't allow gas particle splitting/merging operations
#PARTICLE_EXCISION              # enable dynamical excision (remove particles within some radius)
#MERGESPLIT_HARDCODE_MAX_MASS=(1.0e-6)   # manually set maximum mass for particle merge-split operations (in code units): useful for snapshot restarts and other special circumstances
#MERGESPLIT_HARDCODE_MIN_MASS=(1.0e-7)   # manually set minimum mass for particle merge-split operations (in code units): useful for snapshot restarts and other special circumstances
# --------------------
# ----- MPI & Parallel-FFTW De-Bugging
#USE_MPI_IN_PLACE               # MPI debugging: makes AllGatherV compatible with MPI_IN_PLACE definitions in some MPI libraries
#FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG # MPI debugging
#MPISENDRECV_SIZELIMIT=100      # MPI debugging
#MPISENDRECV_CHECKSUM           # MPI debugging
#DONOTUSENODELIST               # MPI debugging
#NOTYPEPREFIX_FFTW              # FFTW debugging (fftw-header/libraries accessed without type prefix, adopting whatever was
                                #   chosen as default at compile of fftw). Otherwise, the type prefix 'd' for double is used.
#DOUBLEPRECISION_FFTW           # FFTW in double precision to match libraries
# --------------------
# ----- Load-Balancing
#ALLOW_IMBALANCED_GASPARTICLELOAD # increases All.MaxPartSph to All.MaxPart: can allow better load-balancing in some cases, but uses more memory. But use me if you run into errors where it can't fit the domain (where you would increase PartAllocFac, but can't for some reason)
#SEPARATE_STELLARDOMAINDECOMP   # separate stars (ptype=4) and other non-gas particles in domain decomposition (may help load-balancing)
####################################################################################################
















