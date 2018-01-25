//! [SYSTEM]
SYSTEM
# NPM : protein molecules (0, 1) 
# NSM : solvent molecules (>= 0) 
#      NPM      NSM 
         1        0 
END 
//! [SYSTEM]

//! [ENERGYMIN]
ENERGYMIN
# NTEM: 0..1 controls energy minimisation mode 
#       0: do not do energy minimisation (default) 
#       1: steepest-descent minimisation 
# NCYC: >0 number of steps before resetting of conjugate-gradient search direction (not in use !!) 
# DELE: >0.0 energy threshold for convergence 
# DX0: > 0.0 initial step size 
# DXM: > 0.0 maximum step size 
# NMIN > 0 minimum number of minimisation steps 
# FLIM >= 0.0 limit force to maximum value (FLIM > 0.0 is not recommended). 
#     NTEM    NCYC    DELE    DX0     DXM    NMIN    FLIM 
         1       0     0.1   0.01    0.05       1       0 
END 
//! [ENERGYMIN]

//! [STEP]
STEP
#   NSTLIM  >0 number of steps 
#   T       >=0 initial time 
#           -1  read time from configuration file 
#   DT      >0 time step 
# 
#   NSTLIM         T        DT 
       100       0.0     0.005 
END 
//! [STEP]

//! [CONSTRAINT]
CONSTRAINT
#	NTC 
#		1	solvent    solvent only 
#		2	hydrogen   solvent and solute bonds containing hydrogens and  
#		               constraints in the topology's CONSTRAINTS block 
#		3	all        solvent and solute all bonds 
#		4	specified  solvent and constraints in the topology's CONSTRAINTS block 
    3 
#       NTCP: solute algorithm 
#               1        shake 
#               2        lincs 
#               3        flexshake 
#       NTCP 
        1 
#       NTCP0(1)..(3): algorithm options 
#         - shake: tolerance 
#         - lincs: order 
#         - flexshake: tolerance, readin, order 
#       NTCP0(1)   NTCP0(2)   NTCP0(3) 
        0.0001 
#	NTCS: solvent algorithm 
#               1        shake 
#               2        lincs 
#               3        flexshake 
#               4        settle 
#               5        m_shake (only implemented for water and methanol!) 
#               6        gpu_shake 
#       NTCS 
        1 
#       NTCS0(1):  algorithm options 
#         - shake: tolerance 
#         - lincs: order 
#         - flexshake: tolerance, readin, order 
#         - settle: no arguments 
#         - m_shake: tolerance 
#       NTCS0(1) 
        0.0001 
#       NTCG: Number of GPUs 
#       NTCD: Device number of the GPU 
#       NTCG          NTCD 
        1             0 
END 
//! [CONSTRAINT]

//! [PRINTOUT]
PRINTOUT
#  NTPR: print out energies, etc. every NTPR steps 
#  NTPP: =1 perform dihedral angle transition monitoring 
#     NTPR      NTPP 
         0         0 
END 
//! [PRINTOUT]

//! [WRITETRAJ]
WRITETRAJ
# NTWX       controls writing of coordinate trajectory 
#       0: no coordinate trajectory is written (default) 
#      >0: write solute and solvent coordinates every NTWX steps 
#      <0: write solute coordinates every |NTWX| steps 
# NTWSE >= 0 selection criteria for coordinate trajectory writing 
#       0: write normal coordinate trajectory 
#      >0: write minimum-energy coordinate and energy trajectory (based on the 
#          energy entry selected by NTWSE and as blocks of length NTWX) 
#          (see configuration/energy.cc or ene_ana library for indices) 
# NTWV       controls writing of velocity trajectory 
#       0: no velocity trajectory is written (default) 
#      >0: write solute and solvent velocities every NTWV steps 
#      <0: write solute velocities every |NTWV| steps 
# NTWF       controls writing of force trajectory 
#       0: no force trajectory is written (default) 
#      >0: write solute and solvent forces every NTWF steps 
#      <0: write solute forces every |NTWF| steps 
# NTWE >= 0 controls writing of energy trajectory 
#       0: no energy trajectory is written (default) 
#      >0: write energy trajectory every NTWE steps 
# NTWG >= 0 controls writing of free energy trajectory 
#       0: no free energy trajectory is written (default) 
#      >0: write free energy trajectory every NTWG steps 
# NTWB >= 0 controls writing of block-averaged energy trajectory 
#       0: no block averaged energy trajectory is written (default) 
#      >0: write block-averaged energy variables every |NTWB| steps 
#          (and free energies if NTWG > 0) trajectory 
# 
#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB 
      100          0         0         0       100         0       100 
END 
//! [WRITETRAJ]

//! [PRESSURESCALE]
PRESSURESCALE
#	COUPLE:	off(0), calc(1), scale(2) 
#	SCALE:  off(0), iso(1), aniso(2), full(3), semianiso(4) 
#	VIRIAL: none(0), atomic(1), group(2) 
# 
#   COUPLE  SCALE   COMP        TAUP    VIRIAL 
    calc    iso     4.575E-4    0.5     atomic 
#   SEMI (semianisotropic couplings: X, Y, Z) 
#       e.g. 1 1 2: x and y jointly coupled and z separately coupled 
#       e.g. 0 0 1: constant area (xy-plane) and z coupled to a bath 
    1 1 2 
#   reference pressure 
    0.06102     0.00000     0.00000 
    0.00000     0.06102     0.00000 
    0.00000     0.00000     0.06102 
END 
//! [PRESSURESCALE]

//! [BOUNDCOND]
BOUNDCOND
#  NTB: boundary conditions 
#       -1 : truncated octahedral 
#        0 : vacuum 
#        1 : rectangular 
#        2 : triclinic 
#  NDFMIN: number of degrees of freedom subtracted for temperature 
# 
#         NTB    NDFMIN 
            1         0 
END 
//! [BOUNDCOND]

//! [PERTURBATION]
PERTURBATION
#    NTG: 0..1 controls use of free-energy calculation. 
#         0: no free-energy calculation (default) 
#         1: calculate dH/dRLAM 
#  NRDGL: 0,1 controls reading of initial value for RLAM. 
#         0: use initial RLAM parameter from PERTURBATION block 
#         1: read from configuration 
#   RLAM: 0.0..1.0 initial value for lambda 
#  DLAMT: >= 0.0 rate of lambda increase in time. 
# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter 
#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter 
#   NLAM: > 0 power dependence of lambda coupling 
# NSCALE: 0..2 controls use of interaction scaling 
#         0: no interaction scaling 
#         1: interaction scaling 
#         2: perturbation for all atom pairs with scaled 
#            interactions. No perturbation for others. 
# 
#     NTG   NRDGL    RLAM   DLAMT 
        0       0     0.0     0.0 
#  ALPHLJ   ALPHC    NLAM  NSCALE 
      0.0     0.0       1       0 
END 
//! [PERTURBATION]

//! [FORCE]
FORCE
# NTF(1..6): 0,1 determines terms used in force calculation 
#             0: do not include terms 
#             1: include terms 
# NEGR: ABS(NEGR): number of energy groups 
#             > 0: use energy groups 
#             < 0: use energy and force groups 
# NRE(1..NEGR): >= 1 last atom in each energy group 
# NTF(1)    NTF(2)    NTF(3)    NTF(4)    NTF(5)        NTF(6) 
# bonds     angles    improper  dihedral  electrostatic vdW 
  0         1         1         1         1             1 
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR) 
     1        60 
END 
//! [FORCE]

//! [COVALENTFORM]
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential 
#        0: quartic form (default) 
#        1: harmonic form 
# NTBAH: 0,1 controls bond-angle bending potential 
#        0: cosine-harmonic (default) 
#        1: harmonic 
# NTBDN: 0,1 controls torsional dihedral potential 
#        0: arbitrary phase shifts (default) 
#        1: phase shifts limited to 0 and 180 degrees. 
#   NTBBH    NTBAH    NTBDN 
        0        0        0 
END 
//! [COVALENTFORM]

//! [INITIALISE]
INITIALISE
# NTIVEL: 0,1 controls generation of initial velocities. 
#         0: read from configuration (default) 
#         1: generate from Maxell distribution at temperature TEMPI 
# NTISHK: 0..3 controls shaking of initial configuration 
#         0: no intial SHAKE (default) 
#         1: initial SHAKE on coordinates only 
#         2: initial SHAKE on velocities only 
#         3: initial SHAKE on coordinates and velocities 
# NTINHT: 0,1 controls generation of initial Nose-Hoover chain variables 
#         0: read from configuration (default) 
#         1: reset variables to zero. 
# NTINHB: 0,1 controls generation of initial Nose-Hoover (chain) barostat 
#             variables 
#         0: read from strartup file (if applicable) (default) 
#         1: reset variables to zero 
# NTISHI: 0,1 controls initial setting for lattice shift vectors 
#         0: read from configuration (default) 
#         1: reset shifts to zero. 
# NTIRTC: 0,1 controls initial setting of positions and orientations for 
#             roto-translational constraints 
#         0: read from configuration (default) 
#         1: reset based on initial configuraion of startup file 
# NTICOM: 0,1 controls initial removal of COM motion 
#         0: no initial system COM motion removal (default) 
#         1: initial COM translation is removed 
#         2: initial COM rotation is removed 
# NTISTI: 0,1 controls generation of stochastic integrals 
#         0: read stochastic integrals and IG from configuration (default) 
#         1: set stochastic integrals to zero and use IG from here. 
# IG:     random number generator seed 
# TEMPI:  initial temperature 
# 
#  NTIVEL  NTISHK  NTINHT  NTINHB 
        0       0       0       0 
#  NTISHI  NTIRTC  NTICOM 
        0       0       0 
#  NTISTI 
        0 
#      IG   TEMPI 
        0     0.0 
END 
//! [INITIALISE]

//! [COMTRANSROT]
COMTRANSROT
#    NSCM : controls system centre-of-mass (com) motion removal 
#           0: no com motion removal (default) 
#         < 0: com translation and rotation are removed every abs(NSCM) 
#              steps. 
#         > 0: com translation is removed every NSCM steps. 
#     NSCM 
         0 
END 
//! [COMTRANSROT]

//! [HOOMD]
HOOMD
#       PROCESSOR: cpu gpus 
# 
#       PROCESSOR 
    gpus 
END 
//! [HOOMD]

//! [PAIRLIST]
PAIRLIST
#       ALGORITHM  standard(0) (gromos96 like pairlist) 
#                  grid(1) (md++ grid pairlist) 
#                  grid_cell(2) (creates a mask) 
#       NSNB  >0    frequency (number of steps) a pairlist is constructed 
#       RCUTP >0.0  short-range cut-off in twin-range 
#       RCUTL >0.0  intermediate-range cut-off in twin-range 
#       SIZE  >0    grid cell size (or auto = 0.5 * RCUTP) 
#       TYPE  chargegoup(0) (chargegroup based cutoff) 
#             atomic(1)     (atom based cutoff) 
# 
#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE 
        grid            5       0.8     1.4     auto    chargegroup 
# 
END 
//! [PAIRLIST]

//! [CGRAIN]
CGRAIN
# NTCGRAN 0..3 coarse grain selection 
#         0: atomistic (off) 
#         1: coarse-grained using MARTINI model (on) 
#         2: coarse-grained using GROMOS model (on) 
#         3: mixed-grained using GROMOS model (on) 
#     EPS >= 0.0 dielectric constant for coarse grained coulombic interaction 
#    EPSM >= 0.0 dielectric constant for mixed CG-FG coulombic interaction 
# NTCGRAN     EPS     EPSM 
        1      20        1 
END 
//! [CGRAIN]

//! [MULTIBATH]
MULTIBATH
# NTBTYP: temperature coupling algorithm 
#		weak-coupling(0) 
#		nose-hoover(1) 
#		nose-hoover-chains(2)	num 
#		(where num is the number of chains to use) 
#   NTBTYP           NUM 
    nose-hoover-chains	3 
#   NBATHS: number of baths 
    2 
#   TEMP0  TAU 
    300    0.10 
    300    0.10 
#   DOFSET: number of different couplings 
    1 
#   LAST   COM-BATH  IR-BATH 
    60     1         2 
#   (this couples the first 60 atoms com motion to bath 1 and 
#    the internal / rotational motion to bath 2) 
END 
//! [MULTIBATH]

//! [POSITIONRES]
POSITIONRES
#    NTPOR 0..3 controls atom positions re(con)straining. 
#          0: no position re(con)straints (default) 
#          1: restraining with force constant CPOR 
#          2: restraining with force constant CPOR weighted by 
#             atomic B-factors 
#          3: constraining 
#   NTPORB 0,1 controls reading of reference positions and 
#              B-factors 
#          0: read reference positions from startup file. 
#          1: read reference positions and B-factors from 
#             special file 
#   NTPORS 0,1 controls scaling of reference positions upon 
#              pressure scaling 
#          0: do not scale reference positions 
#          1: scale reference positions 
#     CPOR >= 0 position restraining force constant 
# 
#   NTPOR  NTPORB  NTPORS    CPOR 
        0       0       0   2.5E4 
END 
//! [POSITIONRES]

//! [XRAYRES]
XRAYRES
#    NTXR   -2: time-averaged electron density restraints 
#           -1: instantaneous electron density restraints 
#            0: no xray restraints. 
#            1: instantaneous structure factor restraints 
#            2: time-averaged structure factor restraints 
#            3: biquadratic/timeaveraged structure factor restraints 
#    NTXLE   0: do not perform local elevation 
#            1: do perform local elevation 
#    CXR     >= 0 xray restraining force constant 
#    NTWXR   >= 0 write xray data to output file 
#            0: don't write xray data 
#            > 0 write every NTPXRth step 
#    NTWDE   0..3 write density-maps 
#            0: write nothing 
#            1: write electron densitiy map 
#            2: write asymmetric-unit-only electron densitiy map 
#            3: write both 
#    NTWXM   >= 0 write every NTWXMth step electron density map(s) to external file 
#    CXTAU   >=0 xray time-average restraining memory-constant 
#    RDAVG   0/1 read sf-timeaverages (from job to job) 
# 
#   NTXR   NTXLE     CXR   NTWXR   NTWDE   NTWXM   CXTAU   RDAVG 
       0       0     0.0       0       0      0      0.0       0 
END 
//! [XRAYRES]

//! [DISTANCERES]
DISTANCERES
#   NTDIR -2..2 controls distance restraining 
#         0: no distance restraining (default) 
#         1: instantaneous, using force constant CDIR 
#         2: instantaneous, using force constant CDIR x W0 
#        -1: time-averaged, using force constant CDIR 
#        -2: time-averaged, using force constant CDIR x W0 
#  NTDIRA 0,1 controls values for initial distance averages 
#         0: generate initial averages 
#         1: read from configuration 
#    CDIR >= 0.0 force constant for distance restraining 
#    DIR0 > 0.0 distance offset in restraining function 
#  TAUDIR >= 0.0 coupling time for time averaging 
# FORCESCALE 0..2 controls approximation of force scaling 
#         0: approximate d<r>/dr = 1 
#         1: approximate d<r>/dr = (1.0 - exp(-Dt/tau)) 
#         2: use d<r>/dr = (1.0 - exp(-Dt/tau))*(<r>/r)^4 
#    VDIR 0,1 controls contribution to virial 
#         0: no contribution 
#         1: distance restraints contribute to virial 
#  NTWDIR >= 0 write every NTWDIRth step dist. restr. information to external file 
#   NTDIR  NTDIRA    CDIR    DIR0  TAUDIR  FORCESCALE  VDIR  NTWDIR 
        0       0     0.0     1.0     0.0           0     0       0 
END 
//! [DISTANCERES]

//! [DISTANCEFIELD]
DISTANCEFIELD
#   NTDFR 0,1         controls distance field restraining 
#         0: no distance field restraining (default) 
#         1: apply distance field restraining 
#   GRID  > 0.0       grid size for distance field 
#   PROTEINOFFSET >= 0 penalty for distances through the host 
#   PROTEINCUTOFF >= 0 distance to protein atoms to be considered inside 
#   PROTECT >= 0      protect grid points within this radius around the zero-distance 
#                     point from being flagged as protein 
#   UPDATE > 0        update frequency for grid 
#   RL >= 0           linearize forces for distances larger than RL 
#   SMOOTH >= 0       smoothen the protein boundary after grid construction 
#                     by SMOOTH layers 
#   NTWDF >= 0        write every NTWDF step disfield information to external file 
#   PRINTGRID = 0,1   write grid to final configuration file 
# 
#   NTDFR 
        1 
#    GRID   PROTEINOFFSET  PROTEINCUTOFF  PROTECT 
      0.2   15             0.2            0 
#  UPDATE   SMOOTH   RL    NTWDF   PRINTGRID 
      100   1        1.0      50           0 
END 
//! [DISTANCEFIELD]

//! [DIHEDRALRES]
DIHEDRALRES
# NTDLR   0...3 controls dihedral-angle restraining and constraining 
#         0:    off [default] 
#         1:    dihedral restraining using CDLR 
#         2:    dihedral restraining using CDLR * WDLR 
#         3:    dihedral constraining 
# 
# CDLR    >=0.0 force constant for dihedral restraining [kJ/mol/degree^2] 
# PHILIN  >0.0  deviation after which the potential energy function is linearized 
# NTWDLR  >=0   write every NTWDLR step dihedral information to external file 
# 
# NTDLR  CDLR      PHILIN  NTWDLR 
  1      100.0     180.0   100 
END 
//! [DIHEDRALRES]

//! [JVALUERES]
JVALUERES
# NTJVR    -3..2 
#          -3                biquadratic using CJVR * WJVR 
#          -2                time-averaged using CJVR * WJVR 
#          -1                time-avaraged using CJVR 
#           0                no J-value restraints [default] 
#           1                instantaneous using CJVR 
#           2                instantaneous using CJVR * WJVR 
# NTJVRA    0                controls reading of averages from startup file 
#           0                start from initial values of J0 [default] 
#           1                read time averages from startup file (for continuation time-averaged run) 
# CJVR   >= 0                J-value restraining force constant 
#                            (weighted by individual WJVR) 
# TAUJVR >= 0                coupling time for time-averaging 
# NJVRTARS  0,1              omits or includes force scaling by memory decay factor in case of time-averaging 
#           0                omit factor (set (1 - exp(-Dt/tau)) = 1) 
#           1                scale force by (1 - exp(-Dt/tau)) 
# NJVRBIQW  0..2             controls weights (X in Eq. 19 of MD98.17) of the two terms in biquadratic restraining 
#           0                X = 1 
#           1                X = (1 - exp(-Dt/tau)) 
#           2                X = 0 
# LE        0,1              local elevation restraining 
#           0                local elevation off [default] 
#           1                local elevation on 
# NGRID   > 1                number of grid points in local elevation restraining 
# DELTA  >= 0.0              no elevation of potential if J is within DELTA of J0 
# NTWJV  >= 0                write J-value averages and LE grid to special trajectory 
#           0                don't write [default] 
#         > 0                write every NTWJVth step 
# 
#       NTJVR  NTJVRA  CJVR   TAUJVR  NJVRTARS   NJVRBIQW   LE    NGRID   DELTA  NTWJV 
           -3  0       10.0      5.0     0          0       1       16     0.5      0 
END 
//! [JVALUERES]

//! [ORDERPARAMRES]
ORDERPARAMRES
# NTOPR    -2..2 
#          -2                time-averaged using COPR * WOPR 
#          -1                time-averaged using COPR 
#           0                no order-parameter restraints [default] 
#           1                window-averaged using COPR 
#           2                window-averaged using COPR * WOPR 
# NTOPRA    0                controls reading of averages from startup file 
#           0                start from initial values of S0 [default] 
#           1                read time averages from startup file (for continuation time-averaged run) 
# COPR   >= 0.0              order-parameter restraining force constant 
#                            (weighted by individual WOPR) 
# TAUOPR >= 0.0              coupling time for time-averaging 
# UPDOPR  > 0                update average every UPDOPRth step 
# NTWOP  >= 0                write order-parameter to special trajectory 
#           0                don't write [default] 
#         > 0                write every NTWOP step 
# 
#       NTOPR  NTOPRA  COPR   TAUOPR   UPDOPR  NTWOP 
           -2  0       10.0      5.0        1      0 
END 
//! [ORDERPARAMRES]

//! [RDCRES]
RDCRES
# NTRDCR -4..2                 RDC restraining 
#        -4:                   biquadratic using CRDCR * WRDCR 
#        -3:                   biquadratic using CRDCR 
#        -2:                   time averaged using CRDCR * WRDCR 
#        -1:                   time averaged using CRDCR 
#         0:                   no RDC restraints [default] 
#         1:                   instantaneous using CRDCR 
#         2:                   instantaneous using CRDCR * WRDCR 
# NTRDCRA 0,1                  controls reading of average RDCs 
#         0:                   take initial values RDC0 from the RDC restraint file 
#         1:                   read time averages from initial coordinate file 
#                              (for continuation run) 
# 
# NTRDCT  0..2                 Type of alignment representation 
#         0:                   cartesian magnetic field vectors 
#         1:                   alignment tensor 
#         2:                   spherical harmonics 
# NTALR   0,1                  controls reading of values in the chosen representation 
#         0:                   start from values given in RDC restraint file 
#         1:                   read values from initial coordinate file (for continuation run) 
# 
# METHOD  0..2                 Method of updating the magnetic field vectors 
#         0:                   Energy minimisation 
#         1:                   Stochastic dynamics 
#         2:                   Molecular dynamics 
# EMGRAD  > 0.0                (METHOD = 0, EM) stop minimisation if gradient is below EMGRAD 
# EMDX0   > 0.0                (METHOD = 0, EM) initial step size 
# EMNMAX  > 0                  (METHOD = 0, EM) maximum number of minimisation steps 
# SDCFRIC >= 0.0               (METHOD = 1, SD) global friction coefficient gamma 
# TEMP  >= 0.0                 temperature of stochastic bath (SD) and temperature used for initial velocities (MD, SD) 
# DELTA   >= 0                 the flatbottom potential is 2 DELTA wide [ps^-1] 
# CRDCR   >= 0                 RDC restraining force constant [kJ*ps^2] 
#                              (weighted by individual WRDCR) 
# TAU     >= 0                 coupling time for time averaging [ps] 
# NRDCRTARS 0,1                omits or includes force scaling by memory decay factor in case of time-averaging 
#           0                  omit factor (set (1-exp(-dt/tau))=1 ) 
#           1                  scale force by (1-exp(-dt/tau)) 
# NRDCRBIQW 0..2               controls weights of the two terms in biquadratic restraining 
#           0                  X = 1 
#           1                  X = (1 - exp(-dt/tau)) 
#           2                  X = 0 
# NTWRDC   >= 0                write output to special trajectory 
#           0:                 don't write 
#          >0:                 write every NTWRDCth step. 
# 
#      NTRDCR  NTRDCRA  NTRDCT  NTALR  METHOD 
            2        0       0      0       0 
#      EMGRAD  EMDX0  EMNMAX  SDCFRIC    TEMP    DELTA  CRDCR  TAU   NRDCRTARS NRDCRBIQW   NTWRDC   NTWRDC 
        0.001   0.01    1000       20     300      0      1     1           0    0          0        10000 
END 
//! [RDCRES]

//! [PERSCALE]
PERSCALE
# RESTYPE		special energy term to which periodic scaling should 
#            be applied 
#	   0		don't apply periodic scaling 
#	   1		apply periodic scaling to J-value restraints 
# 
# parameters for periodic scaling of J-value restraints 
# KDIH	>= 0		maximum scaling factor for dihedral angle potential 
# KJ	>= 0		maximum scaling factor for J-Value restraint potential 
# T	>= 0		period of cosine scaling function 
# DIFF	>= 0		minimum deviation from target value at which to start 
#                scaling period 
# RATIO	>= 0		minimum fraction of T that needs to be passed before 
#                starting a new scaling period 
# READ	   0,1		controls reading of scaling parameters 
#          0		reset scaling parameters 
#          1		read from configuration 
# 
# RESTYPE 
      0 
#    KDIH      KJ       T   DIFF    RATIO    READ 
      0.1     0.1     0.2    0.7      1.0       0 
END 
//! [PERSCALE]

//! [ROTTRANS]
ROTTRANS
# roto-translational constraints 
# use either centre of mass removal or roto-translational constraints 
# not both! 
# 
#     RTC: 0,1 controls roto-translational constraints 
#          0 don't use roto-translational constraints (default) 
#          1 use roto-translational constraints 
# RTCLAST: last atom to be roto-translationally constrained 
#     RTC  RTCLAST 
        1     1155 
END 
//! [ROTTRANS]

//! [INNERLOOP]
INNERLOOP
# NTILM: 0..4, acceleration method used 
#        0: use standard solvent loops [default] 
#        1: use fast generic solvent loops 
#        2: use solvent loops with hardcoded parameters 
#        3: use solvent loops with tabulated forces and energies 
#        4: use solvent loops with CUDA library 
# NTILS: 0..1, solvent used 
#        0: use topology [default] 
#        1: use SPC 
# NGPUS: number of GPUs to use 
# NDEVG: Which GPU device number to use. If not given driver will determine. 
# NTILM NTILS NGPUS NDEVG 
      4     0     2   0 1 
END 
//! [INNERLOOP]

//! [REPLICA]
REPLICA
#     NRET >= 1 number of replica exchange temperatures 
#    RET() >= 0.0 temperature for each replica 
# LRESCALE 0,1 controls temperature scaling 
#          0 don't scale temperatures after exchange trial 
#          1 scale temperatures after exchange trial 
#   NRELAM >= 1 number of replica exchange lambda values 
#  RELAM() >= 0.0 lambda value for each lambda-replica 
#   RETS() >= 0.0 timestep of each lambda-replica 
# NRETRIAL >= 0 number of overall exchange trials 
#  NREQUIL >= 0 number of exchange periods to equilibrate 
#               (disallow switches) 
#     CONT >= 0 continuation run 
#             0 start from one configuration file 
#             1 start from multiple configuration files 
# 
# NRET 
  10 
# RET(1..NRET) 
  300.0  320.0  340.0 360.0 380.0 
  400.0  420.0  440.0 460.0 480.0 
# LRESCALE 
  1 
# NRELAM 
  10 
# RELAM(1..NRELAM) 
  0.0    0.1    0.2   0.3   0.4 
  0.5    0.6    0.7   0.8   0.9 
# RETS(1..NRELAM) 
  0.002  0.001  0.001 0.001 0.002 
  0.003  0.002  0.001 0.001 0.002 
# NERTRIAL 
  100 
# NREQUIL 
  10 
# CONT 
  0 
END 
//! [REPLICA]

//! [MULTICELL]
MULTICELL
#  NTM: 0,1 switch for multiple-unit-cell simulation. 
#       0 : single-unit-cell simulation [default] 
#       1 : multiple-unit-cell simulation 
#         NTM 
            0 
#  number of subdivisions along axis 
#   NCELLA    NCELLB    NCELLC 
         1         1         1 
#  periodicity checks (relative tolerance) 
#  not available in md++ -> 0.0 
#    TOLPX     TOLPV     TOLPF    TOLPFW 
       0.0       0.0       0.0       0.0 
END 
//! [MULTICELL]

//! [READTRAJ]
READTRAJ
# NTRD  0,1 controls trajectory-reevaluation mode 
#       0: do not use trajectory-reevaluation mode (default) 
#       1: use trajectory-reevaluation mode 
# NTRN  number of files (ignored) 
# NTRB  read box (must be 1) 
# NTSHK 0,1 controls SHAKE on old coordinates 
#       0 perform SHAKE with respect to previous coordinates 
#       1 perform SHAKE with respect to current coordinates 
# 
#   NTRD    NTRN    NTRB   NTSHK 
       0       0       1       0 
END 
//! [READTRAJ]

//! [INTEGRATE]
INTEGRATE
#  NINT 0..1 selects integration method 
#	0: no integration performed 
#	1: leap-frog integration scheme performed (default) 
# 
#    NINT 
        1 
END 
//! [INTEGRATE]

//! [STOCHDYN]
STOCHDYN
# NTSD    0,1 controls stochastic dynamics mode 
#         0: do not do stochastic dynamics (default) 
#         1: do stochastic dynamics 
# NTFR    0..3 defines atomic friction coefficients gamma 
#         0: set gamma to 0.0 (default) 
#         1: set gamma to CFRIC 
#         2: set gamma to CFRIC*GAM0 
#         3: set gamma to CFRIC*w where w approximates the solvent-accessible  
#            surface area as described in the Stochastic Dynamics Chapter in Vol.2 of the manual  
# NSFR    > 0 recalculate gamma every NSFR steps 
# NBREF   > 0 threshold number of neighbour atoms for a buried atom 
# RCUTF   >= 0.0 interatomic distance considered when calculating gamma 
# CFRIC   >= 0.0 global weighting for gamma 
# TEMPSD  >= 0.0 temperature of stochastic bath 
# 
#     NTSD     NTFR     NSFR   NBREF  RCUTF    CFRIC    TEMPSD 
         0        1        0       6    0.3     91.0     300.0 
END 
//! [STOCHDYN]

//! [EWARN]
EWARN
# MAXENER issue a warning if total energy is larger than this value 
# 
# MAXENER 
   100000 
END 
//! [EWARN]

//! [MULTISTEP]
MULTISTEP
#   STEPS calculate non-bonded every STEPSth step. 
#   BOOST 0,1 
#         0: stored forces of STEPSth step are added every step 
#         1: stored forces of STEPSth setp are multiplied by STEPS 
#            and added every STEPSth step (default) 
# 
#   STEPS   BOOST 
        0       0 
END 
//! [MULTISTEP]

//! [CHEMICALMONTECARLO]
CHEMICALMONTECARLO
# 
#     MC  MCSTEPS   MCDLAM 
       0        1      0.5 
END 
//! [CHEMICALMONTECARLO]

//! [POLARISE]
POLARISE
# COS      0,1,2 use polarisation 
#          0: don't use polarisation (default) 
#          1: use charge-on-spring model for dipolar polarisation 
#          2: use charge-on-spring model for dipolar polarisation with off atom site 
# EFIELD   0,1 controls evaluation site for electric field 
#          0: evaluate at atom position 
#          1: evaluate at cos position 
# MINFIELD >0.0 convergence criterium 
# DAMP     0,1 controls polarisability damping 
#          0: don't damp polarisability 
#          1: damp polarisability (with paramaters from topology) 
# WRITE    > 0 write COS positions to special trajectory 
#          0: don't write 
#         >0: write COS positions every WRITEth step 
# 
#     COS  EFIELD MINFIELD    DAMP  WRITE 
        0       0      2.5       0      0 
END 
//! [POLARISE]

//! [RANDOMNUMBERS]
RANDOMNUMBERS
# NTRNG 0,1 random number generator 
#         0 use G96 algorithm (default) 
#         1 use GSL library 
# NTGSL -1.. GSL random number generation algorithm 
#         -1: use default algorithm (mt19937) 
#       >=0 : run contrib/rng_gsl for a list of possible arguments 
# 
#   NTRNG   NTGSL 
        1      -1 
END 
//! [RANDOMNUMBERS]

//! [EDS]
EDS
# EDS        0,1 
#              0: no enveloping distribution sampling (EDS) [default] 
#              1: enveloping distribution sampling 
# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter 
#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter 
# FORM       1-3 
#              1: Single s Hamiltonian 
#              2: Hamiltonian with NUMSTATES*(NUMSTATES-1)/2 (pairwise) S parameters 
#              3: Hamiltonian with (NUMSTATES-1) S parameters 
# NUMSTATES >1  : number of states 
# if NUMSTATES != 3: 
# S         >0.0: smoothness parameter(s) 
# if NUMSTATES == 3: 
# i   j   S     : state pair i j and associated s parameter 
# EIR           : energy offsets for states 
# 
# EDS 
  1 
# ALPHLJ  ALPHC  FORM  NUMSTATES 
  0.0     0.0       2          3 
# S 
  0.2  0.01 0.1 
# EIR 
  0   20   40 
# 
# ---- OR: example for FORM = 3: 
# 
# EDS 
  1 
# ALPHLJ  ALPHC  FORM  NUMSTATES 
  0.0     0.0       3          3 
# i  j  S 
  1  2  0.1 
  2  3  0.5 
# EIR 
  0   20   40 
END 
//! [EDS]

//! [LAMBDAS]
LAMBDAS
# NTIL    off(0), on(1) 
#         0: no special treatment of interactions with individual lambda-values 
#         1: interactions are treated with special individual lambda-values 
# NTLI(1..)  interaction type to treat with individual lambda: 
#            bond(1), angle(2), dihedral(3), improper(4), vdw(5), vdw_soft(6), 
#            crf(7), crf_soft(8), distanceres(9), distancefield(10), 
#            dihedralres(11), mass(12) 
# NILG1, NILG2 energy groups of interactions that are treated with individual 
#              lambda values 
# ALI, BLI, CLI, DLI, ELI polynomial coefficients linking the individual lambda- 
#                         values to the overall lambda-value 
# NTIL 
   1 
# NTLI NILG1  NILG2  ALI   BLI   CLI   DLI   ELI 
    7      1      3    0     0     1     0     0 
END 
//! [LAMBDAS]

//! [PRECALCLAM]
PRECALCLAM
# NRLAM   0  : off 
#         >1 : precalculating energies for NRLAM extra lambda values 
# MINLAM  between 0 and 1: minimum lambda value to precalculate energies 
# MAXLAM  between MINLAM and 1: maximum lambda value to precalculate energies 
# NRLAM	  MINLAM   MAXLAM 
   100      0.0        1.0 
END 
//! [PRECALCLAM]

//! [NONBONDED]
NONBONDED
# NLRELE    1-3 method to handle electrostatic interactions 
#    -1 : reaction-field (LSERF compatibility mode) 
#     0 : no electrostatic interactions 
#     1 : reaction-field 
#     2 : Ewald method 
#     3 : P3M method 
# APPAK     >= 0.0 reaction-field inverse Debye screening length 
# RCRF      >= 0.0 reaction-field radius 
#   0.0 : set to infinity 
# EPSRF     = 0.0 || > 1.0 reaction-field permittivity 
#   0.0 : set to infinity 
# NSLFEXCL  0,1 contribution of excluded atoms to reaction field 
#     0 : contribution turned off 
#     1 : contribution considered (default) 
# NSHAPE    -1..10 lattice sum charge-shaping function 
#    -1 : gaussian 
# 0..10 : polynomial 
# ASHAPE    > 0.0 width of the lattice sum charge-shaping function 
# NA2CALC   0..4 controls evaluation of lattice sum A2 term 
#     0 : A2 = A2~ = 0 
#     1 : A2~ exact, A2 = A2~ 
#     2 : A2 numerical, A2~ = A2 
#     3 : A2~ exact from Ewald or from mesh and atom coords, A2 numerical 
#     4 : A2~ averaged from mesh only, A2 numerical 
# TOLA2     > 0.0 tolerance for numerical A2 evaluation 
# EPSLS      = 0.0 || > 1.0 lattice sum permittivity (0.0 = tinfoil) 
# NKX, NKY, NKZ > 0 maximum absolute Ewald k-vector components 
# KCUT       > 0.0 Ewald k-space cutoff 
# NGX, NGY, NGZ > 0 P3M number of grid points 
# NASORD    1..5 order of mesh charge assignment function 
# NFDORD    0..5 order of the mesh finite difference operator 
#     0 : ik - differentiation 
#  1..5 : finite differentiation 
# NALIAS    > 0 number of mesh alias vectors considered 
# NSPORD        order of SPME B-spline functions (not available) 
# NQEVAL    >= 0 controls accuracy reevaluation 
#     0 : do not reevaluate 
#   > 0 : evaluate every NQEVAL steps 
# FACCUR    > 0.0 rms force error threshold to recompute influence function 
# NRDGRD    0,1 read influence function 
#     0 : calculate influence function at simulation start up 
#     1 : read influence function from file (not yet implemented) 
# NWRGRD    0,1 write influence function 
#     0 : do not write 
#     1 : write at the end of the simulation (not yet implemented) 
# NLRLJ     0,1 controls long-range Lennard-Jones corrections 
#     0 : no corrections 
#     1 : do corrections (not yet implemented) 
# SLVDNS    > 0.0 average solvent density for long-range LJ correction (ignored) 
# 
#   NLRELE 
         1 
#    APPAK      RCRF     EPSRF   NSLFEXCL 
       0.0       1.4      61.0          1 
#   NSHAPE    ASHAPE    NA2CLC      TOLA2    EPSLS 
        -1       1.4         2     0.1E-9      0.0 
#      NKX       NKY       NKZ       KCUT 
        10        10        10      100.0 
#      NGX       NGY       NGZ    NASORD    NFDORD    NALIAS    NSPORD 
        32        32        32         3         2         3         4 
#   NQEVAL    FACCUR    NRDGRD    NWRGDR 
    100000       1.6         0         0 
#    NLRLJ    SLVDNS 
         0      33.3 
END 
//! [NONBONDED]

//! [SASA]
SASA
# NTSASA 
# 0 : not used (default) 
# 1 : use SASA implicit solvent model 
# NTVOL 
# 0 : not used (default) 
# 1 : use VOLUME correction to SASA implicit solvent model (requires NTSASA = 1) 
# P_12 >= 0, <= 1 pair parameter for SASA reduction for first neighbours 
# P_13 >= 0, <= 1 pair parameter for SASA reduction for second neighbours 
# P_1X >= 0, <= 1 pair parameter for SASA reduction for third and higher neighbours 
# SIGMAV >0 scaling parameter for volume energy term (kJ.mol^-1.nm^-3) 
# RSOLV > 0 radius of solvent molecule for SASA calculation (nm) 
# AS1 > 0 an atom with SASA below this contributes to the VOLUME correction (nm^2) 
# AS2 > 0 an atom with SASA above this is not considered for the VOLUME correction (nm^2) 
# atoms with AS1 < SASA < AS2 have a partial contribution determined by a switching function 
#   NTSASA      NTVOL       P_12      P_13     P_1X   SIGMAV  RSOlV    AS1    AS2 
         1          1     0.8875    0.3516   0.3516     -100   0.14   0.01   0.02 
END 
//! [SASA]

//! [LOCALELEV]
LOCALELEV
# NTLES 0,1 controls the use of local elevation. 
#    0 : not used [default] 
#    1 : local elevation is applied 
# NLEPOT >= 0 number of umbrella potentials applied 
# NTLESA 1..2 controls the reading of the potential definition 
#    1 : read from startup file 
#    2 : read from special file (@lud) 
# NTWLE >= 0 write umbrellas to trajectory every NTWLEth step 
# NLEPID[1..NLEPOT] IDs of the umbrella potentials 
# NTLEPFR[1..NLEPOT] 0,1 freeze the umbrella potential 
#    0 : build up 
#    1 : freeze 
# NTLES  NLEPOT  NTLESA  NTWLE 
      1       2       1      0 
# NLEPID NLEPFR 
       1      0 
       2      1 
END 
//! [LOCALELEV]

//! [BSLEUS]
BSLEUS
# 
# The general settings for the B&S-LEUS algorithm 
# BSLEUS:   Dow we use B&S-LEUS? 
#   0:          Don'use it (default) 
#   1:          Use it 
# BUILD:    Are we building? 
#   0:          No 
#   1:          Yes 
# WRITE:    >= 0 Do we write the energies and forces of the Umbrella? 
#   == 0:          No 
#   > 0:           Every nth step 
# 
# BSLEUS    BUILD   WRITE 
  1         1       0 
END 
//! [BSLEUS]

//! [ELECTRIC]
ELECTRIC
# FIELD 0,1 controls the use of applied electric field. 
#    0 : not used [default] 
#    1 : electric field is applied 
# DIPOLE 0,1 controls the calculation of box dipole. 
#    0 : not used [default] 
#    1 : box dipole is calculated and written to special trajectory 
# CURRENT 0,1 controls the calculation of electric currents. 
#    0 : not used [default] 
#    1 : electric current is calculated and written to special trajectory 
# ELECTRIC FIELD COMPONENTS (EF_x, EF_y, EF_z) 
# 0.0 0.0 0.0 
# DIPGRP 0..2 controls the groups considered for box dipole calculation 
#    0 : solute only 
#    1 : solvent only 
#    2 : all 
# NTWDIP >= 0 write dipole box every NTWDIPth step 
# NTWCUR >= 0 write current every NTWDIPth step 
# NCURGRP >=0 number of current groups 
# CURGRP [1..NCURGRP] last atom of the group 
#  FIELD  DIPOLE CURRENT 
       1       1       1 
#   EF_x    EF_y    EF_z 
     0.0     0.0     0.0 
# DIPGRP  NTWDIP 
       0       1 
# NTWCUR  NCURGRP   CURGRP[1]   CURGRP[2] 
       1       2        100        1000 
END 
//! [ELECTRIC]

//! [NEMD]
NEMD
# NEMD 0,1 controls the use of non-equilibrium molecular dynamics. 
#    0 : not used [default] 
#    1 : nemd is used 
# PROPERTY 0- select property to calculate 
#    0 : viscosity 
# METHOD 0- select method of NEMD. 
#    0 : periodic perturbation method (PPM) 
#    1 : internal reservoir method (IRM) 
# SLABNUM >=1 number of slabs used in the discretization along z-direction. 
#             the effective number is 2xSLABNUM due to periodicity 
# PERTFRQ >=1 perturbation frequency: apply perturbation every PERTFRQth timestep 
#             [this flag is ignored by the PPM method, but a value must be provided] 
# AMPLI   >=0 amplitude of applied field 
#             [this flag is ignored by the IRM method, but a value must be provided] 
# STDYAFT >=0 first STDYAFTth steps do not contribute for accumulated averages 
# WRITE >=1 write flux and average velocities to special trajectory every WRITEth timestep 
# NEMD     PROPERTY  METHOD 
    1         0        0 
# SLABNUM  PERTFRQ    AMPLI   STDYAFT   WRITE 
     10       20       10      1000     200 
END 
//! [NEMD]

//! [MULTIGRADIENT]
MULTIGRADIENT
# NTMGRE 0,1 enable multiple gradients 
#    0: disable gradients (default) 
#    1: enable gradients 
# NTMGRP 0..3 print of curves 
#    0: don't print 
#    1: plot the curves 
#    2: print that values of the curves 
#    3: plot and print the curves 
# NTMGRN >= 0 number of gradients 
# MGRVAR: variable name to affect, available are: 
    TEMP0, CPOR, CDIR, RESO, CXR, COPR 
# MGRFRM: functional form of the curve 
#    0: linear interpolation between control points 
#    1: cubic spline interpolation between control points 
#    2: Bezier curve 
#    3: Oscillation: A sin[2Pi/T (t - dt)] + b 
#       Note: MGRNCP is 2. A = MGRCPT[1] T = MGRCPV[1] dt = MGRCPT[2] b = MGRCPV[2] 
# MGRNCP >= 2: number of control points 
# MGRCPT >= 0: time of the control point 
# MGRCPV: value of the control point 
# 
# NTMGRE NTMGRP 
       1      1 
# NTMGRN 
       2 
# MGRVAR MGRFRM MGRNCP 
  TEMP0[0]     0      2 
# MGRCPT MGRCPV 
  0.0    60.0 
  80.0   300.0 
# MGRVAR MGRFRM MGRNCP 
  CPOR        2      4 
# MGRCPT MGRCPV 
  0.0    2.5E5 
  0.0    2.5E1 
 20.0    0.0 
 80.0    0.0 
END 
//! [MULTIGRADIENT]

//! [ADDECOUPLE]
ADDECOUPLE
# ADGR    >= 0 number of adiabatic decoupling groups 
# ADSTART first atom of the adiabatic decoupling group 
# ADEND   last atom of the adiabatic decoupling group 
# SM      scaling factor mass 
# SV      scaling factor potential energy function 
# ST      scaling factor temperature 
# TIR     which temperature bath to scale 
#  1      translational 
#  2      internal-rotatinal 
#  3      both 
# TMF     tau for calculating mean field 
# STAD    printing average to special trajectory 
# ADGR 
      2 
# ADSTART ADEND SM SV  ST TIR 
  1       1500 10   1  0  1 
  1501    3000  1  10  1  3 
# TMF STAD 
  0.1 1000 
END 
//! [ADDECOUPLE]

//! [QMMM]
QMMM
# NTQMMM 0,1 apply QM/MM 
#    0: do not apply QM/MM (default) 
#    1: apply QM/MM 
# NTQMSW 0 QM software package to use 
#    0: MNDO (default) 
#    1: Turbomole 
# RCUTQ >= 0.0 cutoff for inclusion of MM charge groups 
#     0.0: include all atoms 
#    >0.0: include atoms of charge groups closer than RCUTQ 
#          to QM zone. 
# NTWQMMM >= 0 write QM/MM related data to special trajectory 
#    0: do not write 
#   >0: write every NTWQMMMth step 
# 
# NTQMMM  NTQMSW  RCUTQ  NTWQMMM 
       1       0    0.0        0 
END 
//! [QMMM]

//! [SYMRES]
SYMRES
# NTSYM 0..2 apply symmetry restraints 
#    0: do not apply symmetry restraints (default) 
#    1: apply symmetry restraints 
#    2: apply symmetry constraints 
# CSYM >= 0.0 force constants 
# 
# NTSYM     CSYM 
       1     0.0 
END 
//! [SYMRES]

//! [DISTANCERESSPEC]
DISTANCERESSPEC
# DISH, DISC carbon-hydrogen/carbon-carbon distance 
# i,j,k,l  atoms comprising the virtual atom (put 0 if less than four atoms in use) 
# type   virtual atom type 
# r0, w0  target distance and force constant weighting factor 
# rah    form and dimension of the potential 
# full harmonic: 
#     0: x,y,z 
#    10: x,y 
#    20: x,z 
#    30: y,z 
#    40: x 
#    50: y 
#    60: z 
#  subtract or add 1 from these numbers to select a half harmonic 
#  repulsive or attractive potential 
# DISH  DISC 
  0.1   0.153 
# i  j  k  l  type    i  j  k  l  type    r0    w0    rah 
  1  0  0  0  0       10 12 11 13 3       0.2   1.0   0 
END 
//! [DISTANCERESSPEC]

//! [PERTDISRESSPEC]
PERTDISRESSPEC
# DISH  DISC 
  0.1   0.153 
# i  j  k  l  type    i  j  k  l  type n m   A_r0  A_w0  B_r0   B_w0  rah 
  1  0  0  0  0       10 12 11 13 3    1 1    0.2   1.0   0.5    2.0   0 
END 
//! [PERTDISRESSPEC]

//! [DFRESSPEC]
DFRESSPEC
#   DISH H-C bond length for virtual atoms 
#   DISC C-C bond length for virtual atoms 
#   PROTEINATOMS > 0 last atom of the host 
#   K >= 0.0 Force constant 
#   r0 >=0 zero energy distance 
#   TYPE_I Virtual atom type for interaction site I 
#   NUM_I  Number of atoms defining interaction site I 
#   ATOM_I[0..NUM_I] Index numbers of atoms defining interaction site I 
#   TYPE_J Virtual atom type for interaction site J 
#   NUM_J  Number of atoms defining interaction site J 
#   ATOM_J[0..NUM_J] Index numbers of atoms defining interaction site J 
# DISH  DISC 
  0.1   0.153 
# PROTEINATOMS  K    r0 
  1190          500  0.0 
# TYPE_I  NUM_I  ATOM_I[0] .. ATOM_I[NUM_I] 
  -1      7        16  190  249  312  486  632 1208 
# TYPE_J  NUM_J  ATOM_J[0] .. ATOM_J[NUM_J] 
  -1      2      1194 1203 
END 
//! [DFRESSPEC]

//! [PERTDFRESSPEC]
PERTDFRESSPEC
#   DISH H-C bond length for virtual atoms 
#   DISC C-C bond length for virtual atoms 
#   PROTEINATOMS > 0 last atom of the host 
#   A_r0 >=0 reference distance for state A 
#   B_r0 >=0 reference distance for state B 
#   K_A >= 0 force constant state A 
#   K_B >= 0 force constant state B 
#   n >= 0 hidden restraint parameter n 
#   m >= 0 hidden restraint parameter m 
#   TYPE_I Virtual atom type for interaction site I 
#   NUM_I  Number of atoms defining interaction site I 
#   ATOM_I[0..NUM_I] Index numbers of atoms defining interaction site I 
#   TYPE_J Virtual atom type for interaction site J 
#   NUM_J  Number of atoms defining interaction site J 
#   ATOM_J[0..NUM_J] Index numbers of atoms defining interaction site J 
# DISH  DISC 
  0.1   0.153 
# PROTEINATOMS  A_r0  K_A  B_r0  K_B  n  m 
  1190          4.5   500  0.0   500  0  0 
# TYPE_I  NUM_I  ATOM_I[0] .. ATOM_I[NUM_I] 
  -1      7        16  190  249  312  486  632 1208 
# TYPE_J  NUM_J  ATOM_J[0] .. ATOM_J[NUM_J] 
  -1      2      1194 1203 
END 
//! [PERTDFRESSPEC]

//! [MDISRESSPEC]
MDISRESSPEC
# DISH  DISC 
  0.1   0.153 
# N: number of eds states (3 in this example) 
# i  j  k  l  type    i  j  k  l  type    r0[1 ... N]    w0[1 ... N]    rah 
  1  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  1.0  0.0 0.0    0 
  5  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  0.0  1.0 0.0    0 
  8  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  0.0  0.0 1.0    0 
END 
//! [MDISRESSPEC]

