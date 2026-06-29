##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Wed Jun 24 15:26:13 2026

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path holoWT_prot.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 22589 atoms
Valist_getStatistics:  Max atom coordinate:  (144.245, 141.555, 137.061)
Valist_getStatistics:  Min atom coordinate:  (6.67, 16.552, 28.194)
Valist_getStatistics:  Molecule center:  (75.4575, 79.0535, 82.6275)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 75.4575 79.0535 82.6275
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 75.4575 79.0535 82.6275
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 2.64479, 2.64479, 2.64479
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 1.82396, 1.82396, 1.82396
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.689642, 0.689642, 0.689642 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 75.4575 79.0535 82.6275
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 202.407 206.004 209.577
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -51.4925 -47.8965 -44.3225
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 163.007 166.603 170.178
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -12.0925 -8.4965 -4.9225
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 163.007 166.603 170.178
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -12.0925 -8.4965 -4.9225
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 76.7418
Vpbe_ctor2:  solute dimensions = 140.491 x 127.861 x 111.741
Vpbe_ctor2:  solute charge = -9
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (150.355, 137.783, 121.647)
Vclist_setupGrid:  Grid lower corner = (0.28, 10.162, 21.804)
Vclist_assignAtoms:  Have 6056915 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.416176
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.121289e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.403200e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.367240e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 5.221556e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 5.218116e-02
Vprtstp: contraction number = 5.218116e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 6.582254e-03
Vprtstp: contraction number = 1.261424e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.030528e-03
Vprtstp: contraction number = 1.565615e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.982451e-04
Vprtstp: contraction number = 1.923724e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 4.608150e-05
Vprtstp: contraction number = 2.324471e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.259247e-05
Vprtstp: contraction number = 2.732653e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 3.906941e-06
Vprtstp: contraction number = 3.102600e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 1.331601e-06
Vprtstp: contraction number = 3.408296e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 4.840394e-07
Vprtstp: contraction number = 3.635019e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.664657e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.770336e+01
Vpmg_setPart:  lower corner = (-51.4925, -47.8965, -44.3225)
Vpmg_setPart:  upper corner = (202.407, 206.004, 209.577)
Vpmg_setPart:  actual minima = (-51.4925, -47.8965, -44.3225)
Vpmg_setPart:  actual maxima = (202.408, 206.004, 209.578)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 76.7418
Vpbe_ctor2:  solute dimensions = 140.491 x 127.861 x 111.741
Vpbe_ctor2:  solute charge = -9
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (150.355, 137.783, 121.647)
Vclist_setupGrid:  Grid lower corner = (0.28, 10.162, 21.804)
Vclist_assignAtoms:  Have 6056915 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -12.0925, -8.4965, -4.9225
VPMG::focusFillBound -- New mesh maxs = 163.007, 166.603, 170.178
VPMG::focusFillBound -- Old mesh mins = -51.4925, -47.8965, -44.3225
VPMG::focusFillBound -- Old mesh maxs = 202.408, 206.004, 209.578
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.399846
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.183704e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.447200e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 5.870080e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.709062e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 4.352815e-02
Vprtstp: contraction number = 4.352815e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 5.559909e-03
Vprtstp: contraction number = 1.277313e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 8.920647e-04
Vprtstp: contraction number = 1.604459e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.610080e-04
Vprtstp: contraction number = 1.804891e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.263482e-05
Vprtstp: contraction number = 2.026907e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 7.324779e-06
Vprtstp: contraction number = 2.244468e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.857227e-06
Vprtstp: contraction number = 2.535539e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 5.296980e-07
Vprtstp: contraction number = 2.852092e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 6.526179e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 7.534590e+00
Vpmg_setPart:  lower corner = (-12.0925, -8.4965, -4.9225)
Vpmg_setPart:  upper corner = (163.007, 166.603, 170.178)
Vpmg_setPart:  actual minima = (-12.0925, -8.4965, -4.9225)
Vpmg_setPart:  actual maxima = (163.007, 166.603, 170.178)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 3.392896e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Wed Jun 24 15:26:28 2026

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path holoMutant_prot.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 22609 atoms
Valist_getStatistics:  Max atom coordinate:  (145.544, 140.017, 137.501)
Valist_getStatistics:  Min atom coordinate:  (-11.332, 11.084, 24.91)
Valist_getStatistics:  Molecule center:  (67.106, 75.5505, 81.2055)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 67.106 75.5505 81.2055
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 67.106 75.5505 81.2055
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 2.98646, 2.98646, 2.98646
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 2.06563, 2.06563, 2.06563
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.691664, 0.691664, 0.691664 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 67.106 75.5505 81.2055
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 210.456 218.9 224.555
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -76.244 -67.7995 -62.1445
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 166.256 174.701 180.356
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -32.044 -23.5995 -17.9445
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 166.256 174.701 180.356
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -32.044 -23.5995 -17.9445
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 85.1127
Vpbe_ctor2:  solute dimensions = 159.556 x 131.274 x 115.639
Vpbe_ctor2:  solute charge = -8
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (169.656, 141.713, 125.371)
Vclist_setupGrid:  Grid lower corner = (-17.722, 4.694, 18.52)
Vclist_assignAtoms:  Have 5379118 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.566909
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.245959e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.416100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.451810e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 5.260935e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 5.576246e-02
Vprtstp: contraction number = 5.576246e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 6.716345e-03
Vprtstp: contraction number = 1.204456e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 9.674565e-04
Vprtstp: contraction number = 1.440451e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.553613e-04
Vprtstp: contraction number = 1.605874e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.695684e-05
Vprtstp: contraction number = 1.735107e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.991514e-06
Vprtstp: contraction number = 1.851669e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 9.802356e-07
Vprtstp: contraction number = 1.963804e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.309664e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.407220e+01
Vpmg_setPart:  lower corner = (-76.244, -67.7995, -62.1445)
Vpmg_setPart:  upper corner = (210.456, 218.9, 224.555)
Vpmg_setPart:  actual minima = (-76.244, -67.7995, -62.1445)
Vpmg_setPart:  actual maxima = (210.456, 218.9, 224.555)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 85.1127
Vpbe_ctor2:  solute dimensions = 159.556 x 131.274 x 115.639
Vpbe_ctor2:  solute charge = -8
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (169.656, 141.713, 125.371)
Vclist_setupGrid:  Grid lower corner = (-17.722, 4.694, 18.52)
Vclist_assignAtoms:  Have 5379118 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -32.044, -23.5995, -17.9445
VPMG::focusFillBound -- New mesh maxs = 166.256, 174.701, 180.356
VPMG::focusFillBound -- Old mesh mins = -76.244, -67.7995, -62.1445
VPMG::focusFillBound -- Old mesh maxs = 210.456, 218.9, 224.555
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.561451
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.337659e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.415000e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.372900e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.381148e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 4.711019e-02
Vprtstp: contraction number = 4.711019e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 6.088169e-03
Vprtstp: contraction number = 1.292325e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 9.152032e-04
Vprtstp: contraction number = 1.503249e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.588897e-04
Vprtstp: contraction number = 1.736114e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.022915e-05
Vprtstp: contraction number = 1.902525e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 6.459623e-06
Vprtstp: contraction number = 2.136885e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.501779e-06
Vprtstp: contraction number = 2.324871e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 3.784774e-07
Vprtstp: contraction number = 2.520193e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.496202e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.603454e+01
Vpmg_setPart:  lower corner = (-32.044, -23.5995, -17.9445)
Vpmg_setPart:  upper corner = (166.256, 174.701, 180.356)
Vpmg_setPart:  actual minima = (-32.044, -23.5995, -17.9445)
Vpmg_setPart:  actual maxima = (166.256, 174.701, 180.356)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 3.909001e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Wed Jun 24 15:26:43 2026

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path apoMutant.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 22598 atoms
Valist_getStatistics:  Max atom coordinate:  (141.171, 136.281, 137.429)
Valist_getStatistics:  Min atom coordinate:  (2.876, 18.83, 28.597)
Valist_getStatistics:  Molecule center:  (72.0235, 77.5555, 83.013)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 72.0235 77.5555 83.013
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 72.0235 77.5555 83.013
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 2.65729, 2.65729, 2.65729
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 1.83333, 1.83333, 1.83333
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.689926, 0.689926, 0.689926 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 72.0235 77.5555 83.013
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 199.573 205.106 210.563
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -55.5265 -49.9945 -44.537
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 160.024 165.555 171.013
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -15.9765 -10.4445 -4.987
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 160.024 165.555 171.013
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -15.9765 -10.4445 -4.987
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 77.5929
Vpbe_ctor2:  solute dimensions = 141.008 x 119.832 x 111.958
Vpbe_ctor2:  solute charge = -7
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (151.075, 130.231, 121.612)
Vclist_setupGrid:  Grid lower corner = (-3.514, 12.44, 22.207)
Vclist_assignAtoms:  Have 6278836 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.424703
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.115007e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.400500e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.363320e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 5.209069e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 5.173899e-02
Vprtstp: contraction number = 5.173899e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 6.447507e-03
Vprtstp: contraction number = 1.246160e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.021320e-03
Vprtstp: contraction number = 1.584055e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.019187e-04
Vprtstp: contraction number = 1.977035e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 4.654018e-05
Vprtstp: contraction number = 2.304898e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.197333e-05
Vprtstp: contraction number = 2.572686e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 3.296308e-06
Vprtstp: contraction number = 2.753043e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 9.477425e-07
Vprtstp: contraction number = 2.875164e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.484034e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.589466e+01
Vpmg_setPart:  lower corner = (-55.5265, -49.9945, -44.537)
Vpmg_setPart:  upper corner = (199.573, 205.106, 210.563)
Vpmg_setPart:  actual minima = (-55.5265, -49.9945, -44.537)
Vpmg_setPart:  actual maxima = (199.573, 205.106, 210.563)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 77.5929
Vpbe_ctor2:  solute dimensions = 141.008 x 119.832 x 111.958
Vpbe_ctor2:  solute charge = -7
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (151.075, 130.231, 121.612)
Vclist_setupGrid:  Grid lower corner = (-3.514, 12.44, 22.207)
Vclist_assignAtoms:  Have 6278836 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -15.9765, -10.4445, -4.987
VPMG::focusFillBound -- New mesh maxs = 160.024, 165.555, 171.013
VPMG::focusFillBound -- Old mesh mins = -55.5265, -49.9945, -44.537
VPMG::focusFillBound -- Old mesh maxs = 199.573, 205.106, 210.563
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.457660
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.239323e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.419000e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.363290e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.538947e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 4.295402e-02
Vprtstp: contraction number = 4.295402e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 5.686139e-03
Vprtstp: contraction number = 1.323773e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 9.536311e-04
Vprtstp: contraction number = 1.677115e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.837898e-04
Vprtstp: contraction number = 1.927263e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 4.003063e-05
Vprtstp: contraction number = 2.178066e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 9.306162e-06
Vprtstp: contraction number = 2.324760e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 2.376891e-06
Vprtstp: contraction number = 2.554105e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 6.362140e-07
Vprtstp: contraction number = 2.676664e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.120465e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.226236e+01
Vpmg_setPart:  lower corner = (-15.9765, -10.4445, -4.987)
Vpmg_setPart:  upper corner = (160.024, 165.555, 171.013)
Vpmg_setPart:  actual minima = (-15.9765, -10.4445, -4.987)
Vpmg_setPart:  actual maxima = (160.024, 165.555, 171.013)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 3.691391e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Wed Jun 24 15:26:58 2026

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path apoWT.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 22591 atoms
Valist_getStatistics:  Max atom coordinate:  (143.071, 139.391, 137.499)
Valist_getStatistics:  Min atom coordinate:  (2.745, 12.558, 27.043)
Valist_getStatistics:  Molecule center:  (72.908, 75.9745, 82.271)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 72.908 75.9745 82.271
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 72.908 75.9745 82.271
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 2.69375, 2.69375, 2.69375
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 1.85833, 1.85833, 1.85833
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.689869, 0.689869, 0.689869 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 72.908 75.9745 82.271
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 202.208 205.274 211.571
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -56.392 -53.3255 -47.029
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 162.108 165.174 171.471
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -16.292 -13.2255 -6.929
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 162.108 165.174 171.471
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -16.292 -13.2255 -6.929
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 76.7432
Vpbe_ctor2:  solute dimensions = 142.942 x 129.433 x 113.144
Vpbe_ctor2:  solute charge = -7
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (153.106, 139.613, 123.236)
Vclist_setupGrid:  Grid lower corner = (-3.645, 6.168, 20.653)
Vclist_assignAtoms:  Have 5882806 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.480994
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.071954e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.395200e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.342910e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 5.168621e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 5.218905e-02
Vprtstp: contraction number = 5.218905e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 6.431116e-03
Vprtstp: contraction number = 1.232273e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 9.365493e-04
Vprtstp: contraction number = 1.456278e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.533762e-04
Vprtstp: contraction number = 1.637674e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.672097e-05
Vprtstp: contraction number = 1.742185e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.929154e-06
Vprtstp: contraction number = 1.844676e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 9.380183e-07
Vprtstp: contraction number = 1.903001e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.217543e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.322862e+01
Vpmg_setPart:  lower corner = (-56.392, -53.3255, -47.029)
Vpmg_setPart:  upper corner = (202.208, 205.274, 211.571)
Vpmg_setPart:  actual minima = (-56.392, -53.3255, -47.029)
Vpmg_setPart:  actual maxima = (202.208, 205.274, 211.571)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 0.000000e+00
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 76.7432
Vpbe_ctor2:  solute dimensions = 142.942 x 129.433 x 113.144
Vpbe_ctor2:  solute charge = -7
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.124825
Vpbe_ctor2:  Debye length = 8.01121
Vpbe_ctor2:  zkappa2 = 1.22376
Vpbe_ctor2:  zmagic = 6773.76
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (153.106, 139.613, 123.236)
Vclist_setupGrid:  Grid lower corner = (-3.645, 6.168, 20.653)
Vclist_assignAtoms:  Have 5882806 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -16.292, -13.2255, -6.929
VPMG::focusFillBound -- New mesh maxs = 162.108, 165.174, 171.471
VPMG::focusFillBound -- Old mesh mins = -56.392, -53.3255, -47.029
VPMG::focusFillBound -- Old mesh maxs = 202.208, 205.274, 211.571
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 3.455697
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 4.190336e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.443100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 7.236330e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.283077e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 4.380989e-02
Vprtstp: contraction number = 4.380989e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 5.798539e-03
Vprtstp: contraction number = 1.323569e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 9.306953e-04
Vprtstp: contraction number = 1.605051e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.675840e-04
Vprtstp: contraction number = 1.800632e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.319247e-05
Vprtstp: contraction number = 1.980647e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 7.115465e-06
Vprtstp: contraction number = 2.143698e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.712633e-06
Vprtstp: contraction number = 2.406916e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 4.626457e-07
Vprtstp: contraction number = 2.701371e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 8.679065e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 9.943738e+00
Vpmg_setPart:  lower corner = (-16.292, -13.2255, -6.929)
Vpmg_setPart:  upper corner = (162.108, 165.174, 171.471)
Vpmg_setPart:  actual minima = (-16.292, -13.2255, -6.929)
Vpmg_setPart:  actual maxima = (162.108, 165.174, 171.471)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 3.182626e+01
