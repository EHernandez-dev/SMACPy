from time import sleep
import subprocess

''''''

def create_smac_par( path = './', fileName = 'smac_1.par',**kwargs):

    default=dict(eraseFits =  False,
                SIM_FORMAT = 4, 
                USE_KEYS = 1,
                SNAP_BASE = 'snap',
                SNAP_FILE = 'snap_',
                HALO_ID = 'a', 
                SNAP_START = 36,
                SNAP_END   = 2, 
                INITIAL_Z = 0.0, 
                IMAGE_Z   = 0.0,
                R_VIR = 2300.,
                M_VIR = 1.0e15, 
                FILE_FORMAT = 3, 
                FILE_HEADER = 1, 
                KERNEL_TYPE = 0,  
                PART_DISTR = 1, 
                PART_EXCL = 1,
                IMG_XY_UNITS = 2,
                IMG_XY_SIZE = 8000.0,
                IMG_Z_UNITS = 2,
                IMG_Z_SIZE = 8000.0,
                FILL_CARROT = 0.,
                IMG_SIZE = 1024,
                SMOOTH_UNITS = 0,
                SMOOTH_FWHM  = 5.,
                NSIDE = 8,
                MIN_DIST =   5000.,
                MAX_DIST = 420000.,
                OUTPUT_MAP = 1,
                OUTPUT_SUB = 1,
                METAL_SPECIES = 0,
                TEMP_UNIT = 1,
                TEMP_CUT = 0.0,
                XRAY_E0 = 0.5,
                XRAY_E1 = 2.,
                X_TABLES = './Tables',
                FRQ_ZS = 300.,
                ECRP_MIN = 1e9,
                FRQ_Pnu = 1.4e9,
                XCRP = 0.01,
                XCRP_TABLE_PATH = '~/',
                B_TABLE_PATH = '~/',
                KERNEL_TABLE_PATH = '',
                GAM_nu = 1.25,
                XCRP_TABLE_FILE_FMT = '(1F8.6,2X,1F14.6)',
                E_gam = 1e1,
                E_gam_min = 1e-2,
                E_gam_max = 30e3,
                CR_nbins = 48,
                CR_pmin = 10.0,
                CR_pmax = 1.e7,
                CR_subsamples = 10,
                CR_DSlope = 1.0e-6,
                GIVE_MORE_INFO = 0,
                PROJECT = 4,
                CENTER_MOTION = 0,
                REMOVE_LOCAL_GROUP_VEL = 0,
                CENTER = 1,
                CENTER_X = 0.,
                CENTER_Y = 0.,
                CENTER_Z = 0.,
                PERIODIC = 1,
                LIGHTCONE = 0,
                X_ORIGIN = 0.,
                Y_ORIGIN = 0.,
                Z_ORIGIN = 0.,
                OPEN_ANGLE = 1.0,
                MAIN_PROG = '/afs/mpa/project/hydrosims/Hutt/g676/csf/Post/main_prog.a_gas.gv',               
                PREFIX_OUT = 'nonPrefix',
                HUBBLE = 0.6777,
                OMEGA  = 0.307115,
                LAMBDA = 0.692885,
                )

# add default values for non modified keys to write them in the parameter file
    for key, value in default.items():
        if key not in kwargs:
            kwargs[key] = value

    f = open(path + fileName, "w")
    f.write(
    f"""#=============================================================
    # Input file
    #=============================================================
    
    #**** Flag of the data format (integer - reals):
    #         (1) Navarro,
    #         (2) GADGET
    #             intial redshift, (case (1))
    #             omega matter,    (case (1))
    #             Hubble parameter (case (1))
    #             image redshift (set to 0. if obtained from the input data)
    #        (3) TOY-MODELS
    #        (4) GADGET2
    SIM_FORMAT = {kwargs['SIM_FORMAT']}

    #**** Input file
    USE_KEYS = {kwargs['USE_KEYS']}
    OUTPUT_DIR = {kwargs['OUTPUT_DIR']} 
    #SNAP_BASE = {kwargs['SNAP_BASE']}
    #SNAP_FILE = {kwargs['SNAP_FILE']}


    #**** Flag of the loop on the input files (integers):
    #        disabled if SNAP_START = -1 else enabled
    #        first input
    #        last input  (=first_input if not given)
    SNAP_START =  {kwargs['SNAP_START']}
    #SNAP_END  =  {kwargs['SNAP_END']}

    #**** Halo ID
    HALO_ID = {kwargs['HALO_ID']}   
    #**** Intial redshift and image redshift (set IMAGE_Z=0. to obtain it from the input data)
    INITIAL_Z = {kwargs['INITIAL_Z']}
    IMAGE_Z   = {kwargs['IMAGE_Z']}


    #=============================================================
    # KERNEL TO USE
    #=============================================================
    # Flag of kernel used for mapping (integer):
    #         (0) cubic-spline
    #         (1) quintic-spline
    #         (2) WendlandC4
    #         (3) WendlandC6
    KERNEL_TYPE =  {kwargs['KERNEL_TYPE']}

    #=============================================================
    # DISTRIBUTION SCHEME
    #=============================================================

    #**** Flag of particle distribution scheme on the image (integer):
    #         (0) NGC       ->  not implemented
    #         (1) CIC       ->  2D-map
    #         (2) HealPix   ->  full-sky
    #         (3) TSC       ->  2D-map
    PART_DISTR = {kwargs['PART_DISTR']} 
    #**** Flag to exclude particles (integer, default=1):
    #         (0) All the particles are included
    #         (1) Particles with rho > 500 rho_crit(z) * Omega_b AND T < 3e4K are excluded (see Borgani et al. 2002 for details)
    PART_EXCL = {kwargs['PART_EXCL']}

    #=============================================================
    # 2D-MAP SETUP
    #=============================================================

    #**** Flag of the units of the side of image (integer-real):
    #         (1) arcmin
    #         (2) kpc
    #         side of image in the chosen units
    IMG_XY_UNITS = {kwargs['IMG_XY_UNITS']}
    IMG_XY_SIZE  = {kwargs['IMG_XY_SIZE']}

    #**** Flag of the third dimension of the image:
    #        (1) cubic,  same value of the preciding point
    #	 (2) if the value is chosen by the man
    #         side of the third dimension in kpc (case(2))
    #	  minimum filling factor of the carrot [0<->100] (its a percentile)
    IMG_Z_UNITS = {kwargs['IMG_Z_UNITS']}
    IMG_Z_SIZE  = {kwargs['IMG_Z_SIZE']}
    FILL_CARROT = {kwargs['FILL_CARROT']}

    #**** Flag of the number of image side pixels (integers):
    #         nr of pixel
    #         if (=-1) -> given by smoothing
    IMG_SIZE = {kwargs['IMG_SIZE']}

    #**** Flag of the smoothing (FWHM) angle (integer - real):
    #         (0) no smoothing
    #         (1) in arcmin
    #         (2) in kpc
    #         (3) given by npart
    #       smoothing lenght in the choosen units (cases (1), (2))
    SMOOTH_UNITS = {kwargs['SMOOTH_UNITS']}
    SMOOTH_FWHM  = {kwargs['SMOOTH_FWHM']}

    #=============================================================
    # FULL-SKY MAP SETUP
    #=============================================================

    #**** nside HEALPix variable (Mauro)
    #           This variable must be a power of 2 and <= 8192
    #           nside = ( 1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192)
    NSIDE = {kwargs['NSIDE']}

    #**** min/max radius of sphere arround x0,y0,z0 to be used [kpc]
    MIN_DIST = {kwargs['MIN_DIST']}
    MAX_DIST = {kwargs['MAX_DIST']}

    #=============================================================
    # OUTPUT EFFECT
    #=============================================================

    #**** Flag of the type of output map (integers):
    #         (0)  3D electron density
    #         (1)  2D elctron Density
    #         (2)  Temperature
    #           (0) mass weighted
    #           (1) Emission weighted temperature (e.g. former (3))
    #               (set also Energy cutoff bands below!)
    #           (2) Emission weighted temperature (with tabulated cooling function)
    #           (3) Spectroscopic (e.g. former (13))
    #           (4) Emission temperature (by Elena) (e.g. former (14))
    #         (4)  Volume of particles
    #         (5)  age of particles
    #         (6)  X-ray SB of particles
    #           (0) simple sqrt(T)
    #              or  with cooling function from tables
    #           (1)       - energy band [0.1-10] keV (bolometric)
    #           (2)       - energy band [0.1-2.4] keV
    #           (3)       - energy band [0.3-3.5] keV
    #           (4)       - user defined energy band [ XRAY_E0 - XRAY_E1 ] keV (define below)
    #         (7)  Compton y-parameter: tSZ effect
    #             flag of the subtype of output map:
    #           (0) adimensional [DI/I]
    #           (1) adimensional [DT/T]
    #           (2) adimensinal [DT/T] with relativistic correction
    #           (3) n_e T^2, e.g. (1)*T
    #         (8)  Kinematical Sunjaev Zeldovich Effect of particles: kSZ efect
    #             flag of the subtype of output map:
    #           (0) adimensional [DI/I]
    #           (1) adimensional [DT/T]
    #         (9)  Q & U Stokes parameters for transvers particles motion: kpSZ effect
    #	 	      flag of the subtype of output map:
    #           (0) adimensional [Q_I & U_I] - noramlized on the CMB intensity
    #           (1) adimensional [Q_T & U_T] - noramlized on the CMB temperature
    #        (10) Mass weighted velocity
    #        (11) Mass weighted squared velocity
    #        (12) Magnetic field related maps
    #            flag of the subtype of output map:
    #           (0) RM
    #           (1) <B_x,y,z>
    #           (2) <|B|>
    #           (3) Synchrotron Emission
    #          (31) Synchrotron Polarised Emission & Polarisation Angle
    #          (10) RM with relativistic corrections
    #          (40) Gamma Emission [# of photons / cm^2] @ E_gam
    #          (41) Gamma Emission [erg/cm^2  ]			 @ E_gam
    #          (42) Gamma Emission [# of photons / cm^2] integrated over [E_pi_min,E_pi_max]
    #        (15) Metalicity map (of METAL_SPECIES)
    #            flag of the subtype of output map:
    #           (0) Mass weighted
    #           (1) Emission weighted
    #        (16) Rees-Sciama Effect by mass transverse motion [phi_dot]
    #        (17) Mach-Number
    #            flag of the subtype of output map:
    #           (0) Volume weighted
    #           (1) Energy weighted
    #           (2) Mass weighted
    #        (18) Y_X parameter (mass weighted)
    #        (19) Y_X parameter (pasquale weighted)
    #        (20) Projected potential
    #        (21) BP_CR Protons
    #           (0) Pressure
    #           (1) X_crp ( P_crp/P_th )
    #        (22) BP_CR Electrons
    #           (0) Pressure
    #           (1) Number density
    #           (2) Energy density
    #           (3) Synchrotron Emission
    #          (31) Synchrotron Polarised Emission & Polarisation Angle
    #       (100) 3D DM density
    #       (101) 2D DM density
    #       (102) Rees-Sciama Effect by mass transverse motion [projected momentum]
    #             flag of the subtype of output map:
    #           (0) by integral
    #           (1) by FFT
    #       (103) Gravitational lensing: Deflection angle
    OUTPUT_MAP = {kwargs['OUTPUT_MAP']}
    OUTPUT_SUB = {kwargs['OUTPUT_SUB']}

    # Define metal species
    #	 	flag of the subtype of output map:
    #		(0) total
    #               (1) > Fe
    #               (2) C
    #               (3) N
    #               (4) O
    #               (5) Mg
    #               (6) Si
    #               (7) Fe
    METAL_SPECIES = {kwargs['METAL_SPECIES']}

    #**** Unit of output Temperatures
    #        (0) Kelvin
    #        (1) KeV
    TEMP_UNIT = {kwargs['TEMP_UNIT']}

    #**** Ignore particles below this value [in keV]
    TEMP_CUT = {kwargs['TEMP_CUT']}

    #**** Energy cut-off for X-ray (Flag of the type = 6):
    #        ene_a, (case(6), case(2)-subtype 1) lower energy bound for X-ray image [keV]
    #        ene_b, (case(6), case(2)-subtype 1) upper energy bound for X-ray image [keV]
    XRAY_E0 = {kwargs['XRAY_E0']}
    XRAY_E1 = {kwargs['XRAY_E1']}

    # Path to the tables for the cooling function
    X_TABLES = {kwargs['X_TABLES']}

    #**** Input observational frequency for SZ maps [GHz] (real):
    FRQ_ZS = {kwargs['FRQ_ZS']}

    #**** Input for Radio Maps
    #XCRP_TABLE_PATH gives the change of XCRP over radius.
    ECRP_MIN = {kwargs['ECRP_MIN']}
    FRQ_Pnu = {kwargs['FRQ_Pnu']}
    XCRP = {kwargs['XCRP']}
    XCRP_TABLE_PATH = {kwargs['XCRP_TABLE_PATH']}
    B_TABLE_PATH = {kwargs['B_TABLE_PATH']}
    KERNEL_TABLE_PATH = {kwargs['KERNEL_TABLE_PATH']}
    GAM_nu = {kwargs['GAM_nu']}
    #XCRP_TABLE_FILE_FMT = {kwargs['XCRP_TABLE_FILE_FMT']}


    #**** Input for Gamma Maps
    #energies in GeV
    E_gam = {kwargs['E_gam']}
    E_gam_min = {kwargs['E_gam_min']}
    E_gam_max = {kwargs['E_gam_max']}


    #**** Input for LMB_SPECTRAL_CRs
    # Number of momentum bins for CR model
    CR_nbins = {kwargs['CR_nbins']}
    CR_pmin =  {kwargs['CR_pmin']}
    CR_pmax = {kwargs['CR_pmax']}
    CR_subsamples = {kwargs['CR_subsamples']}
    CR_DSlope = {kwargs['CR_DSlope']}


    #**** Set to 1 if you want additional statistical informations (L_x,T,...).
    #         Printed on screen (integer):
    #         (0) disabled
    #         (1) enabled
    GIVE_MORE_INFO = {kwargs['GIVE_MORE_INFO']}

    #=============================================================
    # REAST FRAME - PROJECTION
    #=============================================================

    #**** Flag of the image projection (integer), from - to:
    #         (1) along z, xy plane
    #         (2) along y, xz plane
    #         (3) along x, yz plane
    #         (4) along all 3 axis
    PROJECT       = {kwargs['PROJECT']}

    #**** Flag of the cluster baricenter motion (integer):
    #         (0) disabled
    #         (1) enabled: at the particles velocity is sobstituted the cluster baricenter velocity
    #         (2) enabled: paricles velocity in the cluster baricenter rest frame
    CENTER_MOTION = {kwargs['CENTER_MOTION']}

    #**** Flag to subtract the velocity of the local group (integer):
    #         (0) disabled
    #         (1) enabled
    REMOVE_LOCAL_GROUP_VEL = {kwargs['REMOVE_LOCAL_GROUP_VEL']}

    #**** Flag of the definition of the center (integer):
    #         (0) barycenter of gas particles
    #         (1) selected by user
    #         (2) read by a file
    CENTER = {kwargs['CENTER']}

    #**** Center position/Cluster data for 'select by user' (reals [Mpc]):
    CENTER_X = {kwargs['CENTER_X']}
    CENTER_Y = {kwargs['CENTER_Y']}
    CENTER_Z = {kwargs['CENTER_Z']}

    #**** Flag for PERIODIC boxes (integer):
    #         (0) not periodic
    #         (1) periodic
    PERIODIC = {kwargs['PERIODIC']}

    #**** Flag for building lightcones (integer):
    LIGHTCONE = {kwargs['LIGHTCONE']}
    #**** Position of box in lightcone (reals [Mpc]):
    X_ORIGIN = {kwargs['X_ORIGIN']}
    Y_ORIGIN = {kwargs['Y_ORIGIN']}
    Z_ORIGIN = {kwargs['Z_ORIGIN']}
    #**** Opening angle of the lightcone (degree):
    OPEN_ANGLE = {kwargs['OPEN_ANGLE']}


    #**** File of the centers of the images for 'read by a file':
    MAIN_PROG = {kwargs['MAIN_PROG']}


    #=============================================================
    # OUTPUT FILE
    #=============================================================

    #**** Flag string of output file selection:
    #         (1) binary image (IDL)
    #             flag of the subtype of output map:'
    #             (1) header in txt format'
    #             (2) header in FITS format'
    #         (2) ASCII image 
    #         (3) FITS image
    FILE_FORMAT = {kwargs['FILE_FORMAT']}

    #**** Flag for producing HEADER
    #         (0) no header 
    #         (1) header 
    #             for idl /asi format header will be created as ascii file
    FILE_HEADER = {kwargs['FILE_HEADER']}

    #**** Prefix of the output files:
    PREFIX_OUT = {kwargs['PREFIX_OUT']}

    #=============================================================
    # COSMOLOGY (for distances computations)
    #=============================================================

    #**** Cosmology to use
    HUBBLE = {kwargs['HUBBLE']}
    OMEGA  = {kwargs['OMEGA']}
    LAMBDA = {kwargs['LAMBDA']}
    """)




def read_halo_positons(path = './', snapNum = 36, format = 1, haloID = 0, subfileNumberAndID = [0,0] ):

    import g3read
    import g3matcha as matcha



    # --- ID dictionary ---

    subnum = 2
    haloid = 8
    # --- We read the data from the SubFind ---
    if snapNum < 10:
        groupbase= path + 'groups_00'+str(snapNum)+'/sub_00'+str(snapNum)
    else:
        groupbase= path + 'groups_0'+str(snapNum)+'/sub_0'+str(snapNum)
    
    if format == 1:

        for halo in matcha.yield_haloes(groupbase, 0,  blocks=('MVIR','RVIR','RCRI','R5CC','GPOS')):
            if halo['ihalo'] == haloID:
                halopos=halo['GPOS']
                halo_rvir=halo['RVIR']
                halo_r200c=halo['RCRI']
                halo_r500c=halo['R5CC']
                halo_mvir=halo['MVIR']
                break
    
    if format == 2:
        fof = g3read.read_new(groupbase+'.'+str(subnum), ['GPOS','RVIR','RCRI','R5CC','MVIR'], 0, is_snap=False )
        halopos = fof['GPOS'][haloid] 
        halo_rvir=fof['RVIR'][haloid]
        halo_r200c=fof['RCRI'][haloid]
        halo_r500c=fof['R5CC'][haloid] 
        halo_mvir = fof['MVIR'][haloid]

    else: 
        print('ERROR: WRONG FORMAT')

    return(halopos, halo_mvir, halo_rvir, halo_r200c, halo_r500c )



def plot_haloFromFits(fileName, cmap="viridis", vmin=None, vmax=None, xlabel='[Mpc]', ylabel= '[Mpc]', barlabel='', radius1=0, radius2=0, radius3=0):

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from astropy.io import fits

    # Fits handle
    hdu = fits.open(fileName)[1] # This opens the fits file
    image_data = hdu.data  # This loads the 2D array of image pixels

    # Matplotlib plotting
    fig, ax = plt.subplots(1, 1, figsize=(8.8,6.8), dpi=100) # Opens a Figure instance
    R = hdu.header["BOX_KPC"]/2./1000. # This reads from the header of the fits file the size of the image in physical coordinates
    # The following line produces the actual plot

    im = ax.imshow(image_data[:,:], origin='lower', cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax), aspect='equal', extent=(-R,R,-R,R))
    # fonts
    font1 = {'family':'serif','fontname':'DejaVu Sans','color':'black','size':12}
    font2 = {'family':'serif','fontname':'DejaVu Sans','color':'black','size':10}  
    # labels and plot title
    ax.set_xlabel(xlabel,fontdict=font2)
    ax.set_ylabel(ylabel,fontdict=font2)
    axcbar=fig.colorbar(im, shrink=1, aspect=20) # This add the colorbar
    axcbar.set_label(barlabel)

    # plot virial radius

    circle=plt.Circle((0,0),radius1/1000, color='black', linewidth=0.01, linestyle='--',fill=False)
    ax = fig.gca()
    ax.add_patch(circle)

    # plot r200

    circle=plt.Circle((0,0),radius2/1000, color='r', linewidth=0.01, linestyle='--',fill=False)
    ax = fig.gca()
    ax.add_patch(circle)

    # plot r500

    circle=plt.Circle((0,0),radius3/1000, color='b', linewidth=0.01, linestyle='--',fill=False)
    ax = fig.gca()
    ax.add_patch(circle)

    result1=subprocess.run(["mkdir plots"],shell=True)
    plt.savefig(fileName+'.pdf')
    result2=subprocess.run(["mv *.pdf ../plots/"],shell=True)


def plot_skyMap(path='./',filename="smac_tSZ_032.h.fits",figurename='SkyMap.pdf'):

    import numpy as np
    import math
    import healpy as hp
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap


    NSIDE = 32  # resolution of the map. Usually a power of 2

    print(
        "Approximate resolution at NSIDE {} is {:.2} deg".format(
            NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
        )
    )

    NPIX = hp.nside2npix(NSIDE) # number of pixels in the of the map
    print(NPIX)

    my_map_I = hp.read_map(path+filename)

    hp.visufunc.mollview(
        my_map_I,
        #coord=["E"], # ‘G’, ‘E’ or ‘C’ = Galactic, Ecliptical, Ecuatorial. coord=["G", "E"] transforms from galactic to ecliptical coordinates
        #smac produces the fits file in supergalactic coordinates!
        title="Full Sky Map [SGB]",
        unit="thermal SZ",
        norm='log', #another possibility is "hist"
        min=math.modf(np.amin(my_map_I))[0],
        max=math.modf(np.amax(my_map_I))[0],
        cmap='jet', format='%.2e'
    )
        #just copy from supermuc
    plt.savefig(figurename)


def remove_files(fileNames):
    import os
    os.remove(fileNames)


def  create_smac_runFile(path = './', fileName = 'runSmac_iterat.sh', parameterFile = 'smac_1.inp', executable = 'Smac_6.1_mpi', executablePath = './'):
    f = open(path + fileName, "w")
    f.write(
    f'''module load slurm_setup
module load cfitsio

{executablePath}{executable} \t {parameterFile}''')


def read_fitsfilesName(path = './',parameterFile = 'smac_1.inp'):

    with open(path+parameterFile, "r") as f:
        lines = f.readlines()  #in the matrix lines is all the information of the file

    # -- modify the information of the old param file --
    for i in range(len(lines)):
        if 'PREFIX_OUT' in lines[i]:
            prefix_out = lines[i]
            prefix_out = prefix_out[13:].rstrip()
        if 'SNAP_START' in lines[i]:
            txt = lines[i]
            [int(s) for s in txt.split() if s.isdigit()]
            snap_start = [int(s) for s in txt.split() if s.isdigit()]
            #snap_start = snap_start[0]
            
        if 'PROJECT' in lines[i]:
            if lines[i] == 'PROJECT = 4\n':
                project = 4
            if lines[i] == 'PROJECT = 3\n':
                project = 3
            if lines[i] == 'PROJECT = 2\n':
                project = 2
            if lines[i] == 'PROJECT = 1\n':
                project = 1

    
    if snap_start[0] > 100:
        snap = str(snap_start[0])
    if 9 < snap_start[0] < 100:
        snap ='.0'+ str(snap_start[0])
    if snap_start[0] < 9:
        snap ='.00'+ str(snap_start[0])

    filenames=[]
    if project == 4:
        filenames.append(prefix_out+snap+'.a.x.fits')
        filenames.append(prefix_out+snap+'.a.y.fits')
        filenames.append(prefix_out+snap+'.a.z.fits')
    if project  == 1:
        filenames.append(prefix_out+snap+'.a.z.fits')
    if project  == 2:
        filenames.append(prefix_out+snap+'.a.y.fits')
    if project  == 3:
        filenames.append(prefix_out+snap+'.a.x.fits')

    return(filenames)
    
    
#        for j in range(len(para_name)):
#            if para_name[j] in lines[i]:
#                
#                if parameters[para_name[j]] == None:
#                    continue
#
#                lines[i]= para_name[j] + " = " + str(parameters[para_name[j]]) + "\n"
#    #prefix_out + snapnumber + projection +.fits


def run_smac(path = './',fileName = 'runSmac_iterat.sh'):
    import subprocess
    import os

    fileNames = read_fitsfilesName(fileName)
    os.remove(fileNames)
    result1=subprocess.run(["/bin/sh", path + fileName],stderr=subprocess.PIPE, text=True)
    result2=subprocess.run(["mkdir fitsfiles"], shell=True)
    result3=subprocess.run(["mv *.fits fistfiles"],shell=True)





''' ////////////////////////////////////////////////////////////////'''
''' ============================== MAIN ============================'''
''' ////////////////////////////////////////////////////////////////'''

''' 
data_structure =[ halo_name1, [subfile1, array_position1]),
                  halo_name2, [subfile2, array_position2]),
                .
                .
                . 
                    ]
'''
'''   
cluster_storage = [("Coma", [1, 0]),
                   ("Virgo", [2, 27] ),
                   ("Perseus", [2, 8]),
                   ("Centaurus", [2, 13]),
                   ("A119", [2, 11]),
                   ("A539", [4, 0]),
                   ("A576", [5, 95]),
                   ("A1185", [5, 132]),
                   ("A2256", [8, 7]),
                   ("Ophiucus", [9, 309]),
                   ("A2147", [10, 51]),
                   ("A2877", [10, 209]),
                   ("A3581", [13, 316]),
                   ("Norma", [23, 276]),
                   ("A1367", [18, 169]),
                   ("Fornax", [24, 369]),
                   ("A2319", [38, 426])]



snapnumber = 36
sim_path = '/dss/dssfs02/pr62go/pr62go-dss-0001/Local/1536/agn/'
parameterfile = 'smac_1.par' #default name
parameterfile = 'test.inp'
smacfolder = '/hppfs/work/pr86re/di67map2/1536/SMAC/RunSmac/'
smacfolder = '/gpfs/work/pn68va/di67map/LOCALUNIVERSE/development/SMACPy/test/'

for i in range(len(cluster_storage)):


   haloname = cluster_storage[i][0]
   
   # --- get halo positions ---
   halopos,  halo_mvir, halo_rvir, halo_r200c, halo_r500c = read_halo_positons(snapNum=snapnumber,
                                                                                path = sim_path, 
                                                                                format = 2, 
                                                                                subfileNumberAndID=cluster_storage[i][1])

    
   # --- create new parameterfile with the changes that we want ---

   prefix_out = str(haloname) + '_density_'

   changes = dict(eraseFits=True,
                   CENTER_X = halopos[0]/1000, 
                   CENTER_Y = halopos[1]/1000, 
                   CENTER_Z = halopos[2]/1000,
                   OUTPUT_MAP=1, 
                   PREFIX_OUT = prefix_out
                   )
   
   create_smac_par(path=smacfolder,**changes, fileName=parameterfile)

   # --- run smac ---
   create_smac_runFile(path=smacfolder)

   run_smac(path=smacfolder)

   ## --- produce pdf plot ---
   fitsfiles = read_fitsfilesName()
   for i in range(len(fitsfiles)):
       plot_haloFromFits('fitsfiles/'+fitsfiles[i], cmap="viridis", barlabel= 'density' , radius1=halo_rvir, radius2=halo_r200c,  radius3=halo_r500c)
       
'''
