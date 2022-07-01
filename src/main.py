from time import sleep
import subprocess
from basic_smac_utilities import plot_haloFromFits
from basic_smac_utilities.manage_files import read_fitsfilesName
from basic_smac_utilities.manage_files import remove_files
from basic_smac_utilities.read_halo_positions import read_halo_positions
from basic_smac_utilities.create_smac_par import create_smac_par
from basic_smac_utilities.run_smac import create_smac_runFile
from basic_smac_utilities.run_smac import run_smac


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
   halopos,  halo_mvir, halo_rvir, halo_r200c, halo_r500c = read_halo_positions(snapNum=snapnumber,
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