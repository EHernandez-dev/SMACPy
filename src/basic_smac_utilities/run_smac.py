def  create_smac_runFile(path = './', fileName = 'runSmac_iterat.sh', parameterFile = 'smac_1.inp', executable = 'Smac_6.1_mpi', executablePath = './'):
    f = open(path + fileName, "w")
    f.write(
    f'''module load slurm_setup
module load cfitsio

{executablePath}{executable} \t {parameterFile}''')
    
def run_smac(path = './',fileName = 'runSmac_iterat.sh', parfile = 'smac_1.inp', parfilepath='./'):
    import subprocess
    import os
    from manage_files import read_fitsfilesName

    fileNames = read_fitsfilesName(parfilepath,parfile)
    os.remove(fileNames)
    #result1=subprocess.run(["/bin/sh", path + fileName],stderr=subprocess.PIPE, text=True)
    #result2=subprocess.run(["mkdir fitsfiles"], shell=True)
    #result3=subprocess.run(["mv *.fits fistfiles"],shell=True)
    
run_smac(path =  '/gpfs/work/pn68va/di67map/LOCALUNIVERSE/1536/agn/RunSmac/')