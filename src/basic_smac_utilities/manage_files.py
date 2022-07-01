def remove_files(fileNames):
    import os
    os.remove(fileNames)
    
    
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