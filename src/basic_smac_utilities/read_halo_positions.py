
def read_halo_positions(path = './', snapNum = 36, format = 1, haloID = 0, subfileNumberAndID = [0,0] ):

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