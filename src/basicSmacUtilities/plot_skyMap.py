
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