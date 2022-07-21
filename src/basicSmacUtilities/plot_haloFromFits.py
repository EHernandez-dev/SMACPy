
def plot_haloFromFits(fileName, cmap="viridis", vmin=None, vmax=None, xlabel='[Mpc]', ylabel= '[Mpc]', barlabel='', radius1=0, radius2=0, radius3=0):

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from astropy.io import fits
    import subprocess

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