import netCDF4 as nc            # Reading NetCDF4 files.
import numpy   as np            # Array operations.
import matplotlib.pyplot as plt # Plot data
import sys

from matplotlib import rc       # Globally setup figure options
rc('text',       usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font',       family='serif', size=12)
rc('grid',       linestyle='dotted')
rc('axes',       grid=True)

# At home, screen 27''
rc('figure',     dpi=200)
rc('savefig',    dpi=100)

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 variable list-of-files.")
    quit()

var = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

plt.figure(figsize=(5,4))
axes = plt.gca()

for file in setoffiles:
    data = nc.Dataset(file, 'r')
    data.set_auto_mask(False)

    t=data.variables['t'][:]
    z=data.variables['z'][:]
    f=data.variables[var][:,:] # the first index is time, the second is vertical node

    if np.size(t) == 1:
        plt.plot(f[0,:],z,label=r'{}'.format(file))
    else:
        # its = [ 0 ]                                         # Plot just the first time
        its = list(range(0,np.size(t),int(np.size(t)/5)))  # plot 5 profiles in the given range of times
        for it in its:
            plt.plot(f[it,:],z,label=r'It = {} in {}'.format(it,file))

plt.xlabel(r'{}'.format(var))
plt.ylabel(r'height $z$')
plt.legend(loc='best')
axes.spines['right'].set_visible(False)
axes.spines['left'].set_position(('axes',-0.01))
axes.get_yaxis().tick_left()
axes.spines['top'].set_visible(False)
axes.spines['bottom'].set_position(('axes',-0.01))
axes.get_xaxis().tick_bottom()
plt.tight_layout(pad=0.1)

plt.savefig('figure1.pdf',bbox_inches='tight')
plt.show()
