import scipy as sp
import matplotlib.pyplot as plt
from astropy.io import fits
import glob

F = plt.figure()
ax1 = F.add_subplot(111,projection='mollweide')
colors = ['brown','yellow','k','b','c','m','r','orange','purple','pink','g','xkcd:light purple','xkcd:forest green']
for ii in range(13):
    flist = glob.glob('s{:02d}/*'.format(ii+1))
    for f in flist:
        d = fits.open(f)
        ra = d[1].data['RA']
        dec = d[1].data['DEC']

        ax1.scatter(ra*sp.pi/180 - sp.pi, dec*sp.pi/180,color=colors[ii])

plt.show()
