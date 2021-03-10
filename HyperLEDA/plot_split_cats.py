import scipy as sp
import matplotlib.pyplot as plt
from astropy.io import fits
import glob

F = plt.figure()
ax1 = F.add_subplot(111,projection='mollweide')

colors = ['k','purple','b','c','xkcd:forest green','g','xkcd:lime green','yellow','orange','r','m','xkcd:light purple','pink',
          'k','purple','b','c','xkcd:forest green','g','xkcd:lime green','yellow','orange','r','m','xkcd:light purple','pink']
for ii in range(26):
    flist = glob.glob('s{:02d}/*'.format(ii+1))
    for f in flist:
        try:
            ra,dec = sp.genfromtxt(f,unpack=1,usecols=(5,6),delimiter=',')
            ra *= 180.0/12.0
        except:
            continue
        ax1.scatter(ra*sp.pi/180 - sp.pi, dec*sp.pi/180,color=colors[ii])

plt.show()
