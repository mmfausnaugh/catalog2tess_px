import scipy as sp
import matplotlib.pyplot as plt
from astropy.io import fits
import glob

F = plt.figure()
ax1 = F.add_subplot(111,projection='mollweide')
colors = ['k','purple','b','c','xkcd:forest green','g','xkcd:lime green','yellow','orange','r','m','xkcd:light purple','pink',
          'k','purple','b','c','xkcd:forest green','g','xkcd:lime green','yellow','orange','r','m','xkcd:light purple','pink']
for ii in range(26):
    flist = glob.glob('s{:02d}/*'.format(25-ii))
    for f in flist:
        ra,dec = sp.genfromtxt(f,unpack=1,usecols=(1,3))

        ax1.scatter(ra*sp.pi/180 - sp.pi, dec*sp.pi/180,color=colors[25-ii])

plt.show()
