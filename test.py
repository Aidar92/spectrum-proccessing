# cython: language_level=2
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.table import Table
import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
from scipy import signal
import matplotlib.font_manager as font_manager


bta = EarthLocation(lat=43.65*u.deg, lon=41.3*u.deg, height=2070*u.m)
time = Time('2012-7-12 23:00:00')

plt.rc('font', size=10)
#Azimuth and zenit distance calculate
def AzAlt(star_center, time, location=bta):
    coord = star_center.transform_to(AltAz(obstime=time, location=location))
    return coord

#Radial velocity correction calculate
def RVcorr(star_center, time):
    heliocorr = star_center.radial_velocity_correction(
        'heliocentric', obstime=time, location=bta)
    return heliocorr.to(u.km/u.s)

#Cross-correlation calculate for potasium spectrum on star
def crossCor():
    temlpateFile = fits.open('fits/e665001s.fits')
    tData = pd.DataFrame(temlpateFile[0].data.T)
    trFile = fits.open('fits/e665008s.fits')
    trData = pd.DataFrame(trFile[0].data.T)
    x = np.arange(0, trData.shape[0])
    xi = np.linspace(0, trData.shape[0], trData.shape[0]*50)
    for i in range(trData.shape[1]):
        tf = np.interp(xi, x, tData[i][:])
        df = np.interp(xi, x, trData[i][:])
        corr = signal.correlate(tf, df, mode='same')
        maxind = np.argmax(corr)
        xc = np.linspace(-xi[maxind-1]*2, xi[maxind]*2, trData.shape[0]*50)
        fig, (ax, ax_corr) = plt.subplots(2, 1)
        ax.plot(xc, tf)
        ax.plot(xc, df)
        ax_corr.plot(xc, corr)
        print(xc[maxind], corr[maxind])
        ax_corr.plot(xc[maxind], corr[maxind], 'ro')
        #plt.savefig('fits/img/%d.png' % (i+1), dpi=900)
        plt.show()
        print(i)
        plt.close()


if __name__ == "__main__":
    crossCor()
