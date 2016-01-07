'''
Try converting the IDL procedure to convert a contrast curve to mass 
detection limits to python.
'''
import numpy as np
import pylab as plt
from scipy.interpolate import interp1d
import glob
from scipy.integrate import simps
import os
import cPickle as pickle


class con2mass:

    def __init__(self, title, separation, contrast, sigma, passband, 
                 starmag, distpc, starage):
        '''separation = array of separations (arcsec)
        contrast = contrast curve
        sigma = sigma of the contrast curve (usually 3 or 5)
        passband = passband of the contrast curve (H, Ks, or Lp)
        starmag = apparent magnitude of the star in the passband
        distpc = distance to the star (parsec)
        starage = age of the system (Myr)'''

        # Add stuff
        self.title = title
        good = np.where(np.isfinite(contrast))[0]
        self.sep = separation[good]
        self.con = contrast[good]
        self.starmag = starmag
        self.distpc = distpc
        self.starage = starage
        self.passband = passband
        self.define_passband_stuff()
        
        # Compute stuff
        self.distau = self.distpc * 2.0626e5
        self.planetmag = self.con + self.starmag
        self.absmag = self.planetmag - 5.*(np.log10(self.distpc)-1)
        self.projsep = np.tan(np.deg2rad(self.sep/(60.*60))) * self.distau

        # Get cooling models
        self.cooling_models = np.zeros(0)
        self.get_BTSettl()
        self.get_Spiegel()

        # Compute mass detection limits using IDL interpol function
        self.get_masslim1(BTsettl=1)
        self.get_masslim2(BTsettl=1)
        self.get_masslim1(Spiegel=1)
        self.get_masslim2(Spiegel=1)


        # Make plot
        self.make_plot()


    def define_passband_stuff(self):
        '''Get the central wavelength of the passband in microns and 
        its zeropoint flux density in Jy. Zeropoints are from 
        http://www.astro.princeton.edu/~burrows/warmstart/
        '''
        if self.passband == 'H':
            self.wl = 1.6
            self.fzero = 1040.
        elif self.passband == 'Ks':
            self.wl = 2.2
            self.fzero = 645.
        elif self.passband == 'Lp':
            self.wl = 3.8   
            self.fzero = 252.

        
    def get_BTSettl(self):
        '''Get the BTSettl photometric model at a given age.'''
        self.cooling_models = np.append(self.cooling_models, 'BTSettl')
        d = np.loadtxt('cooling_models/BTsettl%iMyr.dat'%self.starage,
                       skiprows=4)
        self.BTmp = d[1:,0]
        self.BTJ  = d[1:,7]
        self.BTH  = d[1:,8]
        self.BTKs = d[1:,9]
        self.BTLp = d[1:,10]
        self.BTMp = d[1:,11]
        if self.passband == 'H':
            self.BTabsmag = d[1:,8]
        elif self.passband == 'Ks':
            self.BTabsmag = d[1:,9]
        elif self.passband == 'Lp':
            self.BTabsmag = d[1:,10]   
        # Convolve with filter
        self.filtconv_jy(BTsettl=1)


    def get_Spiegel(self):
        '''Get photometric models from Speigel and Burrows.''' 
        self.cooling_models = np.append(self.cooling_models, 'Spiegel')
        fnames = np.array(glob.glob('adamspec/spec_hy1s_mass_*_age_%.4d.txt'%self.starage))
        num = fnames.size
        mps = np.zeros(num)
        absmag = np.zeros(num)
        SpiegJ = np.zeros(num)
        SpiegH = np.zeros(num)
        SpiegKs = np.zeros(num)
        SpiegLp = np.zeros(num)
        SpiegMp = np.zeros(num)
        for i in range(num):
            mps[i] = float(fnames[i].split('mass_')[1].split('_age')[0])
            d = np.loadtxt(fnames[i])
            wls = d[0,1:]
            entropy = d[1:,0]
            index = np.where(entropy == max(entropy))[0][0]
            fnu = d[index,1:]  # mJy
            finterp = interp1d(wls, fnu)
            fnu = finterp(self.wl)
            absmag[i] = -2.5*np.log10(1e-3*fnu/self.fzero) 
            # Get over mags
            SpiegJ[i] = -2.5*np.log10(1e-3*finterp(1.3)/self.fzero)
            SpiegH[i] = -2.5*np.log10(1e-3*finterp(1.6)/self.fzero)
            SpiegKs[i] = -2.5*np.log10(1e-3*finterp(2.2)/self.fzero)
            SpiegLp[i] = -2.5*np.log10(1e-3*finterp(3.8)/self.fzero)
            SpiegMp[i] = -2.5*np.log10(1e-3*finterp(4.8)/self.fzero)
                        
        self.SpiegJ  = SpiegJ
        self.SpiegH  = SpiegH
        self.SpiegKs = SpiegKs
        self.SpiegLp = SpiegLp
        self.SpiegMp = SpiegMp
        self.Spiegmp = mps
        self.Spiegabsmag = absmag
        # Convolve with filter
        self.filtconv_jy(Spiegel=1)        

            
    def filtconv_jy(self, dbnu=0, BTsettl=False, Spiegel=False):
        '''Convolve BTsettl model with filter transmission.'''
        if BTsettl:
            newabsmag = np.zeros(self.BTmp.size)
        if Spiegel:
            newabsmag = np.zeros(self.Spiegmp.size)        

        # J, H, Ks, Lp, Mp
        wl = np.array((1.3,1.6,2.2,3.8,4.8))
        for i in range(newabsmag.size):
            if BTsettl:
                mags = np.array((self.BTJ[i], self.BTH[i], self.BTKs[i], 
                                 self.BTLp[i], self.BTMp[i]))
                fnu = self.fzero * 1e3 * 10**(mags/(-2.5))   # mJy
            if Spiegel:
                mags = np.array((self.SpiegJ[i], self.SpiegH[i], 
                                 self.SpiegKs[i], self.SpiegLp[i], 
                                 self.SpiegMp[i]))
                fnu = self.fzero * 1e3 * 10**(mags/(-2.5)) # mJy

            # Get the filter transmission curve
            d = np.loadtxt('filters/'+self.passband+'band.dat')
            filtf  = d[:,1]
            cutoff = .1
            good = np.where(filtf > cutoff*max(filtf))[0]
            filtwl = d[:,0][good]*1e-3   # microns
            filtf  = filtf[good]
            cenwl  = np.median(filtwl[filtf > cutoff*max(filtwl)])

            # Convert phiwl to phinu
            nnu = filtwl.size
            phi = np.zeros(nnu)
            nufilt = np.zeros(nnu)
            for j in range(nnu):
                nufilt[j] = 2.9979e14 / filtwl[j]
                if dbnu == 0:
                    phi[nnu-j-1] = filtf[j]
                elif dbnu == 1:
                    phi[nnu-j-1] = filtf[j] / nufilt[nnu-j-1]
            nu = 2.9979e14 / wl
            sort = np.argsort(nu)
            finterp = interp1d(nu[sort], fnu[sort], kind='linear')
            fnew = finterp(nufilt)
            # Compare these lines to int_tablulated in IDL
            ftop = simps(fnew*phi, nufilt)
            fbot = simps(phi, nufilt)
            fluxout = ftop/fbot
            newabsmag[i] = -2.5*np.log10(1e-3*fluxout/self.fzero)
        if BTsettl:  # This one works
            self.BTabsmag = newabsmag
        if Spiegel:
            self.Spiegabsmag = newabsmag


    def get_masslim1(self, BTsettl=False, Spiegel=False):
        '''Convert the planet absolute magnitude from the contrast curve
        to mass detection limits using the cooling models and the 
        IDL interpol function (see IDLinterpol.pro).'''
        if BTsettl:
            foridl1 = np.zeros((self.BTabsmag.size, 2))
            foridl1[:,0] = self.BTmp
            foridl1[:,1] = self.BTabsmag
        if Spiegel:
            foridl1 = np.zeros((self.Spiegabsmag.size, 2))
            foridl1[:,0] = self.Spiegmp
            foridl1[:,1] = self.Spiegabsmag

        # Create files for IDL
        np.savetxt('foridl_mass_mag.dat', foridl1, delimiter='\t',
                   fmt='%.6f', header='Mass, Abs Mag')
        np.savetxt('foridl_planetmag.dat', self.absmag, delimiter='\t',
                   fmt='%.6f', header='Planet Abs Mag')

        print 'Instructions for computing the mass limits:\n'
        print 'In a new terminal window open IDL. Compile and run'
        print 'IDLinterpol.pro with no arguments. This procedure'
        print 'will save a file of mass limits to be read in by the'
        print 'python function self.get_masslim2().\n'
        go = 0
        while go != 1:
            try:
                go = input('Enter 1 when IDL function has run successfully. ')
            except SyntaxError:
                go = 0


    def get_masslim2(self, BTsettl=False, Spiegel=False):
        '''Read in the file containing the mass detection limits and 
        save them to the object in MJup. Also clean stuff.'''
        d = np.loadtxt('forpython_planetmass.dat')
        if BTsettl:
            self.BTmasslim = d*1e3
        if Spiegel:
            self.Spiegmasslim = d
        # Clean up after yourself
        os.system('rm for*.dat')

        
    def make_plot(self):
        '''Plot the mass detection limits.'''
        # Plot BTsettl curve
        plt.plot(self.projsep, self.Spiegmasslim, 'b-', lw=2,
                 label='Spiegel & Burrows')
        plt.plot(self.projsep, self.BTmasslim, 'r-', lw=2,
                 label=' BT-Settl')
        
        plt.xlabel('Projected Separation (AU)')
        plt.ylabel(' Planet Mass (M$_{\mathrm{Jup}})$')
        plt.minorticks_on()
        plt.legend(loc='upper right')
        plt.savefig('plots/'+self.title+'.png')
        plt.show()



if __name__ == '__main__':
    d = np.loadtxt('data/hr8799contrast_curve.dat')
    data = con2mass('HR8799', d[:,0], d[:,1], 3, 'Lp', 5.22, 39.4, 30)
    f=open('pickle','wb')
    pickle.dump(data, f)
    f.close()
