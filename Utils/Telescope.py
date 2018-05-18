from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import Bandpass,Sed
from Throughputs import Throughputs

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.constants import *

class Telescope(Throughputs):
    def __init__(self,airmass=1,atmos=True,aerosol=True,**kwargs):
        Throughputs.__init__(self,**kwargs)

        params=['mag_sky','m5','FWHMeff','Tb','Sigmab','zp','counts_zp','Skyb','flux_sky']

        self.data={}
        for par in params:
            self.data[par]={}
        
        self.data['FWHMeff']=dict(zip('ugrizy',[0.92,0.87,0.83,0.80,0.78,0.76]))
         
        self.atmos=atmos

        self.Load_Atmosphere(airmass)
        
    @property
    def FWHMeff(self):
        return self.data['FWHMeff']

    @property
    def mag_sky(self,band='*'):
        if band not in self.data['mag_sky'].keys():
            self.Calc_Inputs(band)
        if band == '*':
            return self.data['mag_sky']
        else:
            return self.data['mag_sky'][band]

    @property
    def m5(self,band='*'):
        if band not in self.data['m5'].keys():
            self.Calc_m5(band)
        if band == '*':
            return self.data['m5']
        else:
            return self.data['m5'][band]

    @property
    def Tb(self,band='*'):
        if band not in self.data['Tb'].keys():
            self.Calc_Inputs(band)
        if band == '*':
            return self.data['Tb']
        else:
            return self.data['Tb'][band]

    @property
    def Sigmab(self,band='*'):
        if band not in self.data['Sigmab'].keys():
            self.Calc_Inputs(band)
        if band == '*':
            return self.data['Sigmab']
        else:
            return self.data['Sigmab'][band]
    @property
    def zp(self,band='*'):
        if band not in self.data['zp'].keys():
            self.Calc_zp(band)
        if band == '*':
            return self.data['zp']
        else:
            return self.data['zp'][band]

    @property
    def ADU_zp(self):
        if band not in self.data['counts_zp'].keys():
            self.Calc_zp(band)
        if band == '*':
            return self.data['counts_zp']
        else:
            return self.data['counts_zp'][band]
            
    @property
    def flux_sky(self,band='*'):
        if band not in self.data['flux_sky'].keys():
            self.Calc_m5(band)
        if band == '*':
            return self.data['flux_sky']
        else:
            return self.data['flux_sky'][band]

    def Calc_Inputs(self,band='*'):
        
        list_band=[band]
        if band == '*':
            list_band=self.filterlist

        for bd in list_band:
            myup=self.Calc_Integ_Sed(self.darksky,self.system[bd])
            self.data['Tb'][bd]=self.Calc_Integ(self.atmosphere[bd])
            self.data['Sigmab'][bd]=self.Calc_Integ(self.system[bd])
            self.data['mag_sky'][bd]=-2.5*np.log10(myup/(3631.*self.data['Sigmab'][bd]))

    def Calc_m5(self,filtre='*'):

        list_filters=[filtre]
        if filtre=='*':
            list_filters=self.filterlist

        for filt in list_filters:
            filter_trans=self.system[filt]
            wavelen_min, wavelen_max, wavelen_step=filter_trans.getWavelenLimits(None,None,None)
                    
            bandpass=Bandpass(wavelen=filter_trans.wavelen, sb=filter_trans.sb)

            flatSedb = Sed()
            flatSedb.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
            flux0b=np.power(10.,-0.4*self.mag_sky[filt])
            flatSedb.multiplyFluxNorm(flux0b)
            photParams = PhotometricParameters(bandpass=filt)
            norm=photParams.platescale**2/2.*photParams.exptime/photParams.gain
            trans=filter_trans

            if self.atmos:
                trans=self.atmosphere[filt]
            self.data['m5'][filt]=SignalToNoise.calcM5(flatSedb,trans,filter_trans,photParams=photParams,FWHMeff=self.FWHMeff[filt])
            adu_int= flatSedb.calcADU(bandpass=trans, photParams=photParams)
            self.data['flux_sky'][filt]=adu_int*norm

    def Calc_zp(self,band='*'):
        
        list_band=[band]
        if band=='*':
            list_band=self.filterlist

        for filt in list_band:
            self.zp_filter(filt)

    def zp_filter(self,filtre):

        photParams=PhotometricParameters(bandpass=filtre)
        Diameter=2.*np.sqrt(photParams.effarea*1.e-4/np.pi) # diameter in meter
        Cte=3631.*np.pi*Diameter**2*2.*photParams.exptime/4/h/1.e36

        self.data['Skyb'][filtre]=Cte*np.power(Diameter/6.5,2.)*np.power(2.*photParams.exptime/30.,2.)*np.power(photParams.platescale,2.)*np.power(10.,0.4*(25.-self.mag_sky[filtre]))*self.Sigmab[filtre]
            
        Zb=181.8*np.power(Diameter/6.5,2.)*self.Tb[filtre]
        mbZ=25.+2.5*np.log10(Zb) 
        filtre_trans=self.system[filtre]
        wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
        bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)
        flatSed = Sed()
        flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0=np.power(10.,-0.4*mbZ)
        flatSed.multiplyFluxNorm(flux0)
        photParams=PhotometricParameters(bandpass=filtre)
        counts = flatSed.calcADU(bandpass, photParams=photParams) #number of counts for exptime
        self.data['zp'][filtre]=mbZ
        #print 'hello',counts/self.photParams.exptime
        self.data['counts_zp'][filtre]=counts/2.*photParams.exptime

    def Calc_Integ(self,bandpass):
        resu=0.
        dlam=0
        for i,wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam=bandpass.wavelen[i+1]-wave
                resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

        return resu

    def Calc_Integ_Sed(self,sed,bandpass,wavelen=None, fnu=None):
      
        use_self = sed._checkUseSelf(wavelen, fnu)
        # Use self values if desired, otherwise use values passed to function.
        if use_self:
            # Calculate fnu if required.
            if sed.fnu is None:
                # If fnu not present, calculate. (does not regrid).
                sed.flambdaTofnu()
            wavelen = sed.wavelen
            fnu = sed.fnu
        # Make sure wavelen/fnu are on the same wavelength grid as bandpass.
        wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
    
        # Calculate the number of photons.
        nphoton = (fnu / wavelen * bandpass.sb).sum()
        dlambda = wavelen[1] - wavelen[0]
        return nphoton * dlambda

    def flux_to_mag(self, flux, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        #print 'zp',zp,band
        m = -2.5 * np.log10(flux) + zp
        return m

    def mag_to_flux(self, mag, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        return np.power(10., -0.4 * (mag-zp))


    def zero_points(self, band):
        return np.asarray([self.zp[b] for b in band])

    def mag_to_flux_e_sec(self,mag,band,trans,sed):

        #this should be debugged at some point
        photrams=PhotometricParameters(bandpass=band)
        E_per_sec = sed.calcADU(bandpass=trans, photParams=photParams)
        e_per_sec/=exptime/photParams.gain
        return e_per_sec
        
