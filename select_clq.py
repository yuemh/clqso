import matplotlib.pyplot as plt
from mylib import *

import os,sys,glob
import numpy as np
from astropy.table import Table,hstack,join,vstack
from astropy.io import fits
#from astrotools.DataMgr import GetData

DR_str='DR6'

data_dir = os.getcwd() + '/../data/DECaLS_'+DR_str


def make_qso_table_DR12(SDSS_file, Legacy_file):
    dummy = 1

#    dr12qso = fits.open(SDSS_file)[1].data
#    decals = fits.open(Legacy_file)[1].data

    dr12qso = Table.read(SDSS_file)
    decals = Table.read(Legacy_file)

    name = dr12qso['SDSS_NAME']
    zqso = dr12qso['Z_VI']
    pmjd = dr12qso['PHOTO_MJD']
    psfmag = dr12qso['PSFMAG']
    psferr = dr12qso['ERR_PSFMAG']

    psfmag = psfmag[:,[1,2,4]]
    psferr = psferr[:,[1,2,4]]

    ff = np.vstack([decals['FLUX_%s'%b] for b in 'GRZ'])
    decfluxes = np.ma.masked_array(ff,ff==0)
    decmag = 22.5 - 2.5*np.ma.log10(decfluxes)
    decmag = decmag.transpose()

    dmag = decmag - psfmag # trim Y band

    objid = np.int64(1e5*decals['BRICKID'] + decals['OBJID'])
    apNum = 3

    ff = np.dstack([decals['APFLUX_%s'%b] for b in 'GRZ'])
    apflux = np.ma.masked_array(ff[:,apNum],ff[:,apNum]==0)
    apMag = 22.5 - 2.5*np.ma.log10(apflux)
    t2 = Table({'name':name,'objid':objid,'zqso':zqso,'pmjd':pmjd,\
                'psfMag':psfmag,'psfMagErr':psferr,\
                'DECAM_MAG':decmag,'DECAM_APMAG':apMag,\
                'dmag':dmag},masked=True)

    t = hstack([decals,t2])
    t.write(data_dir+'/decals-%s-QSOs.fits'%DR_str,overwrite=True)



def clqso_sel(qso_file):
    tab = Table.read(qso_file)
    nobs = np.vstack([np.array(tab['NOBS_%s'%b]) for b in 'GRZ']).transpose()
    has_gr = ( np.all(nobs[:,[0,1]] >= 1,axis=1) )# &
    #np.any(tab['DECAM_NOBS'][:,[1,2]] >= 2,axis=1) )
    good_sdss = ~np.all(tab['psfMag']==0,axis=1)
    fading = ( (tab['dmag'][:,0] > 1.5)  & (tab['dmag'][:,1] > 0.) )
    rising = False
    isgal = ~(tab['TYPE'] == 'PSF')
    lowz = tab['zqso'] < 20.0
    bright = tab['DECAM_MAG'][:,1] < 21.5

    indexes = np.where(good_sdss & has_gr & (fading | rising) & isgal & lowz & bright)[0]


    newtab = tab[indexes]
    newtab.write(data_dir+'/decals-%s-CLQSOs.fits'%DR_str,overwrite=True)

    print(np.max(newtab['RA']))
    print(np.min(newtab['RA']))

    plt.hist(newtab['zqso'])
    plt.show()
    plt.hist(newtab['RA'])
    plt.show()

    t2 = Table({'RA':newtab['RA'], 'DEC':newtab['DEC'], 'name':newtab['name']},\
               masked=True)
    t2.write(data_dir+'/decals-%s-CLQSOs-coord.fits'%DR_str,overwrite=True)

    t3 = Table({'name':newtab['name'], 'RA':newtab['RA'], 'DEC':newtab['DEC'],\
                'MG1':newtab['psfMag'][:,1], 'MG2':newtab['DECAM_MAG'][:,0],\
               'MR1':newtab['psfMag'][:,2], 'MR2':newtab['DECAM_MAG'][:,1],\
                'DMG':newtab['dmag'][:,0],'DMR':newtab['dmag'][:,1]},\
               masked=True)
    t3.write(data_dir+'/decals-%s-CLQSOs-simple.fits'%DR_str,overwrite=True)


    return indexes



def main():

    step = 'select_clqso'

    if step == 'combine_catalog':
        make_qso_table_DR12(data_dir+'/DR12Q.fits',\
                            data_dir+'/survey-dr6-dr12Q.fits')

    elif step == 'select_clqso':
        selected_index = clqso_sel(data_dir+'/decals-%s-QSOs.fits'%DR_str)
        print(selected_index)
        print(len(selected_index))


if __name__=='__main__':
    main()
