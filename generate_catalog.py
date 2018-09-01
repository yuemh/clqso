import numpy as np
import os
from astropy.io import fits

def name_to_coord(name):
    J_index = name.find('J')
    start_index = J_index + 1

    mid_index = np.max([name.find('+'), name.find('-')])

    ra_str = name[start_index:mid_index]
    dec_str = name[mid_index:]

    ra_str = ra_str[:2] + ':'\
            + ra_str[2:4] + ':'\
            +ra_str[4:]
    dec_str = dec_str[:3] + ':'\
            + dec_str[3:5] + ':'\
            +dec_str[5:]

    return (ra_str, dec_str)


def coord_to_name(coord, style='full', header=''):
    dummy = 1

    ra_str, dec_str = coord
    new_ra_str = ra_str.replace(':','')
    new_dec_str = dec_str.replace(':','')

    if not new_dec_str[0] in '+-':
        new_dec_str = '+' + new_dec_str

    if style=='full':
        name = header + 'J' + new_ra_str + new_dec_str

    elif style=='abbr':
        name = header + 'J' + new_ra_str[:4] + new_dec_str[:5]

    else:
        print('Unknown style in coord_to_name(). Use style==full')
        name = header + 'J' + new_ra_str + new_dec_str

    return name


def main():
    data = fits.open('./decals-DR5-CLQSOs.fits')[1].data

    namelist = []
    ralist = []
    declist = []
    pmralist = []
    pmdeclist = []
    maglist = []
    prolist = []
    notelist = []

    for index in range(len(data)):
        origname = data['name'][index]
        origname = 'J' + origname
        coord = name_to_coord(origname)
        print(coord)
        newname = coord_to_name(coord, style='abbr')

        namelist.append(newname)
        namelist.append(newname + '_OFF')

        ralist.append(coord[0])
        ralist.append('')

        declist.append(coord[1])
        declist.append('')

        pmralist.append(0.0)
        pmdeclist.append(0.0)
        pmralist.append(0.0)
        pmdeclist.append(0.0)

        maglist.append(data['DECAM_MAG'][index][0])
        maglist.append(0.0)

        prolist.append(1)
        prolist.append(0)

        notelist.append('J2000.0')
        notelist.append('J2000.0')


    alldat = np.array([namelist, ralist, declist, pmralist,\
                       pmdeclist, maglist, prolist, notelist]).T

    np.savetxt('./mmtcat_myue_20180901.dat', alldat,\
               fmt=['%-15s', '%-11s', '%-11s', '%-4s','%-4s',\
                   '%-8s','%-3s','%-8s'])

main()

