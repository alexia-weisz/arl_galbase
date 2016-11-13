import astropy.io.fits as pyfits
import os
import numpy as np
from pdb import set_trace
import config


def empty_gal_struct(n):
    init_vals = config.INIT_VALS
    columns = config.COLUMNS
    ts = config.COL_TYPES

    dtype = (zip(columns, ts))
    empty = np.recarray((n,),
                        dtype=[('name', object),
                        ('pgc', int),
                        ('alias', object),
                        ('tags', object),
                        ('dist_mpc', float),
                        ('e_dist', float),
                        ('ref_dist', object),
                        ('vhel_kms', float),
                        ('e_vhel_kms', float),
                        ('ref_vhel', object),
                        ('vrad_kms', float),
                        ('e_vrad_kms', float),
                        ('ref_vrad', object),
                        ('vvir_kms', float),
                        ('e_vvir_kms', float),
                        ('ref_vvir', object),
                        ('ra_deg', float),
                        ('dec_deg', float),
                        ('ref_pos', object),
                        ('posang_deg', float),
                        ('e_posang', float),
                        ('ref_posang', object),
                        ('incl_deg', float),
                        ('e_incl', float),
                        ('ref_incl', object),
                        ('log_raxis', float),
                        ('e_log_raxis', float),
                        ('ref_log_raxis', object),
                        ('t', float),
                        ('e_t', float),
                        ('ref_t', object),
                        ('morph', object),
                        ('bar', int),
                        ('ring', int),
                        ('multiple', int),
                        ('ref_morph', object),
                        ('av_sf11', float),
                        ('r25_deg', float),
                        ('e_r25_deg', float),
                        ('ref_r25', object),
                        ('vmaxg_kms', float),
                        ('e_vmaxg_kms', float),
                        ('ref_vmaxg', object),
                        ('vrot_kms', float),
                        ('e_vrot_kms', float),
                        ('ref_vrot', object),
                        ('hi_msun', float),
                        ('ref_hi', object),
                        ('lfir_lsun', float),
                        ('ref_ir', object),
                        ('btc_mag', float),
                        ('ref_btc', object),
                        ('ubtc_mag', float),
                        ('ref_ubtc', object),
                        ('bvtc_mag', float),
                        ('ref_bvtc', object),
                        ('itc_mag', float),
                        ('ref_itc', object)] )

    for i in range(len(empty)):
        for j in range(len(empty[i])):
            empty[i][j] = init_vals[j]

    return empty


def gal_data(name=None, data=None, all=False, data_dir=None, found=None, tag=None):

    if not name and not all and not tag:
        print('Need a name to find a galaxy. Returning empty structure')
        return None


    if not data_dir:
        data_dira = config._GALDATA_DIR
        data_dir = config._GALBASE_DIR


    # READ IN THE DATA
    if data is None:
        dbfile = os.path.join(data_dir, 'gal_base.fits')
        hdulist = pyfits.open(dbfile)
        data = hdulist[1].data
        hdulist.close()


    # ALL DATA ARE DESIRED
    if all:
        return data

    # A SPECIFIC SURVEY IS USED
    if tag is not None:
        n_data = len(data)
        keep = np.ones(n_data)
        survey_file = os.path.join(data_dira, 'survey_' + tag.lower() + '.txt')
        gals = np.genfromtxt(survey_file, dtype='string')

        for i in range(n_data):
            this_tag = data['tags'][i].strip(';').split(';;')
            keep[i] = sum(np.in1d(tag, this_tag))

        if np.sum(keep) == 0:
            print('No targets found with that tag combination.')
            return None

        good_data = data[np.where(keep)]

        return data[np.where(keep)]


    # NAME OR LIST OF NAMES
    alias_vec, name_vec = np.loadtxt(os.path.join(data_dir, 'gal_base_alias.txt'), dtype='string', unpack=True)
    alias_vec = [a.replace(' ', '').upper() for a in alias_vec]
    name_vec = [a.replace(' ', '').upper() for a in name_vec]


    # IDENTIFY THE GALAXY
    name = [name]
    n_names = len(name)
    output = empty_gal_struct(n_names)
    found = np.ones(n_names)

    name = np.asarray([a.replace(' ', '').upper() for a in name])
    name_vec = np.asarray([a.replace(' ', '').upper() for a in name_vec])
    alias_vec = np.asarray([a.replace(' ', '').upper() for a in alias_vec])
    data_name = np.asarray([a.replace(' ', '').upper() for a in data.field('name')])

    for i in range(n_names):
        ind = np.where(alias_vec == name[i])[0]

        if len(ind) == 0:
            print('No match for ' + name[i])
            found[i] = 0
            continue

        data_ind = np.where(data_name == name_vec[ind][0])[0]
        output[i] = data[data_ind]
        found[i] = 1

    return output


