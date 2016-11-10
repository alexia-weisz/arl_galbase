import numpy as np
import astropy.io.fits as pyfits
import config
from pdb import set_trace
import gal_data
import extract_stamp
import warnings

_OUT_DIR = '../cutouts/sings/'
_MIPS_DIR = '/data/tycho/0/leroy.42/ellohess/data/mips/sings/'
_LOS_DIR = '../../ellohess/code/index/'
_KERNEL_DIR = '/data/tycho/0/leroy.42/ellohess/kernels/Low_Resolution/'
_TEST_WRITE_DIR = '/n/home00/lewis.1590/research/galbase_allsky/cutouts/'

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='Create cutouts of a given size around each galaxy center.')
    parser.add_argument('--size', default=30, help='cutout size in arcminutes. Default: 30.')
    parser.add_argument('--cutout', action='store_true')
    parser.add_argument('--copy', action='store_true')
    parser.add_argument('--convolve', action='store_true')
    parser.add_argument('--align', action='store_true')
    parser.add_argument('--write_info', action='store_true', help='write galaxy info and output to file, e.g., # of tiles, time to completion.')
    return parser.parse_args()


def main(**kwargs):

    if kwargs['cutout']:
        warnings.filterwarnings('ignore')

        gals = gal_data.gal_data(tag='SINGS')
        n_gals = len(gals)
        size_deg = kwargs['size'] * 60. / 3600.

        for i in range(n_gals):
            this_gal = np.rec.fromarrays(gals[i], names=list(config.COLUMNS))
            galname = str(this_gal.name).replace(' ', '').upper()
            #print galname
            extract_stamp.galex(band='fuv', ra_ctr=this_gal.ra_deg, dec_ctr=this_gal.dec_deg, size_deg=size_deg, name=galname, write_info=kwargs['write_info'])

            #extract_stamp.galex(band='nuv', ra_ctr=this_gal.ra_deg, dec_ctr=this_gal.dec_deg, size_deg=size_deg, name=galname)

    if kwargs['copy']:
        pass

    if kwargs['convolve']:
        pass

    if kwargs['align']:
        pass



if __name__ == '__main__':
    args = get_args()
    main(**vars(args))
