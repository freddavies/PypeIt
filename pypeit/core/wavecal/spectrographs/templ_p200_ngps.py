""" Generate the wavelength templates for P200/NGPS"""
import os

from pypeit.core.wavecal import templates

##
# NGPS Red Channel
##

# ##############################
def p200_ngps_R(overwrite=False):  # NGPS_R

    binspec = 2 # KEY PARAMETER, WILL NEED TO REDO MANUAL WAVELENGTH SOLUTION WHEN BINNING IS CHANGED
    outroot = 'p200_ngps_R.fits'
    
    # Assume for right now that wavelength solution is roughly equivalent for all 3 slits
    slits = [98, 254, 424] 
    lcut = None

    wpath = os.path.join(templates.template_path, 'P200_NGPS', 'wvcalib_R.fits')
    basefiles = ['wvarxiv_p200_ngps_20250117T1639.fits', 'wvarxiv_p200_ngps_20250117T1639.fits', 'wvarxiv_p200_ngps_20250117T1639.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    templates.build_template(wfiles, slits, lcut, binspec, outroot,
        lowredux=False, normalize=True, overwrite=overwrite)

if __name__ == '__main__':
     p200_ngps_R(overwrite=True)