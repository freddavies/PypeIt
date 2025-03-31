"""
Module for P200/NGPS specific methods.

.. include:: ../include/links.rst
"""
from typing import List, Optional

import numpy as np

from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from astropy.time import Time

from pypeit import msgs
from pypeit import io
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


def flip_fits_slice(s: str) -> str:
    return '[' + ','.join(s.strip('[]').split(',')[::-1]) + ']'


class P200NGPSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle P200/NGPS specific code
    """
    ndet = 2 # Two Detectors (R and I)
    telescope = telescopes.P200TelescopePar()

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='TELRA', required_ftypes=['science', 'standard']) 
        self.meta['dec'] = dict(ext=0, card='TELDEC', required_ftypes=['science', 'standard']) 
        self.meta['target'] = dict(ext=0, card='NAME', compound=True, required_ftypes=['science', 'standard']) 

        self.meta['dispname'] = dict(card=None, compound=True, default='VPH')
        self.meta['decker'] = dict(ext=0, card='SLITW', rtol=1e-2) 
        self.meta['binning'] = dict(card=None, compound=True) # Got by compound meta function

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='SHUTTIME') # Updated: SHUTTIME more accurate than exptime
        self.meta['airmass'] = dict(ext=0, card='AIRMASS', required_ftypes=['science', 'standard'])

        # Extras for config and frametyping
        self.meta['dichroic'] = dict(card=None, compound=True) # Temporary
        self.meta['dispangle'] = dict(card=None, rtol=1e-2, compound=True) # Got by compound meta function
        self.meta['slitwid'] = dict(ext=0, card='SLITW', rtol=1e-2) #  Slit widths (0.5 and 1 in initial sample)
        self.meta['idname'] = dict(ext=0, card='IMGTYPE') 
        self.meta['instrument'] = dict(ext=0, card='INSTRUME') 

        # Detector (R or I) 
        self.meta['mode'] = dict(ext=0, card='SPEC_ID', compound = True)
        
        # Lamps
        self.meta['lampstat01'] = dict(ext=0, card='LAMPBLUC') # Blue Xe
        self.meta['lampstat02'] = dict(ext=0, card='LAMPFEAR') # FeAr
        self.meta['lampstat03'] = dict(ext=0, card='LAMPREDC') # Red Continuum
        self.meta['lampstat04'] = dict(ext=0, card='LAMPTHAR') # ThAr

    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'mjd':
            return Time(headarr[0]['UTSHUT']).mjd
        elif meta_key == 'dispangle':
            try:
                return 9 # CHANGE
            except Exception as e:
                msgs.warn("Could not read dispangle from header:" + msgs.newline() + str(headarr[0]['ANGLE']))
                raise e    
        
        elif meta_key == 'binning':      
            # Always the same binning for DET01 and DET02    
            binspat = headarr[1]['BINSPAT'] 
            binspec = headarr[1]['BINSPEC']
            return f"{binspec},{binspat}"
        
        elif meta_key == 'target':
            if headarr[0]['TARGET'] is not None:
                return headarr[0]['TARGET']
            else:
                return headarr[0]['IMGTYPE']
            

        return None

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object
        """
        return ['binning', 'slitwid']
    
    # ['dispname', 'binning', 'dispangle', 'dichroic'] # CHANGE

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['GRATING', 'ANGLE', 'APERTURE']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys()
    
    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'SCI')
        
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'DOMEFLAT')
        
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)

        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['idname'] == 'FEAR') | (fitstbl['idname'] == 'THAR'))  ####################################

        
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    def get_rawimage(self, raw_file, det):
        """
        Read raw spectrograph image files and return data and relevant metadata
        needed for image processing.

        For P200/NGPS, the ``DATASEC`` and ``OSCANSEC`` regions are read
        directly from the file header and are automatically adjusted to account
        for the on-chip binning.  This is a simple wrapper for
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage` that
        sets ``sec_includes_binning`` to True.  See the base-class function for
        the detailed descriptions of the input parameters and returned objects.
        """

        return super().get_rawimage(raw_file, det, sec_includes_binning=True)
 
class P200NGPSSpectrograph(P200NGPSSpectrograph):
    """
    Child to handle P200/NGPS specific code
    """
    name = 'p200_ngps'
    camera = 'NGPS'
    header_name = 'NGPS'
    supported = True
    comment = 'R and I channels'
    
    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        # Handle dispangle and mjd from superclass method
        retval = super().compound_meta(headarr, meta_key)
        
        # If superclass could not handle the meta key
        if retval is not None:
            return retval
            
        if meta_key == 'binning':
            # Always the same binning for DET01 and DET02    
            binspat = headarr[1]['BINSPAT'] 
            binspec = headarr[1]['BINSPEC']
            return f"{binspec},{binspat}"
        
        # If there is no target keyword, return image type
        elif meta_key == 'target':
            if headarr[0]['TARGET'] is not None:
                return headarr[0]['TARGET']
            else:
                return headarr[0]['IMGTYPE']
        elif meta_key == 'dichroic': 
            return None
        else:
            msgs.error("Not ready for this compound meta: ", meta_key)


    def get_detector_par(self, det: int, hdu: Optional[fits.HDUList] = None):
        """
        Return metadata for the selected detector.
    
        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
    
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning') 

        # Detector 1 (R Channel)
        detector_dict1 = dict(
            binning         = binning,
            det             = 1, # All R channel images assigned to extension 1
            dataext         = 1, # All R channel images assigned to extension 1
            specaxis        = 1,
            specflip        = False, 
            spatflip        = False, 
            platescale      = 0.5, 
            darkcurr        = 0.0,  # e-/pixel/hour (No dark current)
            saturation      = 45000., # ???
            nonlinear       = 40./45.,
            mincounts       = -1e10, # check
            numamplifiers   = 1, 
            gain            = np.atleast_1d(2.8),
            ronoise         = np.atleast_1d(8.5),
            datasec         = np.atleast_1d(flip_fits_slice(hdu[1].header['DATASEC'])),
            oscansec        = np.atleast_1d(flip_fits_slice(hdu[1].header['BIASSEC']))
        )

        # Detector 2 (I Channel)
        detector_dict2 = dict(
            binning         = binning,
            det             = 2, # All I channel images assigned to extension 2
            dataext         = 2, # All I channel images assigned to extension 2
            specaxis        = 1,
            specflip        = False, 
            spatflip        = False, 
            platescale      = 0.5, 
            darkcurr        = 0.0,  # e-/pixel/hour (No dark current)
            saturation      = 45000., # ???
            nonlinear       = 40./45.,
            mincounts       = -1e10, # check
            numamplifiers   = 1, # Updated
            gain            = np.atleast_1d(2.8),
            ronoise         = np.atleast_1d(8.5),
            datasec         = np.atleast_1d(flip_fits_slice(hdu[2].header['DATASEC'])),
            oscansec        = np.atleast_1d(flip_fits_slice(hdu[2].header['BIASSEC']))
        )

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2]
        detectors = detector_container.DetectorContainer(**detector_dicts[det-1])
        return detectors


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 10. # Lower edge tracing thresdhold to catch leftmost slit
       
        # Try to Remove any false edge detections that cross pixels in stray light regions between slits
        #par['calibrations']['slitedges']['rm_slits'] = ['1:1389:332'] # Remove individual slits ['1:1389:332', '1:2092:172']
        #par['calibrations']['slitedges']['exclude_regions'] = ['1:334:349', '1:1392:170'] # Exclude regions between slits from edge tracing altogether

        par['calibrations']['slitedges']['minimum_slit_length'] = 100 # Set minimum slit length 

        # THIS WORKED FOR INITIAL DATASET
        #par['calibrations']['slitedges']['add_slits'] = ['1:2090:35:161', '2:2090:123:250'] # Add slit to R and I channel where the edge tracing failed (leftmost slit in both cases)
        # THIS WORKED FOR 20250201 DATASET
        #par['calibrations']['slitedges']['add_slits'] = ['1:2090:24:153', '2:2090:123:250']

        # Should always work but lose edges of slit 
        par['calibrations']['slitedges']['add_slits'] = ['1:2090:35:153', '2:2090:123:250']
        

        par['scienceframe']['process']['combine'] = 'median'
        par['calibrations']['standardframe']['process']['combine'] = 'median'

        par['scienceframe']['process']['use_overscan'] = True
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5 

        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'

        par['calibrations']['wavelengths']['lamps'] = ['ThAr'] # FeAr and ThAR lamps for NGPS (ThAr)
        par['calibrations']['wavelengths']['method'] = 'full_template' # Use existing wavelength templates (R and I channels)

        par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_p200_ngps_20250131T1227.fits' # R CHANNEL
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_p200_ngps_20250131T1354.fits' # I CHANNEL


        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 1.0 # For now, use high RMS threshold to ensure a wavelength solution is generated 


        # Do not flux calibrate
        par['fluxcalib'] = None

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        #par['sensfunc']['algorithm'] = 'UVIS'
        #par['sensfunc']['UVIS']['polycorrect'] = False
        #par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'
        return par
    
    def get_wavelength_calib(self, det):
        """
        Return the appropriate wavelength template file based on the detector/channel.
        """
        if det == 1:
            return "wvarxiv_p200_ngps_20250131T1227.fits"  # R channel
        elif det == 2:
            return "wvarxiv_p200_ngps_20250131T1354.fits"  # I channel
 
    
    def config_specific_par(self, scifile, inp_par=None): 
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        #channel = self.get_meta_value(self.get_headarr(hdu), 'binning') 

        headarr = self.get_headarr(scifile)
        channel = self.get_meta_value(headarr, 'mode')
        print(channel)

        if channel == 'R':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_p200_ngps_20250131T1227.fits'
        elif channel == 'I':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_p200_ngps_20250131T1354.fits'

        return par


    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed calibration frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """

        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        msgs.info("Using hard-coded BPM for NGPS") # Can create actual BPM later if necessary
        bpm_img[:, 0] = 1

        return bpm_img
    

        # List of NGPS Slits:
        # DET01: 98, 254, 424
        # DET02: 187, 341, 510

