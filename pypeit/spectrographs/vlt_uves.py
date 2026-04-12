"""
Module for VLT/UVES

.. include:: ../include/links.rst
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.table import Table

from pypeit import log
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit.par import parset
from pypeit.images.mosaic import Mosaic
from pypeit.core.mosaic import build_image_mosaic_transform


class UVESMosaicLookUp:
    """
    Provides the geometry required to mosaic VLT/UVES data.
    Similar to :class:`~pypeit.spectrographs.gemini_gmos.GeminiGMOSMosaicLookUp`

    """
    # Original
    geometry = {
        # blue
        'MSC01': {'default_shape': (2048, 4096),
                  'det1': {'shift': (0., 0.), 'rotation': 0.},},
        # red -- Note: Dekker et al. (2000) say that the red arm has a mosaic of two 2k x 4k CCDs, with a gap of 0.96mm,
        # corresponding to 64 pixels (15 um pixels), however, RJC trimmed each detector by 3 pixels on each side in the
        # cross-dispersion direction. Therefore, instead of 2048 it's 2042, and instead of 64 pixel gap, it's 70 pixels
        # (3 pixels on each side of the two detectors, plus the 64 pixel gap).
        # Using the fit_mosaic_parameters (on the 564, 580, 760, 860 setups) in the dev-suite, the best fit for the gap is 104 pixels
        'MSC02': {'default_shape': (2042 * 2 + 104.0, 4096),
                  'det1': {'shift': (0., 0.), 'rotation': 0.},
                  'det2': {'shift': (2042.0 + 104.0, 0.0), 'rotation': 0.0}},
    }

class VLTUVESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle VLT/UVES specific code.

    This spectrograph is not yet supported.
    """

    telescope = telescopes.VLTTelescopePar()
    url = 'https://www.eso.org/sci/facilities/paranal/instruments/uves.html'
    header_name = 'UVES'
    pypeline = 'Echelle'
    ech_fixed_format = False
    supported = False
    # TODO before support = True
    # 1. Implement flat fielding - DONE
    # 2. Test on several different setups - DONE
    # 3. Implement PCA extrapolation into the blue

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA',
            required_ftypes=['science', 'standard'])  # Need to convert to : separated
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        self.meta['dispname'] = dict(card=None, compound=True)
        self.meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR TYPE')
        self.meta['arm'] = dict(card=None, compound=True)
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['echangle'] = dict(card=None, default=0.0)  # There is no header card for this, but it is required
        self.meta['xdangle'] = dict(card=None, compound=True, rtol=0.01)  # There is no tolerance, really, because it's the central wavelength.

    def compound_meta(self, headarr, meta_key):
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
        if meta_key == 'binning':
            try:
                binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            except KeyError:
                log.warning("Cannot determine spatial binning from the header. Setting to 1")
                binspatial = 1
            try:
                binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            except KeyError:
                log.warning("Cannot determine spectral binning from the header. Setting to 1")
                binspec = 1
            # Parse the binning information into a string
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'arm' or meta_key == 'dispname':
            if 'HIERARCH ESO TPL NAME' in headarr[0]:
                tplid = headarr[0]['HIERARCH ESO TPL NAME'].lower()
                if 'blue' in tplid:
                    arm = 'BLUE'
                elif 'red' in tplid:
                    arm = ('RED')
                elif 'HIERARCH ESO INS PATH' in headarr[0]:
                    arm = headarr[0]['HIERARCH ESO INS PATH']
                else:
                    arm = 'None'
            return arm
        elif meta_key == 'xdangle':
            if 'HIERARCH ESO INS GRAT1 WLEN' in headarr[0]:
                cwlen = headarr[0]['HIERARCH ESO INS GRAT1 WLEN']
            elif 'HIERARCH ESO INS GRAT2 WLEN' in headarr[0]:
                cwlen = headarr[0]['HIERARCH ESO INS GRAT2 WLEN']
            else:
                cwlen = 'None'
            return cwlen
        else:
            log.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to construct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['xdangle', 'arm', 'binning']

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        Bias and dark frames are considered independent of a configuration,
        but the DATE-OBS keyword is used to assign each to the most-relevant
        configuration frame group. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations`.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'bias': ['binning', 'arm'], 'dark': ['binning', 'arm']}

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
        return ['HIERARCH ESO SEQ ARM']

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
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'OBJECT')
                | (fitstbl['idname'] == 'OBJECT,POINT')
                | (fitstbl['idname'] == 'SCIENCE')
                | (fitstbl['idname'] == 'STD,TELLURIC')
                | (fitstbl['idname'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'STD,FLUX')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'DARK')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'LAMP,FLAT')
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'LAMP,WAVE')

        log.warning('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    # IS THIS NEEDED??
    def vet_assigned_ftypes(self, type_bits, fitstbl):
        """

        NOTE: this function should only be called when running pypeit_setup,
        in order to not overwrite any user-provided frame types.

        This method checks the assigned frame types for consistency.
        For frames that are assigned both the science and standard types,
        this method chooses the one that is most likely, by checking if the
        frames are within 10 arcmin of a listed standard star.

        In addition, for this instrument, if a frame is assigned both a
        pixelflat and slitless_pixflat type, the pixelflat type is removed.
        NOTE: if the same frame is assigned to multiple configurations, this
        method will remove the pixelflat type for all configurations, i.e.,
        it is not possible to use slitless_pixflat type for one calibration group
        and pixelflat for another.

        Args:
            type_bits (`numpy.ndarray`_):
                Array with the frame types assigned to each frame.
            fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
                The class holding the metadata for all the frames.

        Returns:
            `numpy.ndarray`_: The updated frame types.

        """
        type_bits = super().vet_assigned_ftypes(type_bits, fitstbl)

        # where pixelflat is assigned
        pixelflat_idx = fitstbl.type_bitmask.flagged(type_bits, flag='pixelflat')

        # find configurations where both pixelflat and slitless_pixflat are assigned
        pixflat_match = np.zeros(len(fitstbl), dtype=bool)

        for f, frame in enumerate(fitstbl):
            if pixelflat_idx[f]:
                match_config_values = []
                pixflat_match[f] = np.any(match_config_values)

        # remove pixelflat from the type_bits
        type_bits[pixflat_match] = fitstbl.type_bitmask.turn_off(type_bits[pixflat_match], 'pixelflat')

        return type_bits

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        This routine is only defined for echelle spectrographs, and it is
        undefined in the base class.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning.

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        det = self.get_detector_par(1)
        binspectral, binspatial = parse.parse_binning(binning)

        # Assume no significant variation (which is likely true)
        return np.ones_like(order_vec)*det.platescale*binspatial


# default_pypeit_par, config_specific_par and get_detector_par different for each arm??
class VLTUVESBlueSpectrograph(VLTUVESSpectrograph):
    
    name = 'vlt_uves_blue'
    camera = 'VLT_UVES_blue'
    ndet = 1

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS SLIT2 WID')

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Adjustments to parameters for Keck HIRES
        turn_off_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median')
        par.reset_all_processimages_par(**turn_off_on)
        # Right now we are using the overscan and not biases becuase the
        # standards are read with a different read mode and we don't yet have
        # the option to use different sets of biases for different standards,
        # or use the overscan for standards but not for science frames

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        #par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        #par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames on UVES ??
        par['calibrations']['pixelflatframe']['exprng'] = [None, 120]
        par['calibrations']['traceframe']['exprng'] = [None, 120]
        par['calibrations']['illumflatframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [1, 600]
        par['scienceframe']['exprng'] = [30, None]

        # Slit tracing
        par['calibrations']['slitedges']['edge_thresh'] = 8.0
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3
        par['calibrations']['slitedges']['max_nudge'] = 0.
        par['calibrations']['slitedges']['overlap'] = False
        par['calibrations']['slitedges']['dlength_range'] = 0.25
        par['calibrations']['slitedges']['mask_off_detector'] = False

        par['calibrations']['slitedges']['add_missed_orders'] = True
        par['calibrations']['slitedges']['order_width_poly'] = 2
        par['calibrations']['slitedges']['order_gap_poly'] = 3

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5  # [5, 5, 5] + 12*[7] # + [5]

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.1
        par['calibrations']['wavelengths']['sigdetect'] = 4.
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['n_final'] = 4

        # Setup dependent
        # 346
        # par['calibrations']['wavelengths']['n_final'] = [3] + 31*[4] + [3]
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_346_1x1.fits'
        # 390
        # par['calibrations']['wavelengths']['n_final'] = [3] + 38*[4] + [3]
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_390_1x1.fits'
        # 437
        # par['calibrations']['wavelengths']['n_final'] = [3] + 29*[4] + [3]
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_437_1x1.fits'

        par['calibrations']['wavelengths']['match_toler'] = 1.5
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'echelle'
        par['calibrations']['wavelengths']['cc_shift_range'] = (-80.,80.)
        par['calibrations']['wavelengths']['cc_thresh'] = 0.6
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.25
        par['calibrations']['wavelengths']['reid_cont_sub'] = False

        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 6
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 2.0
        par['calibrations']['wavelengths']['ech_separate_2d'] = False
        par['calibrations']['wavelengths']['bad_orders_maxfrac'] = 0.5

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['global_sky_std'] = False
        par['reduce']['skysub']['sky_sigrej'] = 4.0

        # local sky subtraction operates on entire slit
        par['reduce']['extraction']['model_full_slit'] = True
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3, 3]
        # number of objects
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Assume that there is max two object in each order.
        par['reduce']['findobj']['maxnumber_std'] = 1  # Assume that there is only one object in each order.

        # Sensitivity function parameters
        # par['sensfunc']['algorithm'] = 'IR'
        # par['sensfunc']['polyorder'] = 5 #[9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        # par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_10500_R120000.fits'
        # par['sensfunc']['IR']['pix_shift_bounds'] = (-40.0,40.0)
        
        # Telluric parameters
        # HIRES is usually oversampled, so the helio shift can be large
        # par['telluric']['pix_shift_bounds'] = (-40.0,40.0)
        # Similarly, the resolution guess is higher than it should be
        # par['telluric']['resln_frac_bounds'] = (0.25,1.25)

        # Coadding
        par['coadd1d']['wave_method'] = 'log10'

        return par
        
    def get_detector_par(self, det, hdu=None):
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
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')
        binfact = int(binning.split(',')[0]) * int(binning.split(',')[1])
        ronoise = 3.8/np.sqrt(binfact) if hdu is None else hdu[0].header['HIERARCH ESO DET OUT1 RON']
        gain = 0.54 if hdu is None else hdu[0].header['HIERARCH ESO DET OUT1 GAIN']

        # Detector
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.22,
            darkcurr        = 0.0,  # e-/pixel/hour
            saturation      = 65535.,
            nonlinear       = 0.7, # Website says 0.6, but we'll push it a bit
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d([gain]),
            ronoise         = np.atleast_1d([ronoise]),
            datasec         = np.atleast_1d('[:,51:2098]'), # Any changes to this line should also be made in the mosaic geometry at the top of this file.
            oscansec        = np.atleast_1d('[:,4:50]'),
            )

        return detector_container.DetectorContainer(**detector_dict)

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  Or a row from the
                metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        par = super().config_specific_par(inp, inp_par=inp_par)

        bin_spec, bin_spat = parse.parse_binning(self.get_meta_value(inp, 'binning'))

        # slit edges
        # NOTE: With add_missed_orders set to True and order_spat_range set to the
        # default (None), the code will try to add missing orders over the full
        # range of the detector mosaic!
        par['calibrations']['slitedges']['order_spat_range'] = [10., 2080./bin_spat]

        # wavelength
        par['calibrations']['wavelengths']['fwhm'] = 8.0/bin_spec

        # Return
        return par

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        angle_fits_file = 'vlt_uves_blue_angle_fits.fits'
        composite_arc_file = 'vlt_uves_blue_composite_arc.fits'

        return [angle_fits_file, composite_arc_file]

    # @property
    # def norders(self):
    #     """
    #     Number of orders observed for this spectrograph.
    #     """
    #     # 346
    #     # return 33
    #     # 390
    #     return 40
    #     # 437
    #     # return 31
    #
    # @property
    # def order_spat_pos(self):
    #     """
    #     Return the expected spatial position of each echelle order.
    #
    #     The following lines generated the values below:
    #
    #     .. code-block:: python
    #
    #         from pypeit import edgetrace
    #         edges = edgetrace.EdgeTraceSet.from_file('Edges_A_1_DET01.fits.gz')
    #
    #         nrm_edges = edges.edge_fit[edges.nspec//2,:] / edges.nspat
    #         slit_cen = ((nrm_edges + np.roll(nrm_edges,1))/2)[np.arange(nrm_edges.size//2)*2+1]
    #
    #     """
    #     # 346
    #     # self.slits.spat_id/self.slits.nspat
    #     # return np.array([0.01672502, 0.04038521, 0.06437769, 0.08873001, 0.11342802,
    #     #                  0.13848158, 0.16389396, 0.18967436, 0.2158253 , 0.24235753,
    #     #                  0.26927502, 0.29658814, 0.32430333, 0.3524269 , 0.3809682 ,
    #     #                  0.40993333, 0.43933332, 0.46917632, 0.49946995, 0.53022563,
    #     #                  0.56144813, 0.59315476, 0.62534499, 0.65802905, 0.69121394,
    #     #                  0.72492322, 0.75916269, 0.79394031, 0.82926452, 0.86515239,
    #     #                  0.90162054, 0.93867218, 0.97624925]) #, 1.01457846
    #     # 390
    #     return np.array([0.01015914, 0.02863268, 0.04737307, 0.06639743, 0.08571161,
    #                      0.10532204, 0.12523394, 0.14545357, 0.16598767, 0.18684648,
    #                      0.20803361, 0.22955166, 0.25140925, 0.27361704, 0.29618321,
    #                      0.31911634, 0.34241571, 0.36609522, 0.3901668 , 0.41464044,
    #                      0.43951508, 0.46480507, 0.49052863, 0.51669461, 0.54329339,
    #                      0.57034692, 0.59786838, 0.62586725, 0.65435386, 0.68334213,
    #                      0.71284255, 0.74286562, 0.77342552, 0.80453873, 0.83622065,
    #                      0.86848146, 0.90133283, 0.93479389, 0.96886826, 1.00359619])
    #     # 437
    #     # return np.array([0.01619767, 0.04055738, 0.06533852, 0.09053352, 0.11616566,
    #     #                  0.14224589, 0.16877437, 0.19576498, 0.22323035, 0.25118463,
    #     #                  0.27963696, 0.30860075, 0.33808912, 0.36811271, 0.39868486,
    #     #                  0.42982015, 0.46153456, 0.49384601, 0.52676794, 0.56030554,
    #     #                  0.59448245, 0.62931538, 0.66481626, 0.70100701, 0.73791999,
    #     #                  0.77556252, 0.81395274, 0.85309815, 0.89305741, 0.9338257 ,
    #     #                  0.97538197])
    #
    # @property
    # def order_spat_width(self):
    #     """
    #     Return the expected spatial position of each echelle order.
    #
    #     The following lines generated the values below:
    #
    #     .. code-block:: python
    #
    #         import numpy as np
    #         from pypeit import slittrace
    #         slits = slittrace.SlitTraceSet.from_file('Slits_A_0_DET01.fits.gz')
    #
    #         np.median(slits.right_init-slits.left_init, axis=0)/slits.nspat
    #
    #     """
    #     # 346
    #     # return np.array(33*[0.019389049210670473])
    #     # 390
    #     return np.array(40 * [0.01628965847925201])
    #     # 437
    #     # return np.array(31 * [0.020187482348035246])
    #
    # @property
    # def orders(self):
    #     """
    #     Return the order number for each echelle order.
    #     """
    #     # 346
    #     # return np.array([153, 152, 151, 150, 149, 148, 147, 146, 145, 144, 143, 142, 141,
    #     #                  140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128,
    #     #                  127, 126, 125, 124, 123, 122, 121], dtype=int)
    #     # 390
    #     return np.array([142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130,
    #                      129, 128, 127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 117,
    #                      116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 105, 104,
    #                      103])
    #     # 437
    #     # return np.array([124, 123, 122, 121, 120, 119, 118, 117, 116, 115, 114, 113, 112,
    #     #                  111, 110, 109, 108, 107, 106, 105, 104, 103, 102, 101, 100,  99,
    #     #                  98,  97,  96,  95,  94])
    #
    # @property
    # def spec_min_max(self):
    #     """
    #     Return the minimum and maximum spectral pixel expected for the
    #     spectral range of each order.
    #     """
    #     # 346
    #     # spec_max = np.asarray([3000]*32 + [925])#, 2460
    #     # spec_min = np.asarray([635] + [0]*32)
    #     # 390
    #     spec_max = np.asarray([3000]*38 + [920, 2740])
    #     spec_min = np.asarray([1650, 330] + [0]*38)
    #     # 437
    #     # spec_max = np.asarray([3000]*30 + [2260])
    #     # spec_min = np.asarray([1060] + [0]*30)
    #     return np.vstack((spec_min, spec_max))


class VLTUVESRedSpectrograph(VLTUVESSpectrograph):

    name = 'vlt_uves_red'
    camera = 'VLT_UVES_red'
    ndet = 2

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS SLIT3 WID')

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['rdx']['detnum'] = [(1,2)]

        # Adjustments to parameters for Keck HIRES
        turn_off_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median')
        par.reset_all_processimages_par(**turn_off_on)
        # Right now we are using the overscan and not biases becuase the
        # standards are read with a different read mode and we don't yet have
        # the option to use different sets of biases for different standards,
        # or use the overscan for standards but not for science frames

        # Setup dependent -- This is only temporary until we have the reidentification files for all the red settings
        # par['calibrations']['wavelengths']['n_final'] = [3] + cls().norders*[4] + [3]
        # 564l
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_564l_1x1.fits'
        # par['rdx']['detnum'] = 1
        # 564u
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_564u_1x1.fits'
        # par['rdx']['detnum'] = 2
        # 580l
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_580l_1x1.fits'
        # par['rdx']['detnum'] = 1
        # 580u
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_580u_1x1.fits'
        # par['rdx']['detnum'] = 2
        # 760l
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_760l_1x1.fits'
        # par['rdx']['detnum'] = 1
        # 760u
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_760u_1x1.fits'
        # par['rdx']['detnum'] = 2
        # 860l
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_860l_1x1.fits'
        # par['rdx']['detnum'] = 1
        # 860u
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_uves_860u_1x1.fits'
        # par['rdx']['detnum'] = 2

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        #par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        #par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames on UVES ??
        par['calibrations']['pixelflatframe']['exprng'] = [None, 120]
        par['calibrations']['traceframe']['exprng'] = [None, 120]
        par['calibrations']['illumflatframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [1, 600]
        par['scienceframe']['exprng'] = [30, None]

        # Set default processing for slitless_pixflat
        par['calibrations']['slitless_pixflatframe']['process']['scale_to_mean'] = True

        # Slit tracing
        par['calibrations']['slitedges']['edge_thresh'] = 8.0
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3
        par['calibrations']['slitedges']['max_nudge'] = 0.
        par['calibrations']['slitedges']['overlap'] = True
        par['calibrations']['slitedges']['dlength_range'] = 0.25
        par['calibrations']['slitedges']['mask_off_detector'] = True

        par['calibrations']['slitedges']['add_missed_orders'] = True
        par['calibrations']['slitedges']['order_width_poly'] = 4
        par['calibrations']['slitedges']['order_gap_poly'] = 3

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.1
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['n_final'] = 4

        par['calibrations']['wavelengths']['match_toler'] = 1.5
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'echelle'
        par['calibrations']['wavelengths']['cc_shift_range'] = (-80.,80.)
        par['calibrations']['wavelengths']['cc_thresh'] = 0.6
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.25
        par['calibrations']['wavelengths']['reid_cont_sub'] = False

        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 6
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 2.0
        par['calibrations']['wavelengths']['ech_separate_2d'] = True  # Doesn't seem like there's an offset+rotation that works for VLT/UVES
        par['calibrations']['wavelengths']['bad_orders_maxfrac'] = 0.5

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['global_sky_std'] = False
        # local sky subtraction operates on entire slit
        par['reduce']['extraction']['model_full_slit'] = True
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3, 3]
        # number of objects
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Assume that there is max two object in each order.
        par['reduce']['findobj']['maxnumber_std'] = 1  # Assume that there is only one object in each order.

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 5 #[9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_10500_R120000.fits'
        par['sensfunc']['IR']['pix_shift_bounds'] = (-40.0,40.0)

        # Telluric parameters
        # HIRES is usually oversampled, so the helio shift can be large
        par['telluric']['pix_shift_bounds'] = (-40.0,40.0)
        # Similarly, the resolution guess is higher than it should be
        par['telluric']['resln_frac_bounds'] = (0.25,1.25)

        # Coadding
        par['coadd1d']['wave_method'] = 'log10'

        return par

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  Or a row from the
                metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        par = super().config_specific_par(inp, inp_par=inp_par)

        bin_spec, bin_spat = parse.parse_binning(self.get_meta_value(inp, 'binning'))

        # slit edges
        # NOTE: With add_missed_orders set to True and order_spat_range set to the
        # default (None), the code will try to add missing orders over the full
        # range of the detector mosaic!
        # par['calibrations']['slitedges']['order_spat_range'] = [-50/bin_spat, (2042+50)/bin_spat]
        offset = 50.0  # Extra number of pixels to add to the end of the mosaic to allow for missed orders.
        par['calibrations']['slitedges']['order_spat_range'] = [-offset/bin_spat, (2042 * 2 + 104.0 + offset)/bin_spat]

        # wavelength
        par['calibrations']['wavelengths']['fwhm'] = 8.0/bin_spec

        # Return
        return par

    @property
    def allowed_mosaics(self):
        # TODO :: Move this to the parent?
        """
        Return the list of allowed detector mosaics.

        Only red arm on VLT/UVES requires mosaicing.

        Returns:
            :obj:`list`: List of tuples, where each tuple provides the 1-indexed
            detector numbers that can be combined into a mosaic and processed by
            PypeIt.
        """
        return [(1,), (1,2)]

    @property
    def default_mosaic(self):
        return self.allowed_mosaics[1]

    def get_mosaic_par(self, mosaic, hdu=None, msc_ord=0):
        """
        Return the hard-coded parameters needed to construct detector mosaics
        from unbinned images.

        The parameters expect the images to be trimmed and oriented to follow
        the PypeIt shape convention of ``(nspec,nspat)``.  For returned
        lists, the length of the list is the same as the number of detectors in
        the mosaic, and they are ordered by the detector number.

        Args:
            mosaic (:obj:`tuple`):
                Tuple of detector numbers used to construct the mosaic.  Must be
                one among the list of possible mosaics as hard-coded by the
                :func:`allowed_mosaics` function.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent detector parameters are set to a
                default.  BEWARE: If ``hdu`` is not provided, the binning is
                assumed to be `1,1`, which will cause faults if applied to
                binned images!
            msc_ord (:obj:`int`, optional):
                Order of the interpolation used to construct the mosaic.

        Returns:
            :class:`~pypeit.images.mosaic.Mosaic`: Object with the mosaic *and*
            detector parameters.
        """

        # Validate the entered (list of) detector(s)
        nimg, _ = self.validate_det(mosaic)

        # Index of mosaic in list of allowed detector combinations
        mosaic_id = self.allowed_mosaics.index(mosaic)+1
        # is this needed? Since only red arm MSCO2 needs mosaic
        detid = f'MSC0{mosaic_id}'

        # Get the detectors
        detectors = np.array([self.get_detector_par(det, hdu=hdu) for det in mosaic])
        # Binning *must* be consistent for all detectors
        if any(d.binning != detectors[0].binning for d in detectors[1:]):
            log.error('Binning is somehow inconsistent between detectors in the mosaic!')

        # Collect the offsets and rotations for *all unbinned* detectors in the
        # full instrument, ordered by the number of the detector.  Detector
        # numbers must be sequential and 1-indexed.
        # See the mosaic documentation.
        msc_geometry = UVESMosaicLookUp.geometry
        expected_shape = msc_geometry[detid]['default_shape']
        shift = np.array([(msc_geometry[detid]['det1']['shift'][0], msc_geometry[detid]['det1']['shift'][1]),
                          (msc_geometry[detid]['det2']['shift'][0], msc_geometry[detid]['det2']['shift'][1])])

        rotation = np.array([msc_geometry[detid]['det1']['rotation'], msc_geometry[detid]['det2']['rotation']])

        # The binning and process image shape must be the same for all images in
        # the mosaic
        binning = tuple(int(b) for b in detectors[0].binning.split(','))
        shape = tuple(n // b for n, b in zip(expected_shape, binning))

        msc_sft = [None]*nimg
        msc_rot = [None]*nimg
        msc_tfm = [None]*nimg

        for ii in range(nimg):
            msc_sft[ii] = shift[ii]
            msc_rot[ii] = rotation[ii]
            # binning is here in the PypeIt convention of (binspec, binspat), but the mosaic transformations
            # occur in the raw data frame, which flips spectral and spatial
            msc_tfm[ii] = build_image_mosaic_transform(shape, msc_sft[ii], msc_rot[ii], tuple(reversed(binning)))

        return Mosaic(mosaic_id, detectors, shape, np.array(msc_sft), np.array(msc_rot),
                      np.array(msc_tfm), msc_ord)

    def get_detector_par(self, det, hdu=None):
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
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector base parameters
        detector_base = dict(
            binning         = binning,
            specaxis        = 0,
            specflip        = True,
            spatflip        = True,
            platescale      = 0.135,
            darkcurr        = 0.0,  # e-/pixel/hour
            saturation      = 65535.,
            nonlinear       = 0.7, # Website says 0.6, but we'll push it a bit
            mincounts       = -1e10,
            numamplifiers   = 1,
            # Placeholders, will be updated for each detector -- must all be None, so that  the code will break
            # if it tries to obtain information without passing in the hdu argument.
            det=None,
            dataext=None,
            gain=None,
            ronoise=None,
            datasec=None,
            oscansec=None,
        )
        # Now, depending on the HDU format (which changed at some point into multi-extension fits files),
        # we need to extract information from different HDUs for the two detectors.
        detector_dict1 = detector_base.copy()
        detector_dict2 = detector_base.copy()
        if hdu is not None:
            if len(hdu) == 1:
                # This is the old format, where the two detectors are stored in a single HDU
                # Detector 1
                detector_dict1.update(dict(
                    det=1,
                    dataext=0,
                    gain=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT4 GAIN']]),
                    ronoise=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT4 RON']]),
                    datasec=np.atleast_1d('[:,2201:4242]'), # Any changes to this line requires changes to the mosaic at the top of this file.
                    oscansec=np.atleast_1d('[:,4246:4292]'),
                ))

                # Detector 2
                detector_dict2.update(dict(
                    det=2,
                    dataext=0,
                    gain=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT1 GAIN']]),
                    ronoise=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT1 RON']]),
                    datasec=np.atleast_1d('[:,57:2098]'),
                    oscansec=np.atleast_1d('[:,2102:2148]'),
                ))
            else:
                # This is the new format, where the two detectors are stored in separate HDUs.
                # Detector 1.
                detector_dict1.update(dict(
                    det=1,
                    dataext=2,
                    gain=np.atleast_1d([hdu[2].header['HIERARCH ESO DET OUT1 GAIN']]),
                    ronoise=np.atleast_1d([hdu[2].header['HIERARCH ESO DET OUT1 RON']]),
                    datasec=np.atleast_1d('[:,57:2098]'), # Any changes to this line requires changes to the mosaic at the top of this file.
                    oscansec=np.atleast_1d('[:,2102:]'),
                ))

                # Detector 2.
                detector_dict2.update(dict(
                    det=2,
                    dataext=1,
                    gain=np.atleast_1d([hdu[1].header['HIERARCH ESO DET OUT1 GAIN']]),
                    ronoise=np.atleast_1d([hdu[1].header['HIERARCH ESO DET OUT1 RON']]),
                    datasec=np.atleast_1d('[:,57:2098]'),
                    oscansec=np.atleast_1d('[:,2102:]'),
                ))

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2]
        return detector_container.DetectorContainer( **detector_dicts[det-1])

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        angle_fits_file = 'vlt_uves_red_angle_fits.fits'
        composite_arc_file = 'vlt_uves_red_composite_arc.fits'

        return [angle_fits_file, composite_arc_file]

    # @property
    # def norders(self):
    #     """
    #     Number of orders observed for this spectrograph.
    #     """
    #     return 24  # 564l
    #     return 16  # 564u
    #     return 23  # 580l
    #     return 16  # 580u
    #     return 27  # 760l
    #     return 16  # 760u
    #     return 20  # 860l
    #     return 13  # 860u
    #
#     @property
#     def order_spat_pos(self):
#         """
#         Return the expected spatial position of each echelle order.
#
#         The following lines generated the values below:
#
#         .. code-block:: python
#
#             from pypeit import edgetrace
#             import numpy as np
#             edges = edgetrace.EdgeTraceSet.from_file('Edges_H_7_DET02.fits.gz')
#
#             nrm_edges = edges.edge_fit[edges.nspec//2,:] / edges.nspat
#             slit_cen = ((nrm_edges + np.roll(nrm_edges,1))/2)[np.arange(nrm_edges.size//2)*2+1]
#         """
#         # 564l
#         # return np.array([0.01578222, 0.05138392, 0.08686124, 0.1229205 , 0.15956966,
#         #                 0.19682514, 0.23469824, 0.27320327, 0.31235287, 0.35216309,
#         #                 0.39264916, 0.43382898, 0.47571848, 0.51833697, 0.56169515,
#         #                 0.60581347, 0.65070743, 0.69639906, 0.74290928, 0.79025851,
#         #                 0.83847016, 0.88756732, 0.93757559, 0.98655107])
#         # 564u
#         # return np.array([0.04356585, 0.09744982, 0.15237898, 0.20838549, 0.26550201,
#         #                  0.32375955, 0.38319573, 0.4438481 , 0.50575516, 0.56895762,
#         #                  0.63350089, 0.69943186, 0.76679414, 0.8356444 , 0.90603342,
#         #                  0.98124762])
#         # 580l
#         # return np.array([-0.00666975,  0.04932947,  0.10641487,  0.16463922,  0.22403526,
#         #                  0.28464966,  0.34650731,  0.40965558,  0.47413649,  0.53999441,
#         #                  0.60727987,  0.67603839,  0.74632731,  0.81820221,  0.89171804,
#         #                  0.96691536])
#         # 580u
#         # return np.array([-0.0058766,  0.04932947,  0.10641487,  0.16463922,  0.22403526,
#         #                  0.28464966,  0.34650731,  0.40965558,  0.47413649,  0.53999441,
#         #                  0.60727987,  0.67603839,  0.74632731,  0.81820221,  0.89171804,
#         #                  0.96691536])
#         # 760l
#         # return np.array([0.02183773, 0.05053537, 0.07979247, 0.10962473, 0.14005186,
#         #                  0.17108576, 0.20274669, 0.23504974, 0.26801479, 0.30166046,
#         #                  0.33600903, 0.37108002, 0.40689516, 0.44347488, 0.48084019,
#         #                  0.51902689, 0.55804775, 0.59793554, 0.63871922, 0.68042751,
#         #                  0.72309037, 0.76674139, 0.81141437, 0.85714477, 0.90397296,
#         #                  0.95193807, 0.99918404])
#         # 760u
#         # return np.array([-0.00228510, 0.04915505, 0.10209448, 0.15640136, 0.21213053, 0.26934722,
#         #        0.32809903, 0.38845984, 0.45049678, 0.51428256, 0.57989583,
#         #        0.64742068, 0.71694612, 0.78857003, 0.86239197, 0.93852161])#, 1.01707693])
#         # 860l
#         # return np.array([0.02370136, 0.06174571, 0.10269732, 0.14457767, 0.18740961,
#         #                  0.23122347, 0.27605553, 0.32193436, 0.36890337, 0.4170002 ,
#         #                  0.46626546, 0.51674642, 0.56848357, 0.62151777, 0.67590628,
#         #                  0.73170204, 0.7889566 , 0.84773519, 0.90809881, 0.9701096 ])
#         # 860u
#         return np.array([0.04758201, 0.1148355 , 0.18419335, 0.25561299, 0.32914172,
#                          0.40505491, 0.48333775, 0.56381127, 0.64708008, 0.73297388,
#                          0.82163438, 0.9131291 , 1.00757825])
#
#     @property
#     def order_spat_width(self):
#         """
#         Return the expected spatial position of each echelle order.
#
#         The following lines generated the values below:
#
#         .. code-block:: python
#
#             import numpy as np
#             from pypeit import slittrace
#             slits = slittrace.SlitTraceSet.from_file('Slits_B_1_DET01.fits.gz')
#
#             tmp = np.median(slits.right_init-slits.left_init, axis=0)/slits.nspat
#             print(tmp)
#             print(np.median(tmp))
#         """
#         # 564l and 564u
#         # return np.array([0.03]*self.norders)
#         # 580l and 580u
#         # return np.array([0.0325]*self.norders)
#         # 760l and 760u
#         # return np.array([0.023]*self.norders)
#         # 860l and 860u
#         return np.array([0.034]*self.norders)
#
#     @property
#     def orders(self):
#         """
#         Return the order number for each echelle order.
#         """
#         # 564l
#         # return np.array([132, 131, 130, 129, 128, 127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 117, 116, 115, 114, 113, 112, 111, 110, 109], dtype=int)
#         # 564u
#         # return np.array([107, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92], dtype=int)
#         # 580l
#         # return np.array([128, 127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 117, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106], dtype=int)
#         # 580u
#         # return np.array([105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90], dtype=int)
#         # 760l
#         # return np.array([107, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81], dtype=int)
#         # 760u
#         # return np.array([80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65], dtype=int)
#         # 860l
#         # return np.array([91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72], dtype=int)
#         # 860u
#         return np.array([70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58], dtype=int)
#
#     @property
#     def spec_min_max(self):
#         """
#         Return the minimum and maximum spectral pixel expected for the
#         spectral range of each order.
#         """
#         # 564l
#         # spec_max = np.asarray([4096]*23 + [2050])
#         # spec_min = np.asarray([1500] + [0]*23)
#         # 564u
#         # spec_max = np.asarray([4096]*15 + [2600])
#         # spec_min = np.asarray([2750] + [0]*15)
#         # 580l
#         # spec_max = np.asarray([4096]*22 + [2510])
#         # spec_min = np.asarray([2050] + [0]*22)
#         # 580u
#         # spec_max = np.asarray([4096]*15 + [3000])
#         # spec_min = np.asarray([2700] + [0]*15)
#         # 760l
#         # spec_max = np.asarray([4096]*26 + [1180])
#         # spec_min = np.asarray([250] + [0]*26)
#         # 760u
#         # spec_max = np.asarray([4096]*16)
#         # spec_min = np.asarray([2700] + [0]*15)
#         # 860l
#         # spec_max = np.asarray([4096]*19 + [3800])
#         # spec_min = np.asarray([800] + [0]*19)
#         # 860u
#         spec_max = np.asarray([4096]*12 + [1440])
#         spec_min = np.asarray([0]*13)
#         return np.vstack((spec_min, spec_max))

def indexing(itt, postpix, det=None, xbin=1, ybin=1):
    # TODO :: Is this function even needed?
    """
    Some annoying bookkeeping for instrument placement.

    Parameters
    ----------
    itt : int
    postpix : int
    det : int, optional
        Detector number.
    xbin : int, optional
        The binning in the spectral direction.  This is needed to determine the
        size of the unbinned image and thus the location of the postpix.
    ybin : int, optional
        The binning in the spatial direction.  This is needed to determine the
        size of the unbinned image and thus the location of the postpix.

    Returns
    -------

    """
    log.warning("Currently in the indexing function, which is a bit of a mess.  This should be cleaned up before merging.")
    embed()
    assert False
    # Deal with single chip
    if det is not None:
        tt = 0
    else:
        tt = itt
    ii = int(np.round(2048/xbin))  # TODO :: Before merging, the 2048 should be 2048 for the blue chip, but 2042 for the red chip.
    jj = int(np.round(4096/ybin))
    # y indices
    y1, y2 = 0, jj
    o_y1, o_y2 = y1, y2

    # x
    x1, x2 = (tt%4)*ii, (tt%4 + 1)*ii
    if det is None:
        o_x1 = 4*ii + (tt%4)*postpix
    else:
        o_x1 = ii + (tt%4)*postpix
    o_x2 = o_x1 + postpix

    # Return
    return x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2
