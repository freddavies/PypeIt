import numpy as np

from pypeit import msgs
from pypeit.calibframe import CalibFrame
from pypeit.images import buildimage
from pypeit import spec2dobj

def get_sci_metadata(spectrograph, fitstbl, frame, det):
    """
    Grab the meta data for a given science frame and specific detector

    Args:
        frame (int): Frame index
        det (int): Detector index

    Returns:
        5 objects are returned::
            - str: Object type;  science or standard
            - str: Setup/configuration string
            - `astropy.time.Time`_: Time of observation
            - str: Basename of the frame
            - str: Binning of the detector

    """

    # Set binning, obstime, basename, and objtype
    binning = fitstbl['binning'][frame]
    obstime  = fitstbl.construct_obstime(frame)
    basename = fitstbl.construct_basename(frame, obstime=obstime)
    types  = fitstbl['frametype'][frame].split(',')
    if 'science' in types:
        objtype_out = 'science'
    elif 'standard' in types:
        objtype_out = 'standard'
    else:
        msgs.error('get_sci_metadata() should only be run on standard or science frames.  '
                    f'Types of this frame are: {types}')
    calib_key = CalibFrame.construct_calib_key(fitstbl['setup'][frame],
                                                fitstbl['calib'][frame],
                                                spectrograph.get_det_name(det))
    return objtype_out, calib_key, obstime, basename, binning


def process_one(spectrograph, fitstbl, caliBrate, par, frames, 
                det, bg_frames=None, std_outfile=None, science_path:str=None):
    """
    Reduce Objects in a single exposure/detector pair

    sci_ID and det need to have been set internally prior to calling this method

    Parameters
    ----------
    frames : :obj:`list`
        List of frames to extract; stacked if more than one is provided
    det : :obj:`int`
        Detector number (1-indexed)
    bg_frames : :obj:`list`, optional
        List of frames to use as the background. Can be empty.
    std_outfile : :obj:`str`, optional
        Filename for the standard star spec1d file. Passed directly to
        :func:`~pypeit.specobjs.get_std_trace`.

    Returns
    -------
    global_sky : `numpy.ndarray`_
        Initial global sky model
    sobjs_obj : :class:`~pypeit.specobjs.SpecObjs`
        List of objects found
    sciImg : :class:`~pypeit.images.pypeitimage.PypeItImage`
        Science image
    bkg_redux_sciimg : :class:`~pypeit.images.pypeitimage.PypeItImage`
        Science image before background subtraction
        if self.bkg_redux is True, otherwise None.
        It's used to generate a global sky model without bkg subtraction.
    objFind : :class:`~pypeit.find_objects.FindObjects`
        Object finding speobject

    """
    # Grab some meta-data needed for the reduction from the fitstbl
    objtype, setup, obstime, basename, binning \
            = get_sci_metadata(spectrograph, fitstbl, frames[0], det)

    msgs.info("Image processing begins for {} on det={}".format(basename, det))

    # Is this a standard star?
    std_redux = objtype == 'standard'
    frame_par = par['calibrations']['standardframe'] \
                    if std_redux else par['scienceframe']

    # Build Science image
    sci_files = fitstbl.frame_paths(frames)
    sciImg = buildimage.buildimage_fromlist(
        spectrograph, det, frame_par,
        sci_files, bias=caliBrate.msbias, bpm=caliBrate.msbpm,
        dark=caliBrate.msdark,
        scattlight=caliBrate.msscattlight,
        flatimages=caliBrate.flatimages,
        slits=caliBrate.slits,  # For flexure correction
        ignore_saturation=False)

    # get no bkg subtracted sciImg to generate a global sky model without bkg subtraction.
    # it's a dictionary with only `image` and `ivar` keys if bkg_redux=False, otherwise it's None
    bkg_redux_sciimg = None

    # Background Image?
    if bg_frames is not None and len(bg_frames) > 0:
        # get no bkg subtracted sciImg
        bkg_redux_sciimg = sciImg
        # Build the background image
        bg_file_list = fitstbl.frame_paths(bg_frames)
        # TODO I think we should create a separate self.par['bkgframe'] parameter set to hold the image
        # processing parameters for the background frames.  This would allow the user to specify different
        # parameters for the background frames than for the science frames.  
        bgimg = buildimage.buildimage_fromlist(spectrograph, det, frame_par, bg_file_list,
                                                bpm=caliBrate.msbpm,
                                                bias=caliBrate.msbias,
                                                dark=caliBrate.msdark,
                                                scattlight=caliBrate.msscattlight,
                                                flatimages=caliBrate.flatimages,
                                                slits=caliBrate.slits,
                                                ignore_saturation=False)

        # NOTE: If the spatial flexure exists for sciImg, the subtraction
        # function propagates that to the subtracted image, ignoring any
        # spatial flexure determined for the background image.
        sciImg = bkg_redux_sciimg.sub(bgimg)

    '''
    # Save intermediate files
    if science_path is not None:
        # sciImg
        outfile_sci = science_path / f'sciImg_{basename}.fits'
        sciImg.to_file(outfile_sci, overwrite=True)
        if bkg_redux_sciimg is not None:
            # bkg_redux_sciimg
            outfile_bkg = science_path / f'bkg_{basename}.fits'
            bkg_redux_sciimg.to_file(outfile_bkg, overwrite=True)
    '''

    # Flexure
    spat_flexure = None
    # use the flexure correction in the "shift" column
    manual_flexure = fitstbl[frames[0]]['shift']
    if (objtype == 'science' and par['scienceframe']['process']['spat_flexure_correct']) or \
            (objtype == 'standard' and par['calibrations']['standardframe']['process']['spat_flexure_correct']) or \
                manual_flexure:
        if (manual_flexure or manual_flexure == 0) and not (np.issubdtype(fitstbl[frames[0]]["shift"], np.integer)):
            msgs.info(f'Implementing manual flexure of {manual_flexure}')
            spat_flexure = np.float64(manual_flexure)
            sciImg.spat_flexure = spat_flexure
        else:
            msgs.info(f'Using auto-computed flexure')
            spat_flexure = sciImg.spat_flexure
    #msgs.info(f'Flexure being used is: {spat_flexure}')


    # Return
    return sciImg, bkg_redux_sciimg
