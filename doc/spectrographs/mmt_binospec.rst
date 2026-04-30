************
MMT Binospec
************

Overview
========

This file summarizes several instrument specific
items for the MMTO's Binospec spectrograph.

Wavelength Calibration
++++++++++++++++++++++

Templates were created in 2023 that cover the usable range for each of the 3 gratings. Only the bluest part of the
1000 l/mm grating's range isn't fully covered by the template, but there should be sufficient overlap for the template
to still work at the bluest settings, e.g. a central wavelength of 4000 A. Better test data is needed to confirm this.

Bad pixel mask
++++++++++++++

The static bad pixel mask is loaded from FITS files derived from the IDL
pipeline calibration data (``badpix_binospec.fits`` plus the hard-coded bad
columns and detector trap regions defined in ``bino_mosaic.pro``).  This
adds roughly 12,500 individual bad pixels per detector, along with bad
columns and detector trap region masking, replacing the small set of
hard-coded bad columns previously identified from 2019 flat and bias
observations.  The mask files are distributed with the package as
``static_calibs/mmt_binospec/bpm_binospec_det{1,2}.fits.gz``.

