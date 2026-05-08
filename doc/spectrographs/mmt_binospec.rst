************
MMT Binospec
************

Overview
========

This file summarizes several instrument specific
items for the MMTO's Binospec spectrograph.

Wavelength Calibration
++++++++++++++++++++++

Templates were created in 2023 that cover the usable range for each of the 3
gratings; the 1000 l/mm template was extended in 2026 to support red central
wavelengths.  All wavelengths are vacuum.  The current template coverage is:

==============  ==========================  ==============
Grating         Template file               Range (Å)
==============  ==========================  ==============
270 l/mm        ``mmt_binospec_270.fits``   3825 – 9216
600 l/mm        ``mmt_binospec_600.fits``   4042 – 10047
1000 l/mm       ``mmt_binospec_1000.fits``  3757 – 7211
==============  ==========================  ==============

For the 1000 l/mm grating, central wavelengths from roughly 4500 Å through
6450 Å are covered with full overlap on both detector sides.  Settings near
the blue end of the grating range (e.g. central wavelength ~4000 Å) extend
slightly below the template's blue limit but should still be usable, since
``full_template`` only requires sufficient overlap to anchor the
cross-correlation; better test data are needed to confirm this corner.

Bad pixel mask
++++++++++++++

The bad pixels were identified from flat and bias observations taken in
2019 and need to be verified.

