# synctexnalize

This Matlab script extracts and interpolates pole figures from 2D detector data obtained by synchrotron diffraction experiments. Pole figures are extracted by single peak fitting in user defined intervals (either trapezoidal numerical integration or parametric integration of a fitted Lorentz function).
This script also includes methods for masking, correcting detector artifacts and background corrections.
A minimum dataset needs 2D detector images at known sample rotations. In order to select a hkl pogrammatically, also wavelength and detector-sample distances need to be known.
The script uses functionalities provided by the mtex toolbox https://mtex-toolbox.github.io
