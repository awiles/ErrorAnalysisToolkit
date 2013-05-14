ErrorAnalysisToolkit
====================

Tools to examine image registration and tracking errors.

This "toolkit" started as a collection of code from my PhD, but there is some other recent example code for NDI Aurora communication set up in matlab.  It is not completely up to date, but eventually this tool kit will be able to comletely replicate my code simulations from my thesis and published papers.  This is clearly long overdue but it is still in progress.

matlab/communication
--------------------

This folder contains sample code for communicating with the NDI Aurora system using Matlab.

matlab/math
-----------

This folder contains a collection of math routines for linear transformations, rotations using matrices, Euler angles and quaternions.  There is also some registration routines used in the target_error simulations.

matlab/plotting
---------------

This includes a collection of useful plotting routines to display the results.

matlab/stats
------------

Collection of statistical methods used in the results analysis.

matlab/target_error
-------------------

The target error models and simulations contain three main areas of simulation.

1. Estimation Target Registration Error (TRE) for anisotropic and non-homogenous Fiducial Localizer Error (FLE)

2. Estimation of the Fiducial Localizer Error given an estimate of the recent Fiducial Registration Error

3. Estimation of the Target Tracking Error given an estimate of the uncertainty of the transformation parameters.

All of these experiments are documented in my thesis and some of these results are included in various published papers:

A. Danilchenko, A. D. Wiles, R. Balachandran, and J. M. Fitzpatrick, “Improved method for point-based tracking,” in Medical Image Computing and Computer-Assisted Intervention MICCAI 2010 (T. Jiang, N. Navab, J. P. Pluim, and M. A. Viergever, eds.), vol. 6363 of Lecture Notes in Computer Science, (Beijing, China), September 2010.

A. D. Wiles and T. M. Peters, “Target tracking errors for 5D and 6D spatial measurement systems,” IEEE Transactions on Medical Imaging, vol. 29, pp. 879 – 894, March 2010.

A. D. Wiles and T. M. Peters, “Real-time estimation of FLE statistics for 3D tracking with point-based registration,” IEEE Transactions on Medical Imaging, vol. 28, pp. 1384–1398, September 2009.

A. D. Wiles, A. Likholyot, D. D. Frantz, and T. M. Peters, “A statistical model for point-based target registration error with anisotropic fiducial localizer error,” IEEE Transactions on Medical Imaging, vol. 27, pp. 378–390, March 2008.

A. D. Wiles and T. M. Peters, “Real-time estimation of FLE for point-based registration,” in Proceedings of SPIE, Medical Imaging 2009: Visualization, Image-Guided Procedures, and Modeling (M. I. Miga and K. H.Wong, eds.), vol. 7261, (Lake Buena Vista, FL, USA), p. 72613B, SPIE, February 2009.

A. D. Wiles and T. M. Peters, “Improved statistical TRE model when using a reference frame.,” in Medical Image Computing and Computer-Assisted Intervention MICCAI 2007 (N. Ayache, S. Ourselin, and A. Maeder, eds.), vol. 4791 of Lecture Notes in Computer Science, (Brisbane, QLD, Australia), pp. 442–449, October/November 2007.