When testing the reduction, we found there were a few patches we needed to do to KSkyWizard and KcwiKit in order to have it run correctly.

# KSkyWizard
This file lives in ```/kskywizard/kskywizard/utils.py```.  The patch is specifically geared towards the ```star_telluric()``` function. Here, it queries the PypeIt library to retrieve the intrinsic telluric grid and then queries the Maunakea Atmosphere grid to map it. However, if the name of the standard you used does not exist, then the script throws an error.

The patch does the following: if the name exists in the KSkyWizard, do nothing new. If the name does _not_ exist, PypeIt will assume the standard-star spectrum is a smooth polynomial:

$$
F_{\rm obs} \approx P(\lambda)T(\lambda),
$$

where

$$ 
P(\lambda) = \text{smooth polynomial approximation to the star's intrinsic continuum;}
$$
$$
T(\lambda) = \text{atmospheric transmission drawn from/fitted using the Maunakea TelFit grid.}
$$

So, it asks "What smooth underlying continuum plus what atmospheric absorption best explains this observed standard-star spectrum?", and from this polynomial is able to create the same ```*_tellmodel.fits``` and ```*_tellcorr.fits``` products that KSkyWizard normally produces.

# KcwiKit
This file lives in ```/KcwiKit/kcwikit/kcwi/kcwi.py```. The patch is geared at making sure bad pixels are correctly masked. We found that when it tries to convert NaN -> int, instead of it being 0 it gets assigned to -2147483648, which then results in enormously negative normalized flux values for those specific pixels, even after stacking.

The patch does the following: it first converts all Nan -> 0, then does the int conversion step. The 0 is important here because during the ```kcwi_stack()``` procedure, pixels assigned 0 are flagged as invalid, therefore they do not get passed to be stacked. This is the original intention of the script, and it seems that this bug with NumPy was causing the error.

# Installation

In order for these to be used, you must replace the ```utils.py``` and ```kcwi.py``` files in their respective directories. It may be a good idea to duplicate the file currently there so as to preserve the pre-patch state.
