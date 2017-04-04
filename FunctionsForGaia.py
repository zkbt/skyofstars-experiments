"This modulde contains formulas useful for the analysis of Gaia data"

import numpy as np
from math import pi
from astropy.coordinates import SkyCoord

def toXYZ(RA,DEC,Distance): #Parallax
    """
    This function converts RA, Dec and Distance into
    X,Y and Z coordinates.

    Parameters:
    --------------------------------------------------
    RA = Array of the right ascencion of each star.
    DEC = Array of the declination of each star.
    Distance = Array of the distance from Earth of
    each star, in AU.

    Returns:
    ---------------------------------------------------
    X,Y,Z arrays of the x,y,z coordinates of all the
    stars included, in radians.

    """
    
    c = SkyCoord(RA,DEC,frame='icrs',unit='rad') #use skycoord modules so coordinates are in radians
    longitude = c.galactic.l.rad
    latitude = c.galactic.b.rad
    
    r=Distance*(np.cos(latitude))
    z=Distance*(np.sin(latitude))
    y=r*(np.sin(longitude))
    x=r*(np.cos(longitude))

    return x,y,z

def MagnitudeToFlux(M):
    """
    This function uses the visual magnitude of a group of star
    to find out their energy flux (power/area).

    Parameters:
    -------------------
    M = numpy array of the magnitudes of a group of stars

    Returns:
    -------------------
    Fluxes = numpy array of the fluxes each of the stars emits in units Flx
    (Need flux and magnitude of another star to convert in SI units Watts/m^2 )

    """
    # Use magnitudes formula, need to use the flux of a Star with magnitude 0
    # this will make the formula return fluxes on that unit
    Fluxes = (10**(-M/2.5))
    return Fluxes



# distance measured in parcsecs. Star and its parallax angle relation  d=1/theta
#2x10**5 asrcseconds = 1 radian

def ParallaxToDistance(P):
    """
    This function converts the parallax of a group
    of stars into their distance from Earth.

    Parameters:
    ----------------------------------------------
    P = Array of the parallaxes of a group of stars in arcseconds

    Returns:
    ----------------------------------------------
    Array of the distances from the Earth of a group of stars in AU
    """
    #Calculate distance in parseconds:
    distance = 1/P

    #convert to AU:
    distanceAU= (distance)*(360/(2*pi)*3600)

    return distanceAU

def AbsoluteMagnitude(distances,magnitudes):
    """
    Use aparent magnitudes and distances to calculate the
    absolute magnitude of the stars.

    Parameters:
    ----------------------------------------
    Distances =  Array of the star's distances, calculated from the parallax
    of the Tgas files.
    Magnitudes = Array of magnitudes, taken from the Tgas files. 

    Returns:
    -----------------------------------------
    AbM = Array of the absolute magnitudes of the stars 
    """
    inParsec= distances*(1.0/(360.0/(2.0*pi)*3600.0))
    AbM = magnitudes-5*np.log10(inParsec/10)
    
    return AbM
