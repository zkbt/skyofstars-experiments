'This module contains functions to help load and use a Gaia fits file'

import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline - I'm not sure if a module requires this. Ask ZKBT.
from math import pi
import astropy.io.ascii as a
import astropy.io.fits
import astropy.table

#use module created specially for this:
import FunctionsForGaia as F

def extractStars(filename):

    """
    This function will read in a TGAS FITS file and extract the information that the
    Sky Of Stars project will be working with.
    This function will replace parallaxes that are negative or with an uncertainty greater than 50%
    with the value 10^-6 arcseconds. This is done to use all the stars without having them be too             different than the rest

    Parameters:
    ------------------------------------------------------------------
    Filename = TGAS.fits file

    Returns:
    ------------------------------------------------------------------
    Arrays for X,Y and Z positions and arrays for right ascencion, declination, magnitudes,
    distances and absolute magnitudes.
    label: A specific number that comes from the file name, used to identify each file. 

    """
    #Read the file into a table:
    hdu = astropy.io.fits.open(filename)
    table = astropy.table.Table(hdu[1].data)

    #Create a statement to get rid of stars with parallax error larger than 20% or with negative distances:
#    fraction = np.abs(table['parallax_error']/table['parallax'])
#    ok1 = fraction < 0.2
#    ok2 = table['parallax']>0
#    ok = ok1*ok2
# add [ok] after .data to make use of the info above    
    
    
    
    #From the table, extract parallax and other useful info, take only those rows for which [ok] is true:
    
    Parallax= table['parallax'].data #.data takes only the numbers, getting rid of the tittle of the         column
    
    Parallax = Parallax/1000 # ZKBT: the GAIA parallaxes are in units of milliarcseconds; let's convert       to arcsec
    negative = Parallax<0
    Parallax[negative] = 10**(-6) #replace al negative parallaxes with 10^-6 arcseconds
    uncertainty = (table['parallax_error']/table['parallax']) > 0.5
    Parallax[uncertainty] = 10**(-6) #replace all parallaxes with uncertainty greater than 50% with 10^-6     arcseconds
    Dec = (table['dec'].data)*(pi/180) #change degrees to radians.
    Fluxes= table['phot_g_mean_flux'].data
    Magnitudes = table['phot_g_mean_mag'].data
    
    RA = (table['ra'].data)
    over180 = RA > 180 #want RA values to be from -180 to 180
    RA[over180] = (RA[over180] - 360)
    RA = RA*(pi/180) #change degrees to radians.
    
    #Use formulas to produce other useful arrays:
    Distances = F.ParallaxToDistance(Parallax)
    X,Y,Z= F.toXYZ(RA,Dec,Distances)
    AbsoluteMagnitudes = F.AbsoluteMagnitude(Distances,Magnitudes)

    #Create a label dependant on the file name:
    temp = filename.replace('.fits','.png')
    tempList = temp.split('_')
    label = tempList[-1]

    return X,Y,Z,RA,Dec,Magnitudes,Distances,AbsoluteMagnitudes,Fluxes,label

def PlotStars(X,Y,Z,RA,Dec,Distances,Fluxes,label):
    """
    This function creates demonstraation plots and saves each in a file with name
    dependent on the Gaia Fits file the information it comes from.

    Parameters:
    -------------------------------------------------------------------
    X,Y,Z,RA,DEC,Distances arrays calculated or extracted using the function "extractStars".
    label: string dependant on file name, also obtained with the "extractStars"function.

    Returns and Saves:
    --------------------------------------------------------------------
    4 plots: x vs y position, x vs z position, y vs z position and RA vs Dec.
    The title of each plot is dependant on the label, which is dependant on the file the information came 
    from

    """
    plt.figure(figsize=(9,9))
    plt.scatter(X,Y,s = Fluxes/1000000,alpha=0.1)#make size flux dependant and each point transparent
    plt.xlim(np.min(X),np.max(X))
    plt.ylim(np.min(Y),np.max(Y))
    plt.xlabel('X Coordinates')
    plt.ylabel('y coordinates')
    plt.title('X vs Y for {}'.format(label))
    plt.savefig('XY_{}'.format(label)) #save the graph

    plt.figure(figsize=(9,9))
    plt.scatter(X,Z,s = Fluxes/1000000,alpha=0.1)
    plt.xlim(np.min(X),np.max(X))
    plt.ylim(np.min(Z),np.max(Z))
    plt.xlabel('X Coordinates')
    plt.ylabel('Z coordinates')
    plt.title('X vs Z for {}'.format(label)) #Make the graph title and filename depend on label.
    plt.savefig('XZ_{}'.format(label))

    plt.figure(figsize=(9,9))
    plt.scatter(Y,Z,s = Fluxes/1000000,alpha=0.1)
    plt.xlim(np.min(Y),np.max(Y))
    plt.ylim(np.min(Z),np.max(Z))
    plt.xlabel('Y Coordinates')
    plt.ylabel('Z coordinates')
    plt.title('Y vs Z for {}'.format(label))
    plt.savefig('YZ_{}'.format(label))

    plt.figure(figsize=(9,9))
    plt.subplot(111, projection="aitoff")
    plt.scatter(RA,Dec,color='aqua')
    plt.xlabel('Right Ascension')
    plt.ylabel('Declination')
    plt.grid(True)
    plt.title('RA and Dec for {}'.format(label))
    plt.savefig('RADec_{}'.format(label))

def CreateTxtFile(X,Y,Z,AbsoluteMagnitudes,label):
    """
    This functions create a label dependant, it contains a list of
    X,Y,and Z positions, and absolute magnitude of each star.

    Parameters:
    --------------------------------------------------------
    X,Y,Z,AbsoluteMagnitudes arrays.
    label:string dependant on file name, obtained with the "extractStars"function.

    Saves:
    ---------------------------------------------------------
    "label".txt file with X,Y,Z and absolute magnitudes list for each star.
    """

    t = astropy.table.Table([X, Y,Z,AbsoluteMagnitudes], names=['X', 'Y','Z','Absolute Magnitudes'])
    a.write(t, '{}.txt'.format(label),format='fixed_width')


def convertTgasFile(filename):

    """
    This function uses the functions extractStars, CreateTcxtFile, PlotStars and fileForFiske.

    """
    X,Y,Z,RA,Dec,Magnitudes,Distances,AbsoluteMagnitudes,Fluxes, label = extractStars(filename)
    CreateTxtFile(X,Y,Z,AbsoluteMagnitudes,label)
    PlotStars(X,Y,Z,RA,Dec,Distances,Fluxes,label)
    fileForFiske(filename)
    
    return "Good job Luci"

    
def fileForFiske(filename):
    """
    Create a text file for Fiske.
    Parameters:
    -------------------------------
    filename, of a Gaia Tgas file. 
    
    Returns:
    -------------------------------
    A file with a name determined by the filename.
    Text file contains a column for each: ID number, X,Y and Z positions, dx, dy, dz velocities, BV and       absolute magnitude.
    dx, dy, dz and bv are empty columns but work as place holders.
    """
    X,Y,Z,RA,Dec,Magnitudes,Distances,AbsoluteMagnitudes,Fluxes,label = extractStars(filename)
    ID = np.arange(len(X))
    DX = np.zeros_like(X)
    DY = np.zeros_like(X)
    DZ = np.zeros_like(X)
    BV = np.zeros_like(X)

    # ZKBT: slight change, to save in parsecs instead of AU
    au_per_pc = 180/np.pi*60*60
    X = X/au_per_pc
    Y = Y/au_per_pc
    Z = Z/au_per_pc
    table = astropy.table.Table([ID, X, Y, Z, DX, DY, DZ, BV, AbsoluteMagnitudes],
        names=['ID', 'X', 'Y', 'Z', 'DX', 'DY', 'DZ', 'BV', 'AbsMag'])

    # write out the table
    astropy.io.ascii.write(table, 'GaiaFile{}forFiske.dat'.format(label.split('.')[0]),format='fixed_width',delimiter='\t')
