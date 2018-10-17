#!/usr/bin/python2

# # A quick description:
# 
# This is a script to query the SDSS Galaxy Zoo database for all or a number of galaxies that they are either ellipticals or spirals (according to their classification probability) and given a specific radius.
# This version matches the GALFIT dataset where all galaxies are centered (but with real noise).

import sys

import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

from astroquery.sdss import SDSS
# see details at: https://astroquery.readthedocs.io/en/latest/api/astroquery.sdss.SDSSClass.html#astroquery.sdss.SDSSClass.query_sql

import os, random, sys, glob
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt


# CHECK THESE VARIABLES! #######################################################

path_to_SDSS_database = "/data/deep_learning/SDSS_database/centered/"

path_to_astropy_cache    = "~/.astropy/cache"
path_to_astropy_download = "~/.astropy/cache/download/py2"
path_to_astropy_download = os.path.expanduser(path_to_astropy_download)
# to expand the "~"

# EDIT initial parameters here:

# Select if you want the query to return all galaxies or only a number 
# accepted values: 'all' or number (integer)

select_galaxies = 'all'      # 'all' | 10, 100, etc

# Determine the minimum probablity 
# for their classification as a spiral or elliptical
select_class_prob = 0.8

# Select min and max radius 
# (according to their Petrosian radius at 90% in the g band)
# value is in arcsec
select_min_radius = 10
select_max_radius = 20

# Size of cropped images (around the galaxy's position), in pixels
size = 128
hsize = size/2

# select if you want to clear the directory from previous results or not
# else the script will check which files exist and will continue
clear_all = 'no'         # 'yes' | 'no'

################################################################################


pwd = os.getcwd()

os.chdir(path_to_SDSS_database)


# setting up query --------------------------------------------------------------

if select_galaxies != 'all':
    query_selected_galaxies = "TOP {}".format(select_galaxies)
else:
    query_selected_galaxies = ""

queryZooGalaxies = """SELECT {}
zs.objid, zs.ra, zs.dec, zs.nvote,
zs.p_el,zs.p_cw, zs.p_acw, 
zs.p_edge, zs.p_dk, zs.p_mg, zs.p_cs,
p.objid, p.petroR90_g, p.ra,p.dec, p.run, p.rerun, p.camcol, p.field 
FROM ZooSpec AS zs
JOIN PhotoObj AS p ON p.objid = zs.objid
WHERE zs.p_cw+zs.p_acw > {} OR zs.p_el > {}
  AND p.petroR90_g > {} AND p.petroR90_g < {} """.format(query_selected_galaxies, select_class_prob,select_class_prob,select_min_radius,select_max_radius)

# uncomment that if you want to really see the sql query
#print queryZooGalaxies
#-------------------------------------------------------------------------------

# running query ----------------------------------------------------------------

if not os.path.exists('00_to_be_downloaded.lst'):

    # > Performing query:

    resultsZoo = SDSS.query_sql(queryZooGalaxies)
    print "> Number of objects returned: ", len(resultsZoo)


    # > Storing query results:

    f = open("00_to_be_downloaded.lst", "w")
    
    f.write("# List of objectIDs left to be downloaded\n")
    f.write("#\n")

    for i in range(0,len(resultsZoo)):
        objectID = resultsZoo[i][0]
        f.write("%s\n" % objectID)
	
resultsZoo = SDSS.query_sql(queryZooGalaxies)
# reloading query         

#-------------------------------------------------------------------------------

# purging existing files -------------------------------------------------------
if clear_all == 'yes':
    print "\n! Clearing directory from all previous files !\n"
    os.system("rm *.fits *.png *.npy")
#-------------------------------------------------------------------------------

# loading list of objectIDs to be downloaded -----------------------------------

remaining_objectIDs = np.genfromtxt('00_to_be_downloaded.lst',unpack='True',comments='#')
#-------------------------------------------------------------------------------

# downloading ------------------------------------------------------------------

for i in range(0,len(resultsZoo)):
    objectID = resultsZoo[i][0]
    RA = resultsZoo[i][1]
    Dec = resultsZoo[i][2]
    runID = resultsZoo[i][-4]
    camcolID = resultsZoo[i][-2]
    fieldID = resultsZoo[i][-1]
    ellipt = float(resultsZoo[i][4])
    spiral = float(resultsZoo[i][5])+float(resultsZoo[i][6])

    print "[{}/{}] -- ID {}: ".format(i+1,len(resultsZoo), objectID),
    if ellipt>0.8:
        print "elliptical"
        name_root = "ellipt_{}".format(objectID)
    elif spiral>0.8:
        print "spiral"
        name_root = "spiral_{}".format(objectID)
    else:
        sys.exit("\n ERROR! Something is wrong with the classification")

    print name_root+".npy"
        
    # check if the file is present or not    
    if not objectID in remaining_objectIDs:
        print "\tfile exists! continuing..."
	continue
    else:
        os.system("sed -i '/" + str(objectID) + "/d' 00_to_be_downloaded.lst")
	# object will be downloaded: it is removed from hte list of remaining objectIDs


    # getting g-band image from SDSS
    img = SDSS.get_images(run=runID, camcol=camcolID, field=fieldID, band='g')

    
    """
    # in case of debugging 
    print "image info --> "
    print img[0].info()   
    print img[0][0].header
    img_data = img[0][0].data
    print img_data
    # if you want to save the original image
    img[0].writeto(name_root+str('_original.fits'))
    """

    
    naxis1 = float(img[0][0].header['NAXIS1'])
    naxis2 = float(img[0][0].header['NAXIS2'])

    # WCS <-> pixel coordinate conversion
    w = WCS(img[0][0].header)
    x, y = w.wcs_world2pix(RA, Dec, 0)
        
    if x>5+hsize and y>5+hsize and np.abs(x-naxis1)>5+hsize and np.abs(x-naxis2)>5+hsize:
    # rejecting galaxies falling near the border of the image (cannot be cropped in a square)
    # NOTE: the +5 is arbitrary and compensate for roundings and for the conversion: pixel <-> matrix index
        cropped = img[0][0].data[int(y)-hsize:int(y)+hsize,int(x)-hsize:int(x)+hsize]
        np.save(name_root, cropped)
        print "\t cropping around WCS=(%s,%s) => pix=(%s,%s)" %  (RA, Dec, x, y)

        # selecting a random sample of galaxies to save as images 
        if random.random()>(1.0-100.0/len(resultsZoo)): 
            plt.imshow(cropped, cmap='gray', vmin=np.min(cropped), vmax=0.1*np.max(cropped), origin='lower')
            plt.savefig(name_root+'.png', bbox_inches='tight')
            fits.writeto(name_root+'.fits', cropped, img[0][0].header) 
            
    else:
        print "\t object at a corner - rejecting ... "
        
    os.system("rm /tmp/tmp*")
    # removing temporary file in /tmp
    
    
    n_tmp = len(os.listdir(path_to_astropy_download))
    # numnber of tmp files generated by astropy
    
    if(n_tmp > 1000):
    
        print("|")
        print("Excessive number of tmp files (%s), quitting and removing astropy download folder" % n_tmp)
        print("|")
	
	os.system("rm -rf " + path_to_astropy_download)
    
        exit()

        # Not used:
        # os.system("ls -d -1 " + path_to_astropy_cache + "/download/py2/* | grep -v \"urlmap\" | awk '{printf \"rm %s\\n\", $1}' | /bin/sh")
        # removing temporary file in ~/.astropy/cache/download/py
        # NOTE: sort by date and remove most recent to avoid removing wrong file
    
print "\n### Total Sums ###"
print "- Galaxies:\t", len(glob.glob('*.npy'))
print "- Spirals:\t", len(glob.glob('spiral*.npy'))
print "- Ellipticals:\t", len(glob.glob('ellipt*.npy'))

# creating an image with all show cases
os.system('montage -label %f -tile 10x -geometry +0+0 *.png 00_showcase.jpg')        

os.system("rm -rf " + path_to_astropy_download)

os.chdir(pwd)
