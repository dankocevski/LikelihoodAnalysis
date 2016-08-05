#!/usr/bin/env python
#Testing
print "\nLikelihood Analysis Tool v1.0"
print "Support Contact: Daniel Kocevski (daniel.kocevski@nasa.gov)\n"
print "Importing modules..."
import os
import subprocess
import pyfits
import fileinput
import time
import random
import sys
import glob
import numpy
import subprocess
from functools import wraps
import errno
import os
import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plot
from matplotlib import dates as mdates
import matplotlib.colors as mcolors
import getpass
import make3FGLxml
import traceback
import operator
from math import *
import numpy
import commands
import traceback
import pyfits
import time
import numpy
import aplpy
from astropy.io import fits
from astropy.wcs import WCS
import pickle
import pywcs
from GtApp import GtApp
import pyLikelihood
from UnbinnedAnalysis import *
import BinnedAnalysis
from UpperLimits import UpperLimits
import IntegralUpperLimit as IUL
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

print "Done."


##########################################################################################
def jd_from_MET( met ):

	# Converts from MET to Julian Date
    julianDate = skymaps.JulianDate((skymaps.JulianDate_missionStart().seconds() + met)/skymaps.JulianDate.secondsPerDay)

    return julianDate

##########################################################################################    

def getSunPosition( met ):

    # Environment variable (you need FTOOLS installed and configured)
    os.environ['TIMING_DIR']=os.path.join(os.environ['HEADAS'],"refdata")

    # Get the sun direction
    sun = skymaps.SolarSystem(skymaps.SolarSystem.SUN)
    SunSkyDir = sun.direction(jd_from_MET(met))

    return SunSkyDir


##########################################################################################

def computeEnergyFlux(photonflux, index, emin=100, emax=10000):

    #print 'index= %s '
    if index == -2:
        index = -2.000001
    elif index == -1:
        index = -1.000001
        pass

    energyflux = photonflux*(1.+index)/(2+index)*(pow(emax,index+2)-pow(emin,index+2))/(pow(emax,index+1)-pow(emin,index+1))
    print 'index= %s , <E>= %s' %(index,energyflux/photonflux)

    return energyflux#*MeV2erg

##########################################################################################

def BOOL(string):

	if 'True' in string:
		return True
	else:
		return False

##########################################################################################

class TimeoutError(Exception):
	pass

class timeout:
	def __init__(self, seconds=1, error_message='Timeout'):
		self.seconds = seconds
		self.error_message = error_message
	def handle_timeout(self, signum, frame):
		raise TimeoutError(self.error_message)
	def __enter__(self):
		signal.signal(signal.SIGALRM, self.handle_timeout)
		signal.alarm(self.seconds)
	def __exit__(self, type, value, traceback):
		signal.alarm(0)



##########################################################################################

def AngularSeperation(ra1,dec1,ra2,dec2):

	d1 = radians(90.-dec2)
	d2 = radians(90.-dec1)
	d3 = radians(ra2-ra1)

	sep = cos(d1) * cos(d2)

	sep += sin(d1) * sin(d2)*cos(d3)

	sep = degrees( acos(sep))
	return sep # in degrees


##########################################################################################

def GetCatalogFluxValues(like):

	print 'Extracting catalog flux values...'
	Sources = like.sourceNames()
	catalogFluxValues = {}

	for source in Sources:

		if (('gll_iem' in source) or ('iso_source' in source)):
			print 'Skipping %s' % source
		else:
			fluxValue = like.flux(source)
			catalogFluxValues[source] = fluxValue

			print ' %s %s' % (source, fluxValue)

	return catalogFluxValues


##########################################################################################

def FreezeAllSources(like):

	print 'Freezing all sources...'

	Sources = like.sourceNames()	
	for source in Sources:
		print ' %s' % source

		parameters = like.freePars(source)

		for parameter in parameters:
			parameterIndex = like.par_index(source,parameter.getName())

			# Freeze the parameter
			like.freeze(parameterIndex)

	print '\nDone.'


##########################################################################################

def FreeSourceParameters(like, srcname, TS=1, TSlim=0):

	"""Free source parameters depending on cataloged TS """

	#print "\n  >> Working on source %s with 2FGL_TS=%.1f, freeing its norm"%(srcname,TS)
	
	#src.isFree='yes' # change the flag
	#TSlim = 25 # this is the limit TS to free also the index


	# --- now THAW the normalization, free the index only if the TS>TSlim
	spec_type = GetSpectralType(like,srcname)

	if(spec_type=="PowerLaw"):
		norm_idx = like.par_index(srcname,'Prefactor')
		spec_idx = like.par_index(srcname,'Index')

		like.thaw(norm_idx)

		par_val ={'Prefactor':like[norm_idx].value()}


		if(TS>TSlim):
			like.thaw(spec_idx)
			par_val ={'Prefactor':like[norm_idx].value(),
					  'Index':like[spec_idx].value()
					  }
			 

	elif(spec_type=="LogParabola"):
		norm_idx = like.par_index(srcname,'norm')
		spec_idx = like.par_index(srcname,'alpha')
		beta_idx = like.par_index(srcname,'beta')
	
		like.thaw(norm_idx)

		par_val ={'norm':like[norm_idx].value()}


		if(TS>TSlim):
			like.thaw(spec_idx)
			like.thaw(beta_idx)

			par_val ={'norm' : like[norm_idx].value(),
					  'alpha': like[spec_idx].value(),
					  'beta' : like[beta_idx].value()
					  }


	elif(spec_type=="PowerLaw2"):
		norm_idx = like.par_index(srcname,'Integral')
		spec_idx = like.par_index(srcname,'Index')

		like.thaw(norm_idx)
		par_val ={'Integral':like[norm_idx].value()}

		if(TS>TSlim):
			like.thaw(spec_idx)
			par_val ={'Integral': like[norm_idx].value(),
					  'Index'   : like[spec_idx].value()
					  }
			
			
	elif(spec_type=="ExpCutoff"):
		norm_idx = like.par_index[srcname,'Prefactor']
		spec_idx = like.par_index[srcname,'Index']

		like.thaw(norm_idx)
		par_val ={'Prefactor' : like[norm_idx].value()}


		if(TS>TSlim):
			like.thaw(spec_idx)
			par_val ={'Prefactor' : like[norm_idx].value(),
					  'Index'     : like[spec_idx].value()
					  }


	elif(spec_type=="PLSuperExpCutoff"):
		norm_idx = like.par_index(srcname,'Prefactor')
		spec_idx = like.par_index(srcname,'Index1')

		like.thaw(norm_idx)
	 
		par_val ={'Prefactor' : like[norm_idx].value()}
 

		if(TS>TSlim):
			like.thaw(spec_idx)

			par_val ={'Prefactor' : like[norm_idx].value(),
					  'Index1'    : like[spec_idx].value()
					  }

		
	else:
		
		print ">>FreeNextBrightSource::ERROR spectral model %s not defined"%(spec_type)

	# === return the parameter value before the fit
	return par_val


##########################################################################################

def GetSpectralType(like,srcname):

	string = str(like[srcname])

	tok1 = string.find("Spectrum")
	tok2 = tok1+string[tok1:].find("\n")

	spec_type = string[tok1+10:tok2]

	if(spec_type!="PowerLaw2" and spec_type!="PowerLaw" and
		spec_type!="LogParabola" and  spec_type!="PLSuperExpCutoff"):

		print " ERROR >>GetSpectralType:: spectrum `%s' not known for src %s" % (spec_type,srcname)
		sys.exit()

	return spec_type

##########################################################################################


def astro_query(options, store='store', verbose=True):
	command = "/afs/slac/u/gl/glast/astroserver/prod/astro"
	for item in options:
		if options[item] is not None:
			command += ' --%s %s' % (item, options[item])
	command += ' ' + store
	if verbose:
		print command
	subprocess.call(command, shell=True)
	for item in glob.glob('*_fix_checksums.sh'):
		os.remove(item)

##########################################################################################

def getFT1(ra, dec, radius, tmin, tmax, emin=100., emax= 3e5,
		   evclassname="Source", outfile='FT1.fits', verbose=True):
#    options = {'event-sample' : 'P7.6_P120_ALL',
#    options = {'event-sample' : 'P7.6_P130_ALL',
	options = {'event-sample' : 'P7_P202_ALL',
			   'output-ft1' : outfile,
			   'minEnergy' : emin,
			   'maxEnergy' : emax,
			   'minTimestamp' : tmin,
			   'maxTimestamp' : tmax,
			   'ra' : ra,
			   'dec' : dec,
			   'radius' : radius}
	astro_query(options, store='store', verbose=verbose)


##########################################################################################

def getFT2(tmin, tmax, outfile='FT2.fits', verbose=True):
	options = {'event-sample' : 'P7_P202_ALL',
			   'output-ft2-30s' : outfile,
			   'minTimestamp' : tmin,
			   'maxTimestamp' : tmax}
	astro_query(options, store='storeft2', verbose=verbose)


##########################################################################################

def make_ds9_reg(ra, dec, outfile='CandidateSource.reg'):
	tpl = """# Region file format: DS9 version 4.0
# Filename: cmap_sum.fits
global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source
fk5
point(%.3f,%.3f) # point=circle
"""
	output = open(outfile, 'w')
	output.write(tpl % (ra, dec))
	output.close()


##########################################################################################

def writeXML(self, ra, dec, galpropModel, isotropicModel, outfile='SourceModel.xml'):
	output = open(outfile, 'w')
	output.write(Source_xml % (ra, dec, galpropModel, isotropicModel))
	output.close()


##########################################################################################

def AddCandidateSource(ra, dec, xmlModel, fixIndex=True):

	print "\nAdding a candidate source at RA=%s Dec=%s to the xml model..." % (ra, dec) 

	if fixIndex == True:
		print "Spectral index set fixed."
		freeFlag = "0"
	else:
		print "Spectral index set free."		
		freeFlag = "1"

	# GRBAnalysis
# 	CandidateSource = """
# <source name="CandidateSource" type="PointSource">
# 	<spectrum type="PowerLaw2">
# 	  <parameter free="1" max="1000000" min="1e-10" name="Integral" scale="1e-06" value="0.01" />
# 	  <parameter free="1" max="0.01" min="-6" name="Index" scale="1" value="-2.2" />
# 	  <parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="100.0"/>
# 	  <parameter free="0" max="300000.0" min="20.0" name="UpperLimit" scale="1.0" value="1e5"/>
# 	</spectrum>
# 	<spatialModel type="SkyDirFunction">
# 	  <parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="%.3f"/>
# 	  <parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="%.3f"/>
# 	</spatialModel>
#   </source>
# </source_library>
# """ % (float(ra), float(dec))

# 	# pyapp.py
# 	CandidateSource = """
# <source name="CandidateSource" type="PointSource">
# 	<spectrum type="PowerLaw2">
# 	  <parameter free="1" max="1000" min="1e-06" name="Integral" scale="1e-08" value="1" />
# 	  <parameter free="1" max="0.01" min="-6" name="Index" scale="1" value="-2.2" />
# 	  <parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="100.0"/>
# 	  <parameter free="0" max="300000.0" min="20.0" name="UpperLimit" scale="1.0" value="1e5"/>
# 	</spectrum>
# 	<spatialModel type="SkyDirFunction">
# 	  <parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="%.3f"/>
# 	  <parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="%.3f"/>
# 	</spatialModel>
#   </source>
# </source_library>
# """ % (float(ra), float(dec))

	# Custom with frozen spectral index
	CandidateSource = """
<source name="CandidateSource" type="PointSource">
	<spectrum type="PowerLaw2">
	  <parameter free="1" max="1e10" min="1e-20" name="Integral" scale="1e-08" value="1" />
	  <parameter free="%s" max="0.01" min="-6" name="Index" scale="1" value="-2.1" />
	  <parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="100.0"/>
	  <parameter free="0" max="300000.0" min="20.0" name="UpperLimit" scale="1.0" value="1e5"/>
	</spectrum>
	<spatialModel type="SkyDirFunction">
	  <parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="%.3f"/>
	  <parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="%.3f"/>
	</spatialModel>
  </source>
</source_library>
""" % (freeFlag, float(ra), float(dec))


	xmlModelContents = open(xmlModel).read()
	xmlModelContents = xmlModelContents.replace('</source_library>',CandidateSource)
	xmlModelContentsAppended = open(xmlModel,'w')
	xmlModelContentsAppended.write(xmlModelContents)
	xmlModelContentsAppended.close()
	print 'Done.'



##########################################################################################

def RemoveCandidateSource(xmlModel, xmlModelModified, candidateSource='CandidateSource', RemoveSource=True, FixSources=False):

	# Open the files
	infile = open(xmlModel,'r')
	outfile = open(xmlModelModified,'w')

	# Read in the lines from the original flare file
	lines = infile.readlines()

	# Loop through each line 
	doFix = True
	eraseLine = False

	for line in lines:

		newline = line

		# Fix all non-diffuse model components
		if 'type="DiffuseSource">' in newline:
			print 'Keeping Component:'
			print newline.rstrip()
			doFix = False

		# Make sure the CandidateSource is not affected
		if 'name="%s"' % candidateSource in newline:
			print 'Keeping Component:'
			print newline.rstrip()
			doFix = False

		if '<source name=' in line:
			SourceName = line

		if 'free="1"' in line and doFix == True and FixSources == True:
			print 'Fixing Component:'
			print SourceName.rstrip()
			newline = line.replace('free="1"','free="0"')

		if '</source>' in line and doFix == False:
			doFix = True

		# Remove the candidate source from the model    
		if candidateSource in line and RemoveSource == True:
			print 'Removing Component:'
			print candidateSource
			newline = ''
			eraseLine = True

		if '</source>' in line and eraseLine == True:
			newline = ''
			eraseLine = False

		if eraseLine == True:
			newline = ''

		# Save the modified line
		outfile.write(newline)

	print "\n\nWriting modified xml file to: %s" % xmlModelModified
	infile.close()
	outfile.close()
	print "Done."


##########################################################################################

def ModifySourceModel(xmlModel, xmlModelModified, candidateSource='CandidateSource', RemoveSource=True, FixSources=True):

	# Open the files
	infile = open(xmlModel,'r')
	outfile = open(xmlModelModified,'w')

	# Read in the lines from the original flare file
	lines = infile.readlines()

	# Loop through each line 
	doFix = True
	eraseLine = False

	for line in lines:

		newline = line

		# Fix all non-diffuse model components
		if 'type="DiffuseSource">' in newline:
			print 'Keeping Component:'
			print newline.rstrip()
			doFix = False

		# Make sure the CandidateSource is not affected
		if 'name="%s"' % candidateSource in newline:
			print 'Keeping Component:'
			print newline.rstrip()
			doFix = False

		if '<source name=' in line:
			SourceName = line

		if 'free="1"' in line and doFix == True and FixSources == True:
			print 'Fixing Component:'
			print SourceName.rstrip()
			newline = line.replace('free="1"','free="0"')

		if '</source>' in line and doFix == False:
			doFix = True

		# Remove the candidate source from the model    
		if candidateSource in line and RemoveSource == True:
			print 'Removing Component:'
			print candidateSource
			newline = ''
			eraseLine = True

		if '</source>' in line and eraseLine == True:
			newline = ''
			eraseLine = False

		if eraseLine == True:
			newline = ''

		# Save the modified line
		outfile.write(newline)

	print "Writing modified xml file to:\n%s" % xmlModelModified
	infile.close()
	outfile.close()
	print "\nDone.\n"

##########################################################################################

def ExtractSources(xmlModel, sourceNames):

	# Open the files
	infile = open(xmlModel,'r')
	outfile = open(sourceNames,'w')

	# Read in the lines from the original flare file
	lines = infile.readlines()

	# Go through each line and extract the source names
	SourceNames =[]
	Values = []
	for line in lines:
		if '<source name=' in line:

			# Extract the source name
			SourceName = line
			SourceName = SourceName.replace('  <source name="','')
			SourceName = SourceName[0:SourceName.find('"')]
			SourceName = SourceName.replace('_2FGL', '2FGL_')

			# Extract the url value to the passed
			Value = SourceName
			Value = SourceName.replace('+','%2B')

			# Add the source names and values to their respective arrays
			SourceNames.append(SourceName)
			Values.append(Value)



	# Remove the candidate source, and the galactic, and isotropic components
	SourceNames = SourceNames[1:-2]

	# Add a none source
	SourceNames.insert(0,'none')

	# Write the source names to a file
	for Value, Source in zip(Values, SourceNames):
		line = '<option value="' + Value + '">' + Source + '</option>\n'
		outfile.write(line)

	print "Writing source names file to:\n%s" % sourceNames
	infile.close()
	outfile.close()
	print "Done."



##########################################################################################

def SetPfilesDirectory(pfile_dir):
	"""each thread/job which uses FTOOLS must have its own
	PFILES dir"""
	
	try:
		hea_pdir = os.getenv("HEADAS")+"/syspfiles/"
		gt_pdir  = os.getenv("INST_DIR")+"/syspfiles/"
	except:
		print '\n*** Science tools has not be inialized! ***\n'
		print 'Exiting!'
		sys.exit()
	

	if(os.path.isdir(pfile_dir)==False):
		print "\nCreating custom pfile directory:\n%s" % pfile_dir        
		cmd = "mkdir -p "+pfile_dir
		os.system(cmd)

	# --- now copy all the .par files
	cmd = "cp %s/*.par %s"%(hea_pdir,pfile_dir)
	os.system(cmd)
	cmd = "cp %s/*.par %s"%(gt_pdir,pfile_dir)
	os.system(cmd)
	
	# ---Set the new environmental PFILES    
	#os.putenv("PFILES",pfile_dir)
	os.environ['PFILES'] = pfile_dir

	# --- test
#   print "Testing: temporary PFILES location "
#   print os.listdir(pfile_dir)

	return pfile_dir



##########################################################################################	

def customRALabel(deg):
	if (deg == 360):
		return ''
	return r'%s$^{\circ}$' % deg

##########################################################################################

def customDecLabel(deg):
	return r'%s$^{\circ}$' % deg

##########################################################################################	

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

##########################################################################################

# def MakeButterflyPlot(sourceName):
# 	E1 = like2.model[sourceName].funcs['Spectrum'].getParam('LowerLimit').value()
#     E2 = like2.model[sourceName].funcs['Spectrum'].getParam('UpperLimit').value()
# 	gamma = like2.model[sourceName].funcs['Spectrum'].getParam('Index').value()
# 	I = like2.model[sourceName].funcs['Spectrum'].getParam('Integral').value()
# 	cov_gg = like2.covariance[16][16]
# 	cov_II = like2.covariance[15][15]
# 	cov_Ig = like2.covariance[15][16]
# 	print "Index: " + str(gamma) + " +/- " + str(math.sqrt(cov_gg))
# 	print "Integral: " + str(I) + " +/- " + str(math.sqrt(cov_II))
# 	print E1,E2
# 	epsilon = (E2/E1)**(1-gamma)
# 	logE0 = (math.log(E1) - epsilon*math.log(E2))/(1-epsilon) + 1/(gamma-1) + cov_Ig/(I*cov_gg)
# 	E0 = math.exp(logE0)
# 	print E0


	
##########################################################################################

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

##########################################################################################

def RemoveSources(DuplicateSources, Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique):

	for Source in DuplicateSources:
		i = numpy.where(Names_SourcesUnique == Source)
		Names_SourcesUnique = numpy.delete(Names_SourcesUnique,i)
		RA_SourcesUnique = numpy.delete(RA_SourcesUnique,i)
		DEC_SourcesUnique = numpy.delete(DEC_SourcesUnique,i)
		pass

	return Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique

##########################################################################################

def ExtractCoordinates(xmlModel, Source):

	# Define some default values
	sourceFound = False
	sourceRA = None
	sourceDEC = None

	# Loop through the xml file and look for the requested source 
	for line in fileinput.input([xmlModel]):
		if Source in line:
			sourceFound = True

		# Extract the RA and Dec for extended sources
		if ('spatialModel file=' in line) and (sourceFound == True):
			lowerIndex = line.find('"')+1
			upperIndex = lowerIndex + line[lowerIndex:].find('"')
			SpatialModel = line[lowerIndex:upperIndex]

			# Open the spatial model and extract the coordinates
			fitsfile = pyfits.open(SpatialModel)
			header = fitsfile[0].header
			sourceRA = float(header['CRVAL1'])
			sourceDEC = float(header['CRVAL2'])

		# Extract the RA for a point source
		if ('name="RA"' in line) and (sourceFound == True):
			sourceRA = line[line.find('value='):]
			sourceRA = sourceRA.replace('/>','')
			sourceRA = sourceRA.replace('value=','')
			sourceRA = sourceRA.replace('"','')
			sourceRA = float(sourceRA)

		# Extract the Dec for a point source		
		if ('name="DEC"' in line) and (sourceFound == True):
			sourceDEC = line[line.find('value='):]
			sourceDEC = sourceDEC.replace('/>','')			
			sourceDEC = sourceDEC.replace('value=','')
			sourceDEC = sourceDEC.replace('"','')
			sourceDEC = float(sourceDEC)

		# Stop recording values
		if ('</source>' in line):
			sourceFound = False

	fileinput.close()

	return sourceRA, sourceDEC


##########################################################################################

def plotImage(image, wcs_astropy, filename=None, region=None, colorbarLabel=None, ROI=None, stretchColorValue=None, maxValue=100):

	# Convert the astropy wcs objec to a pywcs wcs object
	header_astropy = wcs_astropy.to_header()

	header_pyfits = pyfits.Header()

	for key in header_astropy.keys():
		header_pyfits[key] = header_astropy[key]

	# Create a hdu object
	hdu = pyfits.PrimaryHDU(header=header_pyfits)

	# Create a new wcs object from the hdu
	wcs = pywcs.WCS(hdu.header)

	# Add the dimensions to the wcs object
	xsize = int(image.shape[0])-1
	ysize = int(image.shape[1])-1
	wcs.naxis1 = xsize
	wcs.naxis2 = ysize

	# Create a new plot figure
	fig = plot.figure(figsize=(8.5, 8.28) )

	# Create an image using the wcs object
	fig_wcs = aplpy.FITSFigure(wcs, figure=fig)

	if stretchColorValue != None:
		image[0][0] = stretchColorValue

	if maxValue != None:
		i = numpy.where(image >= maxValue)
		image[i] = maxValue

	# Plot the data
	plot.imshow(image, interpolation='nearest')

	# Add a colorbar
	cbar = plot.colorbar(pad=0.01, aspect=25, shrink=0.84)

	# Add a label to the colorbar
	if colorbarLabel != None:
		cbar.set_label(colorbarLabel)

	# Add grid lines and labels
	fig_wcs.add_grid()
	fig_wcs.axis_labels.set_font(size='large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')

	# Tick labels:
	fig_wcs.tick_labels.show()
	fig_wcs.tick_labels.show_x()
	fig_wcs.tick_labels.show_y()

	try:
		fig_wcs.tick_labels.set_xformat('ddd:dd')
		fig_wcs.tick_labels.set_yformat('ddd:dd')
	except:
		fig_wcs.tick_labels.set_xformat('%0.2f')
		fig_wcs.tick_labels.set_yformat('%0.2f')


	# Overlay a region file, if available
	if region != None and os.path.isfile(region) == True:

		# Calculate the corner coordinates
		ra_max = wcs_astropy.all_pix2world(0,0,0)[0]
		dec_min = wcs_astropy.all_pix2world(0,0,0)[1]
		ra_min = wcs_astropy.all_pix2world(xsize-1,ysize-1,0)[0]
		dec_max = wcs_astropy.all_pix2world(xsize-1,ysize-1,0)[1]

		print "\nSuperimposing ds9 region..."
		with open(region, 'r') as file:
			for line in file:
				line=line.strip()
				if 'point' not in line:
					continue
				if ";" in line:
						coord=line.split("#")[0].split(";")[0]
						ra, dec=line.split("#")[0].split(";")[1].replace("point(", "").replace(")", "").split(",")
						name=line.split("{")[-1].replace(" ", "").replace("}", "")
						if coord!="J2000":
							print "pyapp.fitstopng: region support is experimental, only J2000 ra-dec is implemented."
							fig.save(outfile)
							return
				else:
					ra, dec=line.split("#")[0].replace("point(", "").replace(")", "").split(",")
					name=line.split("{")[-1].replace(" ", "").replace("}", "")
					xmin, xmax, ymin, ymax=0., 360., -90., 90.      #temporary hack for different format of region files

				# Make sure the coordinates are floats
				ra = float(ra)
				dec = float(dec)

				# Only plot sources that are within the image area
				if ra > ra_min and ra < ra_max and dec > dec_min and dec < dec_max:
					print "adding source", name, " at (ra=%f, dec=%f) deg" % (ra, dec) 
					fig_wcs.show_markers(ra, dec, edgecolor='lightgray', facecolor='lightgray', marker='x', alpha=0.6, s=40, linewidths=1.5, clip_on=True)
					fig_wcs.add_label(ra, dec-0.35, name, color='lightgray', clip_on=True, alpha=0.6)

	if ROI != None:

		ra_center = wcs_astropy.to_header()['CRVAL1']
		dec_center = wcs_astropy.to_header()['CRVAL2']
		# x_center = wcs_astropy.wcs_world2pix(ra_center,dec_center,0)[0].item()
		# y_center = wcs_astropy.wcs_world2pix(ra_center,dec_center,0)[1].item()

		if type(ROI) is list:
			for roi in ROI:
				fig_wcs.show_circles(ra_center, dec_center, roi, linestyle='dashed', edgecolor='lightgray', facecolor='none')
		else:
			fig_wcs.show_circles(ra_center, dec_center, ROI, linestyle='dashed', edgecolor='lightgray', facecolor='none')


	# Transform the image variable name to a filename
	if filename == None:
		for key, value in list(locals().iteritems()):
			if value is image:
				filename = key + '.png'

	# Save the figure
	print '\nSaving image to: %s\n' % filename
	fig_wcs.save(filename)


##########################################################################################

def SourceAnalysis(sourceName, ra, dec, tmin, tmax, emin=100, emax=1e5, tsMin=25, irfs='P8R2_SOURCE_V6', ROI=12, zmax=105, association='none', dra=7, ddec=7, tsmapBinSize=0.15, getData=True, generateFiles=True, performLikelihoodFit=True, makeSourceMap=False, maketsmap=False, makeSummaryMap=False, makeModelMap=False, cleanup=True, cleanupAll=False, nuke=False, justGetData=False, fixedModel=True, statistic='UNBINNED', optimizer='MINUIT', skipDiffuseResponse=False, resultsPrefix='Results', performRefinedFit=False, removeWeakSources=False, plotFit=True, makeModel=True, fixIndex=True):

	# Import additional libraries
	import traceback
	import pyfits
	import time
	import numpy

	# Make sure the required arguments are floats
	ra = float(ra)
	dec = float(dec)
	tmin = int(float(tmin))
	tmax = int(float(tmax))
	emin = float(emin)
	emax = float(emax)
	tsMin = float(tsMin)

	print "\nPerforming analysis on: %s" % sourceName
	print "tmin = %s, tmax = %s" % (tmin, tmax)
	print "RA = %s, Dec = %s" % (ra, dec)

	if maketsmap == True:
		print '\nTS map parameters:'
		print "dRA = %s, dDec = %s, binsize=%s" % (dra, ddec, tsmapBinSize)


	# Define the home directory
	LikelihoodDirectory = '/nfs/slac/g/ki/ki08/kocevski/Likelihood'

	# Use the working directory if the likelihood directory isn't available.
	if(os.path.isdir(LikelihoodDirectory) == False):
		LikelihoodDirectory = os.getcwd()


	# Define an output directory
	OutputDirectory = "%s/%s/%s" % (LikelihoodDirectory, resultsPrefix, sourceName)

	# Check if the user has write access to this location.  If not, create a subdirectory in the currnet working directory
	if os.access(OutputDirectory, os.W_OK) == False:
		OutputDirectory = os.getcwd() + '/%s' % sourceName


	# Get the job id, if it exists
	JobID = os.environ.get('LSB_JOBID')

	# Wait a period of time before starting in order to not crash the asf/nsf disks at SLAC, but only if this is a batch job
	if JobID != None:
		try:
			sourceNumber = int(sourceName.replace('Source'))
			waitTime = sourceNumber * 20
		except:
			waitTime = random.random() * 600

		print "\nWaiting %i seconds before starting..." % waitTime
		time.sleep(waitTime)
		print "Proceeding."

		# Defind the temporary scratch directories
		Username = getpass.getuser()
		ScratchDirectoryJobSpecific = "/scratch/%s/%s" % (Username, JobID)
		ScratchDirectoryUserSpecific = "/scratch/%s" % (Username)
		ScratchDirectory = ScratchDirectoryJobSpecific

		# Save the name of the original output directory
		OutputDirectory_NSF = OutputDirectory

		# Set the scratch directory to be the output directory
		OutputDirectory = ScratchDirectory


	# Create the output directory
	print "\nCreating custom output directory:\n%s" % OutputDirectory        
	cmd = "mkdir -p " + OutputDirectory
	os.system(cmd)

	# Extended source templates
	extendedSourcesDirectory = LikelihoodDirectory + '/ExtendedSources/Templates/'

	# Catalog files
	Catalog2FGL = LikelihoodDirectory + '/Catalogs/gll_psc_v08.fit'
	Catalog3FGL = LikelihoodDirectory + '/Catalogs/gll_psc_v16.fit'

	# Diffuse and isotropic models
	galpropModel = LikelihoodDirectory + '/DiffuseModels/template_4years_P8_V2_scaled_trim.fits'

	# Select the appropriate isotropic model
	if 'P8R2_SOURCE_V6' in irfs:
		isotropicModel = LikelihoodDirectory + '/DiffuseModels/iso_P8R2_SOURCE_V6_v06.txt'
	elif 'P8R2_TRANSIENT020_V6' in irfs:
		isotropicModel = LikelihoodDirectory + '/DiffuseModels/iso_P8R2_TRANSIENT020_V6_v06.txt'
	elif 'P8R2_TRANSIENT010_V6' in irfs:
		isotropicModel = LikelihoodDirectory + '/DiffuseModels/iso_P8R2_TRANSIENT010_V6_v06.txt'
	else:
		print '\n***Error: Unrecognized IRF ***'
		print 'Current options:'
		print 'P8R2_TRANSIENT020'
		print 'P8R2_TRANSIENT010'		
		print 'P8R2_SOURCE_V6'
		print '\nExiting.'

	# Alternative diffuse models
	# galpropModel = '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/gll_iem_v06.fits'
	# isotropicModel = '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/iso_P8R2_SOURCE_V6_v06.txt'


	# Move into the output directory 
	os.chdir(OutputDirectory)

	# Define the pfile directory
	if JobID == None:
		PFILESDirectory = "%s/pfiles_%s/" % (OutputDirectory, sourceName)	
	else:
		PFILESDirectory = "%s/pfiles/" % ScratchDirectory

	# Remove any pre-existing pfiles
	if(os.path.isdir(PFILESDirectory)==True):
		cmd = "rm -r %s" % PFILESDirectory
		os.system(cmd)			

	# Set the new pfiles directory
	SetPfilesDirectory(PFILESDirectory)

	# Import the necessary gtapps	
	gtselect = GtApp('gtselect')
	gtmktime = GtApp('gtmktime')
	gtexpmap = GtApp('gtexpmap')
	gtbin = GtApp('gtbin')
	gtltcube = GtApp('gtltcube')
	gtexpcube2 = GtApp('gtexpcube2')
	gtdiffrsp = GtApp('gtdiffrsp')
	gtlike = GtApp('gtlike')
	gttsmap = GtApp('gttsmap')
	gtsrcmaps = GtApp('gtsrcmaps')
	gtmodel = GtApp('gtmodel')
	gtfindsrc = GtApp('gtfindsrc')

	# Setup the data type
	if 'P8R2_TRANSIENT010_V6' in irfs:
		dataclass_FT1 = 'P8_P302_ALL'		# Data class
		dataclass_FT2 = 'P8_P302_ALL'		# Data class
		evclass = 64						# Source Class	
		evclassname = 'Transient'			# The base event class. Source class can be filtered later using evclass in gtselect
		evclsmin = 0 						# Transient Class		
		evclsmax = 64 						# Source Class
		evtype = "INDEF"					# Event Type

	if 'P8R2_TRANSIENT020_V6' in irfs:
		dataclass_FT1 = 'P8_P302_ALL'		# Data class
		dataclass_FT2 = 'P8_P302_ALL'		# Data class
		evclass = 16						# Source Class	
		evclassname = 'Transient'			# The base event class. Source class can be filtered later using evclass in gtselect
		evclsmin = 0 						# Transient Class		
		evclsmax = 16 						# Source Class
		evtype = "INDEF"					# Event Type

	elif 'P8R2_SOURCE_V6' in irfs:
		dataclass_FT1 = 'P8_P302_BASE'		# Data class
		dataclass_FT2 = 'P8_P302_BASE'		# Data class
		evclass = 128						# Source Class
		evclassname = 'Source'				# The base event class. Source class can be filtered later using evclass in gtselect
		evclsmin = 0 						# Transient Class		
		evclsmax = 128 						# Source Class
		evtype = "INDEF"					# Event Type

	else:
		print '\n***Error: Unrecognized IRF ***'


	# Setup data extraction parameters
	radius = ROI						# Degrees		
	astroServerRadius = 30 				# Degrees
	binsize = 0.25						# Degrees
	dt = 86400.0 						# Light curve binning (1 day)

	# Specify the radius outside of which the xml model is fixed to the catalog values
	if fixedModel == True:
		radLim = 0.001					# Degrees
	else:
		radLim = 10						# Degrees


	# Define the working files
	ft1file = "%s/ft1_%s.fits" % (OutputDirectory, sourceName)
	ft2file = "%s/ft2_%s.fits" % (OutputDirectory, sourceName)	
	filteredEvents = '%s/ft1_filtered_%s.fits' % (OutputDirectory, sourceName)
	filteredEventsGTI = '%s/ft1_filteredGTI_%s.fits' % (OutputDirectory, sourceName)
	filteredEventsMKT = '%s/ft1_filteredMKT_%s.fits' % (OutputDirectory, sourceName)
	cmap = '%s/cmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	cmap_HighResolution	= '%s/cmap_HighResolution_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	ccube = '%s/ccube_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	ltcube = '%s/ltcube_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	expmap = '%s/expmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	expmapFlat = '%s/expmapFlat_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	bexpmap = '%s/bexpmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	bexpmap_allSky = '%s/bexpmap_allSky_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	srcmap = '%s/srcmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	likelihoodResults = '%s/likelihoodResults_%s_%s_%s.txt' % (OutputDirectory, tmin, tmax, sourceName)
	tsmap = '%s/tsmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	tsmapPNG = '%s/tsmap_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	regionFile = '%s/3FGL.reg' % (OutputDirectory)		
	modelMap = '%s/modelmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	gtfindsrcResults = '%s/gtfindsrc_%s_%s_%s.txt' % (OutputDirectory, tmin, tmax, sourceName)
	modelMapPlot = '%s/modelmap_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	modelMapPlotAnnotated = '%s/modelmapAnnotated_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	lightCurvePlot = '%s/lightcurve_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	lightCurveData = '%s/lightcurve_%s_%s_%s.txt' % (OutputDirectory, tmin, tmax, sourceName)	
	xmlModel = '%s/Model_%s.xml' % (OutputDirectory, sourceName)
	xmlModelFirstPass = '%s/Model_%s_FisrtPass.xml' % (OutputDirectory, sourceName)	
	xmlModelFit = '%s/Model_%s_Fit.xml' % (OutputDirectory, sourceName)
	xmlModelFixed = '%s/Model_%s_Fixed.xml' % (OutputDirectory, sourceName)	
	xmlModel_WithoutCandidateSource = '%s/Model_%s_WithoutCandidateSource.xml' % (OutputDirectory, sourceName)
	xmlModelFixed_WithoutCandidateSource = '%s/ModelFixed_%s_WithoutCandidateSource.xml' % (OutputDirectory, sourceName)
	xmlModelFitMaxTS = "%s/Model_%s_MaxTS.xml" % (OutputDirectory, sourceName)
	likelihoodResults = "%s/likelihoodResults_%s.txt" % (OutputDirectory, sourceName)
	recursiveLikelihoodResults = "%s/recursiveLikelihoodFit_%s_%s_%s.txt" % (OutputDirectory, tmin, tmax, sourceName)
	summaryMapTS = "%s/summaryMapTS_%s_%s_%s.png" % (OutputDirectory, tmin, tmax, sourceName)
	summaryMapFluxRatio = "%s/summaryMapFluxRatio_%s_%s_%s.png" % (OutputDirectory, tmin, tmax, sourceName)
	apertureLightCurvePlot = "%s/apertureLightCurvePlot_%s_%s_%s.png" % (OutputDirectory, tmin, tmax, sourceName)
	apertureLightCurvePlotNormalized = "%s/apertureLightCurvePlotNormalized_%s_%s_%s.png" % (OutputDirectory, tmin, tmax, sourceName)
	likelihoodFitPlot = '%s/likelihoodFit_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)

	# Get the data from the astroserver at SLAC
	if getData == True:

		# Setup the FT1 parameters
		FT1_options = {'event-sample' : dataclass_FT1,
				   'output-ft1' : ft1file,
				   'minEnergy' : emin,
				   'maxEnergy' : emax,
				   'minTimestamp' : tmin,
				   'maxTimestamp' : tmax,
				   'ra' : ra,
				   'dec' : dec,
				   'radius' : astroServerRadius}

		# Query the astroserver and get the FT1 file		  
		print '\nGetting the FT1 file:' 
		astro_query(FT1_options, store='store', verbose=True)


		# Setup the FT2 parameters
		FT2_options = {'event-sample' : dataclass_FT2,
				   'output-ft2-30s' : ft2file,
				   'minTimestamp' : tmin,
				   'maxTimestamp' : tmax}

		# Query the astroserver and get the FT2 file
		print '\nGetting the FT2 file:' 			   
		astro_query(FT2_options, store='storeft2', verbose=True)

	# Generating the neccessary working files
	if generateFiles == True:

		# Select the photons
		print '\nSelecting the photons:'
		gtselect.run(infile=ft1file,
					outfile=filteredEventsMKT,
					ra=ra, dec=dec, rad=radius,
					emin=emin, emax=emax,
					tmin=tmin, tmax=tmax,
					zmax=zmax, evclass=evclass, 
					evclsmax=evclsmax, convtype=-1,
					phasemin=0.0, phasemax=1.0)


		# Check to see if any photons survived the cut			 
		ft1 = pyfits.open(filteredEventsMKT)
		if ft1['EVENTS'].header['NAXIS2'] == 0:
			print 'No photons survived the gtmktime cut'
			return		


		# Select the good time intervals
		print '\nSelecting the good time intervals:'		
		gtmktime.run(scfile=ft2file,
					filter="IN_SAA!=T && LIVETIME>0 && (ANGSEP(RA_ZENITH,DEC_ZENITH,%s,%s) + 12.0 < 105.0)" % (ra, dec),
					#filter='DATA_QUAL>0 && LAT_CONFIG==1',
					roicut='yes',
					evfile=filteredEventsMKT,
					outfile=filteredEvents,
					overwrite='yes')



		# Check to see if any photons survived the cut			 
		ft1 = pyfits.open(filteredEvents)
		if ft1['EVENTS'].header['NAXIS2'] == 0:
			print 'No photons survived the gtmktime cut'
			return				


		# Generate the livetime cube
		print '\nGenerating the livetime cube:'		
		gtltcube.run(evfile=filteredEvents,
					scfile=ft2file,
					outfile=ltcube,
					dcostheta=0.025,
					binsz=1)


		# Generate an exposure map.  This is for unbinned likelihood only!	
		print '\nGenerating the exposure map:'				
		gtexpmap.run(evfile=filteredEvents,
					scfile=ft2file,
					expcube=ltcube,
					outfile=expmap,
					irfs=irfs,
					srcrad=astroServerRadius,
					nlong=88, nlat=88, nenergies=10, coordsys='CEL')
					# nlong=120, nlat=120, nenergies=37, coordsys='CEL')


		# Generate an exposure map.  This is for unbinned likelihood only!	
		# print '\nGenerating the exposure map:'				
		# gtexpmap.run(evfile=filteredEvents,
		# 			scfile=ft2file,
		# 			expcube=ltcube,
		# 			outfile=expmapFlat,
		# 			irfs=irfs,
		# 			srcrad=astroServerRadius,
		# 			nlong=120, nlat=120, nenergies=2, coordsys='CEL')	
		# 			#nlong=88, nlat=88, nenergies=10, coordsys='CEL')



	# Create the source model
	if makeModel == True:

		print '\nGenerating the xml model:'
		mymodel=make3FGLxml.srcList(Catalog3FGL, filteredEvents, xmlModel)
		mymodel.makeModel(galpropModel, 'GAL_V02', isotropicModel, 'EG_v02', radLim=radLim, extDir=extendedSourcesDirectory, psForce=False, makeRegion=True)

		# Determine the object of interest
		if 'none' in association:

			# Save a copy of this model before adding a candidate source 
			cmd = "cp %s %s" % (xmlModel, xmlModel_WithoutCandidateSource)
			os.system(cmd)	

			# Add a condidate point source at the inital FAVA position (regardless of whether there is an association)
			sourceOfInterest = 'CandidateSource'
			AddCandidateSource(float(ra), float(dec), xmlModel, fixIndex=fixIndex)

		else:

			# Generate the source name
			sourceOfInterest = association

			# Check to see if our associated source is in the xml model.  If not, add a candidate point source.
			SourceRA, SourceDec = ExtractCoordinates(xmlModel, sourceOfInterest)

			# Add a point source to the xml model if there is no association
			if SourceRA == None:
				print '\nCould not locate %s in the xml model.' % sourceOfInterest
				return


	# Generating the neccessary working files
	if generateFiles == True:

		# Make a counts map from the event data
		print '\nCreating the counts map:'
		gtbin.run(evfile=filteredEvents,
					scfile=ft2file,
					outfile=cmap,
					algorithm='CMAP',
					nxpix=160, nypix=160, binsz=binsize, coordsys='CEL',
					#nxpix=200, nypix=200, binsz=0.2, coordsys='CEL',
					#nxpix=300, nypix=300, binsz=0.2, coordsys='CEL',					
					xref=ra, yref=dec, axisrot=0, proj='AIT')


		print '\nCreating the high resolution counts map:'
		gtbin.run(evfile=filteredEvents,
					scfile=ft2file,
					outfile=cmap_HighResolution,
					algorithm='CMAP',
					nxpix=300, nypix=300, binsz=0.1, coordsys='CEL',					
					xref=ra, yref=dec, axisrot=0, proj='AIT')


		# Make a counts cube from the event data
		print '\nCreating the counts cube:'
		gtbin.run(evfile=filteredEvents,
					scfile=ft2file,
					outfile=ccube,
					algorithm='CCUBE',
					nxpix=160, nypix=160, binsz=binsize, coordsys='CEL',
					#nxpix=200, nypix=200, binsz=0.2, coordsys='CEL',
					xref=ra, yref=dec, axisrot=0, proj='AIT',
					emin=emin, emax=emax, enumbins=30)


		# Generate an exposure map.  This is for binned likelihood only, but is also used in the srcmap generation
		print '\nGenerating the binned exposure map:'				
		gtexpcube2.run(infile=ltcube,
					cmap=ccube,
					outfile=bexpmap,
					irfs=irfs,
					nxpix=400, nypix=400, binsz=binsize, coordsys='CEL',
					xref=ra, yref=dec, axisrot=0, proj='AIT',
					emin=emin, emax=emax, nenergies=30)



		# Generate an all sky exposure map.  This is for binned likelihood only, but is also used in the srcmap generation
		if makeModelMap == True or makeSourceMap == True:
			print '\nGenerating the all sky binned exposure map:'				
			gtexpcube2.run(infile=ltcube,
						cmap=ccube,
						outfile=bexpmap_allSky,
						irfs=irfs,
						nxpix=int(360/binsize), nypix=int(180/binsize), binsz=binsize, coordsys='CEL',
						xref=ra, yref=dec, axisrot=0, proj='AIT',
						emin=emin, emax=emax, nenergies=30)


	# Quite here if the user only wants the data and associated products
	if justGetData == True:
		print 'Done.\n'
		return

	# Continue generating the neccessary working files
	if generateFiles == True:

		if skipDiffuseResponse == False:
			# Compute the diffuse response
			print '\nComputing the diffuse response:'				
			gtdiffrsp.run(evfile=filteredEvents,
						scfile=ft2file,
						srcmdl=xmlModel,
						irfs=irfs,
						#evclass=evclass,
						clobber='yes',
						convert='yes')
		else:
			print '\nSkipping diffuse response'





	# Run an likelihood fit at the initial FAVA position
	if performLikelihoodFit == True:

		# print '\n\nPerforming the likelihood fit...'
		print ''
		print "\nFitting xml model:\n%s" % xmlModel

		try:
			# Setup the unbinned likelihood object for the first pass
			print '\nCreating the likelihood object (first pass)...'
			if 'UNBINNED' in statistic:
				obs = UnbinnedObs(filteredEvents, ft2file, expMap=expmap, expCube=ltcube, irfs=irfs)
				like = UnbinnedAnalysis(obs, xmlModel, optimizer='MINUIT')
			else:				
				obs = BinnedAnalysis.BinnedObs(filteredEvents, ft2file, expMap=bexpmap, expCube=ltcube, irfs=irfs)
				like = BinnedAnalysis.binnedAnalysis(obs, xmlModel, optimizer='MINUIT')				
			print 'Done.'

			# Set the tolerance type to "relative"
			like.tolType = 1

			# Print out the likelihood parameters
			print '\nLikelihood Fit Parameters:'
			print like
			print 'Statistic: %s' % statistic
			print 'Fixed model: %s' % fixedModel

			# Set the optimizer
			if 'NEWMINUIT' in optimizer:
				optObject = pyLike.NewMinuit(like.logLike)
			else:
				optObject = pyLike.Minuit(like.logLike)

			# Performing the fit
			print '\nFitting likelihood model...'
			logL = like.fit(verbosity=1, covar=True, tol=0.001, optObject=optObject)
			print 'Done.\n'
			print "%s Return Code: %s" % (optimizer, optObject.getRetCode())			
			print 'logL = %s' % logL				

			# Write out the first pass model
			like.logLike.writeXml(xmlModel)

			# Plot the results
			#like.plot()


			if performRefinedFit == True:

				print '\nCreating the likelihood object (second pass)...'
				if 'UNBINNED' in statistic:
					obs = UnbinnedObs(filteredEvents, ft2file, expMap=expmap, expCube=ltcube, irfs=irfs)
					like = UnbinnedAnalysis(obs, xmlModel, optimizer=optimizer)
				else:				
					obs = BinnedAnalysis.BinnedObs(filteredEvents, ft2file, expMap=bexpmap, expCube=ltcube, irfs=irfs)
					like = BinnedAnalysis.binnedAnalysis(obs, xmlModel, optimizer=optimizer)				
				print 'Done.'

				# Set the tolerance type to "relative"
				like.tolType = 1

				# Print out the likelihood parameters
				print '\nLikelihood fit parameters:'
				print like

				# Set the optimizer
				if 'NEWMINUIT' in optimizer:
					optObject = pyLike.NewMinuit(like.logLike)
				else:
					optObject = pyLike.Minuit(like.logLike)
			

				# Performing the fit
				print '\nFitting likelihood model...'
				logL = like.fit(verbosity=1, covar=True, tol=0.1, optObject=optObject)
				print 'Done.\n'
				print "%s Return Code: %s" % (optimizer, optObject.getRetCode())			
				print 'logL = %s' % logL			


			# Write out the second pass model
			like.logLike.writeXml(xmlModel)


			sourceDetails = {}
			print "\nSource: Final TS"
			for source in like.sourceNames():
				ts_source = like.Ts(source)
				print "%s: %s" % (source, ts_source)
				sourceDetails[source] = ts_source


			if removeWeakSources  == True:

				# Delete unconstrained sources and refit
				print "\nDeleting unconstrained sources..."

				for source, ts_source in sourceDetails.iteritems():
					if (ts_source < 2):
						if sourceOfInterest not in source:
							print "Deleting %s" % source
							like.deleteSource(source)

				# Set the optimizer
				if 'NEWMINUIT' in optimizer:
					optObject = pyLike.NewMinuit(like.logLike)
				else:
					optObject = pyLike.Minuit(like.logLike)

				# Performing the fit
				print '\nRe-fitting likelihood model...'
				logL = like.fit(verbosity=1, covar=True, tol=0.1, optObject=optObject)
				print 'Done.\n'
				print "%s Return Code: %s" % (optimizer, optObject.getRetCode())			
				print 'logL = %s' % logL	


			# Go back and extract the TS and flux values for each of the fit sources
			print '\nExtracting fit parameters...'
			ts = like.Ts(sourceOfInterest)
			flux = like.flux(sourceOfInterest)

			# Display the best fit parameters if ts >= tsMin, otherwise calculate the upper limits
			if float(ts) >= float(tsMin): 

				# Extract the likelihood fit results for the source of interest
				try:
					print like.model[sourceOfInterest]
					TSValue = like.Ts(sourceOfInterest)
				except:
					TSValue = 'NA'

				# Extract the photon flux
				try:
					PhotonFlux = "%.2e" % like.flux(sourceOfInterest, emin=emin, emax=emax)
					PhotonFluxError = "%.2e" % like.fluxError(sourceOfInterest, emin=emin, emax=emax)
				except:
					PhotonFlux = 'NA'
					PhotonFluxError = 'NA'

				# Extract the photon index
				try:
					Index = like.par_index(sourceOfInterest, 'Index')
					PhotonIndex = "%.2f" % like[Index].value()
					PhotonIndexError = "%.2f" % like[Index].error()
				except:
					PhotonIndex = 'NA'
					PhotonIndexError = 'NA'

				# Print the maximum likelihood results
				print "\nTS: %s" % like.Ts(sourceOfInterest)
				print "Ra Dec: %.3f %.3f" % (float(ra), float(dec))
				print "Photon Index: %s +/- %s" % (PhotonIndex, PhotonIndexError)
				print "Photon Flux: %s +/- %s" % (PhotonFlux, PhotonFluxError)			

				# Save the likelihood results
				output = open(likelihoodResults, 'w')
				output.write("Likelihood Results\n")
				output.write("TS: %s\n" % like.Ts(sourceOfInterest))
				output.write("RaDec: %.3f %.3f +/- %.3f\n" % (float(ra), float(dec), 0))
				output.write("PhotonIndex: %s +/- %s\n" % (PhotonIndex, PhotonIndexError))
				output.write("PhotonFlux: %s +/- %s\n" % (PhotonFlux,PhotonFluxError))
				output.close()

			else:

				# Profile likelihood upper limits (?)
				# ul=UpperLimits(like)
				# ul[sourceOfInterest].compute()
				# print ul[sourceOfInterest].results

				# Bayesian upper limits
				# upper_limit95, results = IUL.calc_int(like, sourceOfInterest, cl=0.95, freeze_all=True)				
				photonFluxUpperLimit95, results = IUL.calc_int(like, sourceOfInterest, cl=0.95)

				# Convert the photon flux upper limit to an energy flux upper limit
				energyFluxUpperLimit95 = computeEnergyFlux(photonFluxUpperLimit95, '-2.1', emin=float(emin), emax=float(emax))

				

				# Print the maximum likelihood results
				print "\nTS: %s" % like.Ts(sourceOfInterest)
				print "Ra Dec: %.3f %.3f" % (float(ra), float(dec))
				print "Photon Flux Upper Limit: %s ph cm-2 s-1" % photonFluxUpperLimit95
				print "Energy Flux Upper Limit (index=-2.1): %s ph cm-2 s-1" % energyFluxUpperLimit95

				# Save the likelihood results
				print '\nSaving results to:\n%s' % likelihoodResults
				output = open(likelihoodResults, 'w')
				output.write("Likelihood Results\n")
				output.write("TS: %s\n" % like.Ts(sourceOfInterest))
				output.write("RaDec: %.3f %.3f +/- %.3f\n" % (float(ra), float(dec), 0))
				output.write("PhotonIndex: %s +/- %s\n" % (PhotonIndex, PhotonIndexError))
				output.write("PhotonFluxUpperLimit: %s\n" % photonFluxUpperLimit95)
				output.write("EnergyFluxUpperLimit: %s\n" % energyFluxUpperLimit95)				
				output.close()


			# Save the final xml file
			like.writeXml(xmlFile=xmlModel)



		except Exception, message:
			print message
			print traceback.format_exc()



	if plotFit == True:

		# Make sure matplotlib has been successfully imported
		import matplotlib.pylab as plot

		# Make the plot
		plot.figure(figsize=(10, 6.39))

		# Set the plot limits
		plot.ylim((1e-2,1e3))
		plot.xlim((emin,emax))

		# Set the plot labels
		plot.xlabel('Energy (MeV)')
		plot.ylabel('Counts')

		# Extract the energy range
		E = (like.energies[:-1] + like.energies[1:])/2.

		# Extract and plot the model contributions
		sum_model = numpy.zeros_like(like._srcCnts(like.sourceNames()[0]))

		numberOfSources = 0
		for sourceName in like.sourceNames():
			sum_model = sum_model + like._srcCnts(sourceName)
			plot.loglog(E,like._srcCnts(sourceName), label=sourceName[1:])
			numberOfSources = numberOfSources + 1

		plot.loglog(E,sum_model, label='Total Model')

		# Plot the errors
		y = like._Nobs()
		yerr = numpy.array(numpy.sqrt(like._Nobs()))
		yerr2 = yerr
		yerr2[yerr>=y] = y[yerr>=y]*.999999
		# plot.errorbar(E,like._Nobs(), yerr=numpy.sqrt(like._Nobs()), fmt='o', label='Counts')
		plot.errorbar(E, y, yerr=yerr, fmt='o', label='Counts')

		# Add a legend
		if numberOfSources < 30:
			plot.legend(loc=1, fontsize='xx-small')

		# Save the plot
		print '\nSaving likelihood fit plot to:\n%s' % likelihoodFitPlot
		plot.savefig(likelihoodFitPlot, dpi=72)



	if maketsmap == True:

		if 'UNBINNED' in statistic:
			countsMap_tsMap = cmap
		else:
			countsMap_tsMap = ccube

		# Generate a tsmap
		print "\nCalling gttsmap..."
		gttsmap.run(evfile=filteredEvents,
					scfile=ft2file,
					expmap=expmap,
					expcube=ltcube,
					cmap=countsMap_tsMap,
					bexpmap=bexpmap,
					srcmdl=xmlModel_WithoutCandidateSource,
					outfile=tsmap,
					irfs=irfs,
					statistic=statistic,
					optimizer=optimizer,
					#ftol=1e-5,
					nxpix=int(float(dra)/float(tsmapBinSize)), nypix=int(float(ddec)/float(tsmapBinSize)), binsz=tsmapBinSize, coordsys='CEL',
					xref=ra, yref=dec, proj='AIT')

		# Plot the ts map
		fig=aplpy.FITSFigure(tsmap)
		fig.show_colorscale(cmap = "jet", stretch='linear')

		# Add a title
		fig.add_label(0.5, 1.05, sourceName, relative=True, size='large')

		# Add a color
		fig.add_colorbar()
		fig.colorbar.set_location('right')
		fig.colorbar.set_width(0.2)  
		fig.colorbar.set_pad(0.05) 
		fig.add_grid()
		fig.axis_labels.set_font(size='large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')

		# Tick labels:
		fig.tick_labels.show()
		fig.tick_labels.show_x()
		fig.tick_labels.show_y()
		fig.tick_labels.set_xformat('ddd:dd')
		fig.tick_labels.set_yformat('ddd:dd')

		if os.path.isfile(regionFile):

			# Get the center of the image. Temporary fix for display regions, since we draw labels only ifsource are inside the image
			hdulist = pyfits.open(tsmap)
			xpix_size=abs(float(hdulist[0].header['CDELT1']))
			xpix_ref=float(hdulist[0].header['CRPIX1'])
			xpix_ref_x=float(hdulist[0].header['CRVAL1'])
			npix_x=float(hdulist[0].header['NAXIS1'])
			ypix_size=abs(float(hdulist[0].header['CDELT2']))
			ypix_ref=float(hdulist[0].header['CRPIX2'])
			ypix_ref_y=float(hdulist[0].header['CRVAL2'])
			npix_y=float(hdulist[0].header['NAXIS2'])
			xmin=xpix_ref_x-xpix_size*(xpix_ref)
			xmax=xpix_ref_x+xpix_size*(npix_x-xpix_ref)
			ymin=ypix_ref_y-ypix_size*(ypix_ref)
			ymax=ypix_ref_y+ypix_size*(npix_y-ypix_ref)
			hdulist.close()

			#read the ds9 region file manually and plot the stuff
			#see http://matplotlib.org/api/markers_api.html#module-matplotlib.markers for marker styles
			print "superimposing ds9 region.."
			with open(regionFile, 'r') as file:
				for line in file:
					l=line.strip()
					if 'point' not in l:
						continue
					if ";" in l:
						coord=l.split("#")[0].split(";")[0]
						ra, dec=l.split("#")[0].split(";")[1].replace("point(", "").replace(")", "").split(",")
						name=l.split("{")[-1].replace(" ", "").replace("}", "")
						if coord!="J2000":
							print "pyapp.fitstopng: region support is experimental, only J2000 ra-dec is implemented."
							fig.save(outfile)
							return
					else:
						ra, dec=l.split("#")[0].replace("point(", "").replace(")", "").split(",")
						name=l.split("{")[-1].replace(" ", "").replace("}", "")
						xmin, xmax, ymin, ymax=0., 360., -90., 90.      #temporary hack for different format of region files

					# ll, bb = utils.equ2gal(float(ra), float(dec))
					#only marker inside the figure will be drawn
					# if (ll>xmin) and (ll<xmax) and (bb>ymin) and (bb<ymax):
					if (float(ra)>xmin) and (float(ra)<xmax) and (float(dec)>ymin) and (float(dec)<ymax):
						print "adding source", name, " at (ra=%f, dec=%f) deg" %(ra, dec) 
						fig.show_markers(ra, dec, edgecolor='black', facecolor='white', marker='x', alpha=1, s=40, linewidths=1.5)
						fig.add_label(ra, dec-0.2, name, color='black')
					#                    else:
					#                        print "point (%f, %f) Galcoord is outside the map" %(ll, bb)
					#                        print "xmin, xmax, ymin, ymax", xmin, xmax, ymin, ymax

		print "computing 95% (2sigma) contour"
		#compute contours levels
		maxts=numpy.amax(numpy.array(fig._data))
		loglikestep=5.99        ## 95% 2 sigma. for p=70% loglikestep=2.41 #check here for info http://seal.web.cern.ch/seal/documents/minuit/mnerror.pdf 
		contour_levels=[maxts-loglikestep]  #contlvls=[maxts-2.41]
		lstyles=["dashdot"]  #, "dashed", "dashdot", "solid", "dotted"]
		fig.show_contour(tsmap, levels=contour_levels, colors="black", linestyles=lstyles, interpolation='cubic')

		### save and quit
		fig.save(tsmapPNG)
		fig.close()
		print "TSMap image saved to", tsmapPNG



	# Generate a source map of the best fit likelihood model
	if makeSourceMap == True:

		print '\nGenerating a source map:'
		try:
			gtsrcmaps.run(scfile=ft2file,
							cmap=ccube,
							expcube=ltcube,
							srcmdl=xmlModel,
							outfile=srcmap,
							bexpmap=bexpmap_allSky,
							irfs=irfs,
							emapbnds='no')
		except Exception, message:
			print message
			print traceback.format_exc()

	if makeModelMap == True:

		# Generate a source map of the best fit likelihood model
		print '\n\nGenerating a source map:'
		try:
			gtsrcmaps.run(scfile=ft2file,
							cmap=ccube,
							expcube=ltcube,
							srcmdl=xmlModelBest,
							outfile=srcmap,
							bexpmap=bexpmap_allSky,
							irfs=irfs,
							emapbnds='no')
		except Exception, message:
			print message
			print traceback.format_exc()


		# Generate model map of the best fit likelihood model	
		print '\nGenerating a model map:'	
		try:	
			gtmodel.run(srcmaps=srcmap,
							srcmdl=xmlModelBest,
							outfile=modelMap,
							irfs=irfs,
							expcube=ltcube,
							bexpmap=bexpmap_allSky)
		except Exception, message:
			print message
			print traceback.format_exc()
		

		# Plot the model map without catalog annotations
		try:

			import matplotlib
			import matplotlib.pylab as plot
			import pyfits
			import traceback
			from matplotlib.ticker import FuncFormatter

			fits = pyfits.open(modelMap)
		
			# Extract the wcs data from the fits file
			xReferencePixel = fits[0].header['CRPIX1']
			xReferencePixelRA = fits[0].header['CRVAL1']
			xIncrementPerPixel = fits[0].header['CDELT1']
			yReferencePixel = fits[0].header['CRPIX2']
			yReferencePixelDEC = fits[0].header['CRVAL2']
			yIncrementPerPixel = fits[0].header['CDELT2']

			# Calculate the extent of the image
			xmin = xReferencePixelRA - (xReferencePixel* xIncrementPerPixel)
			xmax = xReferencePixelRA + (xReferencePixel* xIncrementPerPixel)
			ymin = yReferencePixelDEC - (yReferencePixel* yIncrementPerPixel)
			ymax = yReferencePixelDEC + (yReferencePixel* yIncrementPerPixel)
			xRange = [xmin,xmax]
			yRange = [ymin,ymax]

			# Make sure that we don't have any ra values below zero or greater than 360, they should wrap ra instead.
			for i in range(len(xRange)):
				if xRange[i] < 0:
					xRange[i] = xRange[i] + 360.0
				if xRange[i] > 360:
					xRange[i] = xRange[i] - 360.0

			# Make sure that we don't have any dec values below or above +/- 90, they should instead wrap in both ra and dec.
			for i in range(len(yRange)):
				if yRange[i] < -90:
					yRange[i] = ((yRange[i] + 90) + 90)*-1
					#xRange[i] = xRange[i] + 180.0
				if yRange[i] > 90:
					yRange[i] = 90 - (yRange[i] - 90)
					#xRange[i] = xRange[i] + 180


			# Extract the model map
			array = fits[0].data

			# Explicitly create the figure.  This helps with forcing a specific aspect ratio later.
			fig = plot.figure()
			ax = fig.add_subplot(111)		

			# Plot the model map using a Lambert Azimuthal Equal Area Projection
			sys.path.append("/nfs/slac/g/ki/ki08/kocevski/LATBA/lib/python_rhel6-64/")
			from mpl_toolkits.basemap import Basemap

			# Create a base map on which to plot the results
			m = Basemap(height=4.3e6,width=4.3e6, projection='laea', lon_0 = float(ra)*-1, lat_0 = float(dec), resolution ='l',area_thresh=1000., celestial=True)

			# plot the image map
			m.imshow(numpy.log10(array), origin='lower', cmap=matplotlib.cm.get_cmap('seismic'), )

			# Setup the map grid
			m.drawmapboundary(fill_color='#ffffff')
			m.drawparallels(numpy.arange(-90,95,5),labels=[1,0,0,0], fmt=customDecLabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
			m.drawmeridians(numpy.arange(0,365,5),labels=[0,0,0,1],fmt=customRALabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
			m.suppress_ticks = False

			# Force the aspect ratio to be 1:1
			try:
				# Force the aspect ratio to be 1:1
				forceAspect(ax,aspect=1)
			except Exception, message:
				print message
				print traceback.format_exc()

			# Annotate the plot with the original FAVA location
			x_FAVA, y_FAVA = m(float(ra), float(dec))
			#m.scatter(x_FAVA, y_FAVA, marker='+',color='skyblue', alpha=1, clip_on=True)
			#m.scatter(x_FAVA, y_FAVA, marker='+', s=75, facecolors='none', edgecolors='w')

			# Generate the annotation label
			if 'CandidateSource' in sourceOfInterest:
				annotationLabel = 'Candidate Source'
			else:
				annotationLabel = association

			# Convert the max ts ra and dec to map coordinates
			(mMaxRA, mMaxDec) = m(float(ra_sourceOfInterest)+0.5,float(dec_sourceOfInterest)+0.5)

			# Add the candidate source annotation
			plot.annotate(annotationLabel, xy=(mMaxRA,mMaxDec), xytext=(-25,25), textcoords='offset points', ha='center', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='black'))

			# Setup the plot
			plot.xlabel('RA')
			plot.ylabel('Dec')
			plot.gca().xaxis.labelpad = 20
			plot.gca().yaxis.labelpad = 20

			# Save the model map
			print "\nSaving model map plot to: %s" % modelMapPlot
			plot.savefig(modelMapPlot, bbox_inches='tight', dpi=100)
			plot.close()

		except Exception, message:
			print message
			print traceback.format_exc()


		# Plot the model map with catalog annotations
		try:

			fig = plot.figure(figsize=(6.2,6.0))
			fig.clip_on = True
			ax = fig.add_subplot(111)	

			# Create a base map on which to plot the results
			m = Basemap(height=4.3e6,width=4.3e6, projection='laea', lon_0 = float(ra)*-1, lat_0 = float(dec), resolution ='l',area_thresh=1000., celestial=True)

			# plot the image map
			m.imshow(numpy.log10(array), origin='lower', cmap=matplotlib.cm.get_cmap('seismic'))

			# Setup the map grid
			m.drawmapboundary(fill_color='#ffffff')
			m.drawparallels(numpy.arange(-90,95,5),labels=[1,0,0,0], fmt=customDecLabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
			m.drawmeridians(numpy.arange(0,365,5),labels=[0,0,0,1],fmt=customRALabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
			m.suppress_ticks = False

			# Force the aspect ratio to be 1:1
			try:
				# Force the aspect ratio to be 1:1
				forceAspect(ax,aspect=1)
			except Exception, message:
				print message
				print traceback.format_exc()

			# Setup the plot
			plot.xlabel('RA')
			plot.ylabel('Dec')
			plot.gca().xaxis.labelpad = 20
			plot.gca().yaxis.labelpad = 20

			# Extract 2FGL sources
			try:
				print '\nAdding 2FGL sources...'
				LAT2FGL = pyfits.open('/nfs/slac/g/ki/ki08/kocevski/LATBA/Catalogs/2FGL.fits')

				print '\nAdding monitored sources...'
				MonitoredSourceList = pyfits.open('/nfs/slac/g/ki/ki08/kocevski/LATBA/Catalogs/MonitoredSourceList.fits')
				
				print '\nAdding ATel sources...'
				Name_ATels = numpy.array([])
				Name_2FGL_ATels = numpy.array([])
				RA_ATels = numpy.array([])
				DEC_ATels = numpy.array([])
				ATelNumber = numpy.array([])
				ATelCatalog = '/nfs/slac/g/ki/ki08/kocevski/LATBA/Catalogs/ascd_atels_feb142013.dat'
				for line in fileinput.input([ATelCatalog]):
					if 'RA (J2000.0)' not in line:
						LineContents = line.split('|')
						Name_ATels = numpy.append(Name_ATels,LineContents[1].strip())
						Name_2FGL_ATels = numpy.append(Name_2FGL_ATels,LineContents[2].strip())
						RA_ATels = numpy.append(RA_ATels,float(LineContents[3].strip()))
						DEC_ATels = numpy.append(DEC_ATels,float(LineContents[4].strip()))
						ATelNumber = numpy.append(ATelNumber,LineContents[5].strip())
				print 'Done.'
				fileinput.close()
				
			except Exception, message:
				print message
				print traceback.format_exc()

			# 2FGL Sources
			RA_2FGL = LAT2FGL[1].data.RAJ2000
			DEC_2FGL = LAT2FGL[1].data.DEJ2000
			Flux_2FGL = LAT2FGL[1].data.Flux100_300
			SourceName_2FGL = LAT2FGL[1].data.Source_Name
			ASSOC1 = LAT2FGL[1].data.ASSOC1

			# Determine the best 2FGL source name
			BestName_2FGL = SourceName_2FGL
			#BestName_2FGL = numpy.array([])
			#for x,y in zip(ASSOC1, SourceName_2FGL):
			#		if len(x) > 0:
			#				BestName_2FGL = numpy.append(BestName_2FGL, x)
			#				pass
			#		else:
			#				BestName_2FGL = numpy.append(BestName_2FGL, y)
			#				pass
			#		pass

			# Extract the monitored sources
			Name_MonitoredSources = MonitoredSourceList[1].data['NAME']
			RA_MonitoredSources = MonitoredSourceList[1].data['RA']
			DEC_MonitoredSources = MonitoredSourceList[1].data['DEC']

			# Keep only the unique entries
			Name_MonitoredSourcesUnique, index = numpy.unique(Name_MonitoredSources,return_index=True)
			RA_MonitoredSourcesUnique = RA_MonitoredSources[index]
			DEC_MonitoredSourcesUnique = DEC_MonitoredSources[index]

			# Combine all sources
			RA_Sources = numpy.concatenate((RA_2FGL,RA_MonitoredSourcesUnique,RA_ATels))
			DEC_Sources = numpy.concatenate((DEC_2FGL,DEC_MonitoredSourcesUnique,DEC_ATels))
			Name_Sources = numpy.concatenate((BestName_2FGL,Name_MonitoredSourcesUnique,Name_ATels))

			# Remove any duplicate sources
			Names_SourcesUnique, index = numpy.unique(Name_Sources, return_index=True)
			RA_SourcesUnique = RA_Sources[index]
			DEC_SourcesUnique = DEC_Sources[index]

			# Remove individual sources
			DuplicateSources1 = ['2FGL J1745.6-2858','BL LAC', 'CGRABS J1848+3219', 'CGRABS J1849+6705','CGRABS J0211+1051','LS I+61 303','Mkn 421','PKS 0727-11','PKS 1424-41']
			DuplicateSources2 = ['PKS 1510-08','S5 1803+78','SUN','0FGLJ1641.4+3939','0235+164','S5 0716+71','GB6 J0742+5444','0827+243','PKS B0906+015','0FGL J0910.2-5044']
			DuplicateSources3 = ['1150+497','TON 599','PKS B1222+216','4C +21.35','J123939+044409','PSRB1259-63','1510-089','FERMI J1532-1321','PKS B 1622-297','4C +38.41']
			DuplicateSources4 = ['1633+382','J1717-5156','1730-130','3EG J2033+4118']
			DuplicateSources = DuplicateSources1 + DuplicateSources2 + DuplicateSources3 + DuplicateSources4
			Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique = RemoveSources(DuplicateSources, Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique)

			# Add the sources
			sort = [numpy.argsort(RA_SourcesUnique)]
			RA_SourcesUnique = RA_SourcesUnique[sort]
			DEC_SourcesUnique = DEC_SourcesUnique[sort] 
			Names_SourcesUnique = Names_SourcesUnique[sort]

			# Annotate the plot with the original FAVA location
			x_FAVA, y_FAVA = m(float(ra), float(dec))
			#m.scatter(x_FAVA, y_FAVA, marker='+',color='skyblue', s=50, alpha=1, clip_on=True)
			m.scatter(x_FAVA, y_FAVA, marker='+', s=75, facecolors='none', edgecolors='w')

			# Annotate the plot with the source names
			x_Sources, y_Sources = m(RA_SourcesUnique, DEC_SourcesUnique)
			m.scatter(x_Sources, y_Sources, marker='x',color='skyblue', alpha=1, clip_on=True)
			for i,j,k in zip(x_Sources, y_Sources, Names_SourcesUnique):
				plot.text(i+5e4,j,k,clip_on=True, alpha=1,size=8,color='skyblue')
				pass

			Annotate the plot with the sun position if the sun is within 25 degrees of ROI center
			try:
				SkyDir = getSunPosition( (float(tmin) + float(tmax)) / 2.0 )
				SunRa = SkyDir.ra()
				SunDec = SkyDir.dec()

				if maketsmap == True and MaxTS >= 25:
					SunSeperation = AngularSeperation( SunRa, SunDec, float(MaxRa), float(MaxDec))
				else:
					SunSeperation = AngularSeperation( SunRa, SunDec, float(ra_sourceOfInterest), float(dec_sourceOfInterest))

				if SunSeperation <= 15:
					x_Sun, y_Sun = m(SunRa, SunDec)
					m.scatter(x_Sun, y_Sun, marker='x',color='skyblue', alpha=1, clip_on=True)

					plot.text(x_Sun+5e4,y_Sun,'Sun',clip_on=True, alpha=1,size=8,color='skyblue')

			except Exception, message:
				print message				
				print traceback.format_exc()


			# Generate the annotation label
			if 'CandidateSource' in sourceOfInterest:
				annotationLabel = 'Candidate Source'
			else:
				annotationLabel = association

			# Convert the max ts ra and dec to map coordinates
			(mMaxRA, mMaxDec) = m(float(ra_sourceOfInterest)+0.5,float(dec_sourceOfInterest)+0.5)

			# Add the candidate source annotation
			plot.annotate(annotationLabel, xy=(mMaxRA,mMaxDec), xytext=(-25,25), textcoords='offset points', ha='center', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='black'))

			# Save the model map
			print "\nSaving annotated model map plot to: %s" % modelMapPlotAnnotated
			plot.savefig(modelMapPlotAnnotated, bbox_inches=matplotlib.transforms.Bbox([[0.28, 0.13], [5.67, 5.5]]), dpi=100)
			plot.close()

		except Exception, message:
			print message
			print traceback.format_exc()


	# Clean up
	print "\n\nCleaning up..."

	# Remove any dumped cores
	cores = glob.glob('core.*')

	if len(cores) > 0:
		cmd = "rm %s/core.*" % OutputDirectory
		os.system(cmd)	


	# Delete the pfiles directory
	if(os.path.isdir(PFILESDirectory)==True):
		cmd = "rm -R %s" % PFILESDirectory
		os.system(cmd)


	# Delete all files except those 
	if nuke == True:

		print '** Nuclear option initiated! **'
		print "rm %s/*.*" % OutputDirectory
		os.system("rm %s/*.*" % OutputDirectory)
		os.system("rmdir %s" % OutputDirectory)

		print '\n\nAnalysis Complete.\n'
		return

	elif cleanupAll == True:

		# Create a list of files to exclude from the deletion
		excludeList = [likelihoodResults]		

		# Delete all other files in the OutputDirectory
		files = glob.glob("%s/*" % OutputDirectory)
		for file in files:
			if file not in excludeList:
				print "rm %s" % file
				os.system("rm %s" % file)

	# Delete working files
	elif cleanup == True:

		#os.system("rm %s/ft1_*%s.*" % (OutputDirectory, sourceName))
		#os.system("rm %s/ft2_*%s.*" % (OutputDirectory, sourceName))
		if (os.path.isfile("%s/LC_*%s.*" % (OutputDirectory, sourceName))): os.system("rm %s/LC_*%s.*" % OutputDirectory)
		if (os.path.isfile("%s/expmap_*%s.*" % (OutputDirectory, sourceName))): os.system("rm %s/expmap_*%s.*" % OutputDirectory)
		if (os.path.isfile("%s/bexpmap_*%s.*" % (OutputDirectory, sourceName))): os.system("rm %s/bexpmap_*%s.*" % OutputDirectory)
		if (os.path.isfile("%s/ltcube_*%s.*" % (OutputDirectory, sourceName))): os.system("rm %s/ltcube_*%s.*" % OutputDirectory)
		if (os.path.isfile("%s/ccube_*%s.*" % (OutputDirectory, sourceName))): os.system("rm %s/ccube_*%s.*" % OutputDirectory)
		if (os.path.isfile("%s/cmap_*%s.*" % (OutputDirectory, sourceName))): os.system("rm %s/cmap_*%s.*" % OutputDirectory)
		if (os.path.isfile("%s/srcmap_*%s.*" % (OutputDirectory, sourceName))): os.system("rm %s/srcmap_*%s.*" % OutputDirectory)
		if (os.path.isfile("%s/modelmap_*%s.fits" % (OutputDirectory, sourceName))): os.system("rm %s/modelmap_*%s.fits" % OutputDirectory)
		if (os.path.isfile("%s/Model%s_Fit.xml" % (OutputDirectory, sourceName))): os.system("rm %s/Model%s_Fit.xml" % OutputDirectory)
		if (os.path.isfile("%s/Model%s_Fit_Modified.xml" % (OutputDirectory, sourceName))): os.system("rm %s/Model%s_Fit_Modified.xml" % OutputDirectory)


	# Move things back from the scratch directory if this is a batch job
	if JobID != None:

		print "\nMoving results from scratch disk..."
		if os.path.isdir(OutputDirectory_NSF):
			print "mv %s/* %s" % (OutputDirectory, OutputDirectory_NSF)
			os.system("mv %s/* %s" % (OutputDirectory, OutputDirectory_NSF))
		else:
			print "mv %s %s" % (OutputDirectory, OutputDirectory_NSF)
			os.system("mv %s %s" % (OutputDirectory, OutputDirectory_NSF))

		# Delete the scratch space directory
		if(os.path.isdir(OutputDirectory)==True):
			cmd = "rm -R %s/*" % OutputDirectory
			print cmd
			os.system(cmd)	


	print '\n\nAnalysis Complete.\n'
	return


##########################################################################################

def tsmap(sourceName, ra, dec, tmin, tmax, dra=7, ddec=7, binsize=0.15, emin=100, emax=1e5, tsMin=25, irfs='P8R2_SOURCE_V6', association='none', getData=True, generateFiles=True, performLikelihoodFit=False, makeSourceMap=False, cleanup=True, justGetData=False, statistic='UNBINNED', skipDiffuseResponse=False):

	SourceAnalysis(sourceName, ra, dec, tmin, tmax, dra=dra, ddec=ddec, binsize=binsize, emin=emin, emax=emax, tsMin=tsMin, irfs=irfs, association=association, getData=getData, generateFiles=generateFiles, performLikelihoodFit=False, makeSourceMap=False, maketsmap=True, makeSummaryMap=False, makeModelMap=False, cleanup=True, cleanupAll=False, nuke=False, justGetData=False, fixedModel=True, statistic='UNBINNED', skipDiffuseResponse=False, resultsPrefix='Results')


##########################################################################################

def dtsmap(sourceName, ra, dec, tmin, tmax, dra=7, ddec=7, binsize=0.15, emin=100, emax=1e5, action=None, tsMin=25, association='none', irfs='P8R2_SOURCE_V6', getData=True, generateFiles=True, performLikelihoodFit=True, nuke=True, fixedModel=True, statistic='UNBINNED', verbose=False, test=False, maxJobs=200, batch=True, interpolation='nearest', resubmit=True, plotMaps=True, pickleResults=False, region=None, removeWeakSources=True, ROI=None, maxValue=None, fixIndex=True):

	# Define the home directory
	LikelihoodDirectory = '/nfs/slac/g/ki/ki08/kocevski/Likelihood'

	# Define the script directory
	ScriptDirectory = LikelihoodDirectory + '/Scripts'

	# Define an output directory
	OutputDirectory = "%s/Results/%s" % (LikelihoodDirectory, sourceName)

	# Check if the user has write access to this location.  If not, create a subdirectory in the currnet working directory
	if os.access(OutputDirectory, os.W_OK) == False:
		OutputDirectory = os.getcwd() + '/%s' % sourceName


	# Define the log directory
	LogDirectory = OutputDirectory + '/dtsmap'
	JobsDirectory = OutputDirectory + '/dtsmap'

	# Generate the output file names
	tsMapFigure = OutputDirectory + '/dTSmap_%s.png' % sourceName
	ulMapFigure = OutputDirectory + '/dULmap_%s.png' % sourceName
	tsMapPickle = OutputDirectory + '/dTSmap_%s.pickle' % sourceName
	ulMapPickle = OutputDirectory + '/dULmap_%s.pickle' % sourceName
	upperLimitResultsFile = OutputDirectory + '/dULResults_%s.txt' % sourceName

	# Create the output directory if necessary
	if 'submit' in action:

		print "\nCreating custom output directory:\n%s" % OutputDirectory        
		cmd = "mkdir -p " + OutputDirectory
		os.system(cmd)

		# Create the output directory
		print "\nCreating log/job directory:\n%s" % LogDirectory        
		cmd = "mkdir -p " + LogDirectory
		os.system(cmd)

	# Convert everything to a float
	ra = float(ra)
	dec = float(dec)
	dra = float(dra)
	ddec = float(ddec)
	binsize = float(binsize)

	# Calculate the pixel scale
	xsize = int( dra / (binsize) )
	ysize = int( ddec / (binsize) )

	# Create the wcs string
	headerString = """
	NAXIS   =                    2
	NAXIS1  =                   %s
	NAXIS2  =                   %s
	CTYPE1  = 'RA---AIT'
	CRPIX1  =                   %s
	CRVAL1  =                   %s
	CDELT1  =                   %s
	CUNIT1  = 'deg     '
	CTYPE2  = 'DEC--AIT'          
	CRPIX2  =                   %s
	CRVAL2  =                   %s
	CDELT2  =                   %s
	CUNIT2  = 'deg     '
	CROTA2  =                    0
	""" % (xsize, ysize, (xsize+1)/2.0, ra, -1*binsize, (ysize+1)/2.0, dec, binsize)

	# Create a fits header using the wcs string
	target_header = fits.Header.fromstring(headerString, sep='\n')

	# Create the wcs object from the fits header
	wcs = WCS(target_header)

	# Loop through all combinations of Ra and Dec and submit a job for each one
	RaDecPairs = {}
	binNumbers = []
	binNumber = 0

	# Loop through each x, y pair
	for x in range(xsize):
		for y in range(ysize):

			# Get the ra and dec for given pixel
			raStep = wcs.all_pix2world(x,y,0)[0].item()
			decStep = wcs.all_pix2world(x,y,0)[1].item()

			# Record the coordinate pair for the bin number
			RaDecPairs[binNumber] = [raStep, decStep]

			# Increment the bin number
			binNumber = binNumber + 1

	# Calculate the corner coordinates
	raDec_lowerLeft = wcs.all_pix2world(0,0,0)
	raDec_upperLeft = wcs.all_pix2world(0,ysize-1,0)
	raDec_lowerRight = wcs.all_pix2world(xsize-1,0,0)
	raDec_upperRight = wcs.all_pix2world(xsize-1,ysize-1,0)

	# Display the coordinate grid information
	print "\nTS map grid size: %s x %s" % (xsize, ysize)
	print "Bin size: %s deg" % binsize
	print " "
	print "(%.2f, %.2f)   (%.2f, %.2f)" % (raDec_upperLeft[0].item(), raDec_upperLeft[1].item(), raDec_upperRight[0].item(), raDec_upperRight[1].item())
	print " -----------------------------"
	print " |                           |"
	print " |                           |"
	print " |                           |"
	print " |                           |" 
	if len(str("%.2f" % float(ra))) == 6:
		if len(str("%.2f" % float(dec))) == 4:
			print " |      (%.2f, %.2f)       |" % (float(ra), float(dec))
		else:
			print " |      (%.2f, %.2f)     |" % (float(ra), float(dec))
	elif len(str("%.2f" % float(ra))) == 5:
		if len(str("%.2f" % float(dec))) == 4:
			print " |      (%.2f, %.2f)       |" % (float(ra), float(dec))
		else:		
			print " |      (%.2f, %.2f)       |" % (float(ra), float(dec))
	else:
		if len(str("%.2f" % float(dec))) == 4:
			print " |      (%.2f, %.2f)       |" % (float(ra), float(dec))
		else:		
			print " |      (%.2f, %.2f)       |" % (float(ra), float(dec))			
	print " |             x             |"
	print " |                           |"
	print " |                           |"
	print " |                           |"
	print " |                           |"					
	print " -----------------------------"
	print "(%.2f, %.2f)    (%.2f, %.2f)" % (raDec_lowerLeft[0].item(), raDec_lowerLeft[1].item(), raDec_lowerRight[0].item(), raDec_lowerRight[1].item())

	# Warn the user that they haven't specified an action
	if action == None:
		print "\nError: No action specified.\n\nOptions include:\n'submit' - Submit jobs to the batch farm\n'collect' - Collect the results from submitted jobs"
		return

	# Submit the jobs
	if 'submit' in action or 'both' in action:

		# Print the total number of jobs
		print "\nSubmitting a total of %s Jobs\n" % (xsize * ysize)

		# Keep track of the jobs submitted. This allows the user to specify how many jobs should be running at any given time.
		JobsInQueue = 0

		# Get the bin numbers
		binNumbers = RaDecPairs.keys()

		for binNumber in binNumbers:

			# Get the ra and dec step
			raStep = RaDecPairs[binNumber][0]
			decStep = RaDecPairs[binNumber][1]

			# Check to see if the number of jobs in the queue is greater than the max.  If so, wait for the jobs to leave the queue before submitting more.
			while JobsInQueue >= int(maxJobs):

				print "\nMaximum number of submitted jobs (%s) reached.  Waiting..." % maxJobs
				print "Total number of remaining jobs to submit: %s" % remainingJobs

				# Wait 60 seconds before polling the job statuses
				if test == False:
					time.sleep(60)

				# Get the number of jobs actively in the queue
				command = "bjobs -g %s | wc" % JobsDirectory	
				process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
				lines = process.stdout.readlines()
				
				# Extract the number of jobs that are running
				if 'No' in lines[0].split()[0]:
					JobsInQueue = 0
				else:
					JobsInQueue = int(lines[0].split()[0])
					
			# Setup the job
			logfile = LogDirectory + "/dtsmap_bin%s.log" % binNumber

			# Make the job name
			jobName = "dtsmap_b%s" % binNumber				

			# Generate a unique source name
			uniqueSourceName = '%s_b%s' % (sourceName, binNumber)

			# Construct the command		
			command = """LikelihoodAnalysis.py %s %s %s %s %s irfs=%s fixedModel=%s nuke=%s removeWeakSources=%s fixIndex=%s""" % (uniqueSourceName, raStep, decStep, tmin, tmax, irfs, fixedModel, nuke, removeWeakSources, fixIndex)

			# Specify where to find the python script
			if ScriptDirectory != None:
				command = ScriptDirectory + "/" + command 

			# Construct the process call
			process = 'bsub -oo ' + logfile + ' -J ' + jobName + ' -W 2880 -R rhel60 -g ' + JobsDirectory + ' "' + command + '"'
		
			# Display the command
			if verbose == True:
				print process

			# Start the process
			if test == False:
				if batch == True:

					# Wait 2 seconds between job submissions
					time.sleep(2)

					# Execute the command
					subprocess.call(process, shell=True)

				else:

					# Execture the command
					os.system(command)


			# Increment the bin number
			JobsInQueue = JobsInQueue + 1

			# Get the number of remaining jobs
			remainingJobs = ((xsize * ysize) - binNumber)


		print "\nTotal number of jobs submitted: %s" % binNumber
		print "All jobs submitted."


		# # Check to see if the number of jobs in the queue is greater than the max.  If so, wait for the jobs to leave the queue before submitting more.
		# while JobsInQueue > 0:

		# 	print "\nTotal number of remaining jobs: %s" % remainingJobs

		# 	# Wait 60 seconds before polling the job statuses
		# 	if test == False:
		# 		time.sleep(60)

		# 	# Get the number of jobs actively in the queue
		# 	command = "bjobs -g %s | wc" % JobsDirectory	
		# 	process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
		# 	lines = process.stdout.readlines()
			
		# 	# Extract the number of jobs that are running
		# 	if 'No' in lines[0].split()[0]:
		# 		JobsInQueue = 0
		# 	else:
		# 		JobsInQueue = int(lines[0].split()[0])

		# if 'submit' in action:
		# 	return 


	if 'collect' in action or 'both' in action:

		print '\nCollecting results...'

		failedJobs = 0
		binNumbersRead = []
		TSs = []
		photonFluxUpperLimits = []

		# Create dictionaries to store the results
		LikelihoodResults = {}
		UpperLimitResults = {}

		numberOfReturnCode102 = 0
		numberOfReturnCode0 = 0
		numberOfReturnCode1 = 0
		numberOfReturnCode2 = 0

		print '\nReading logs from:\n%s\n' % LogDirectory

		# Get the bin numbers
		binNumbers = RaDecPairs.keys()

		# Calculate the total number
		totalNumberOfBins = xsize * ysize

		# Loop through each bin numbers
		for binNumber in binNumbers:

			# print binNumber, float(binNumber), int(float(binNumber)), totalNumberOfBins, int(float(binNumber))/totalNumberOfBins

			sys.stdout.write('\r')
			sys.stdout.write("Progress: %d%%" % ( ( float(binNumber) / float(totalNumberOfBins) ) * 100 ) )
			sys.stdout.flush()

			# Get the ra and dec step for each bin number
			raStep = RaDecPairs[binNumber][0]
			decStep = RaDecPairs[binNumber][1]

			# Generate a unique source name
			uniqueSourceName = '%s_b%s' % (sourceName, binNumber)

			# Define the log file
			logfile = "%s/dtsmap_bin%s.log" % (LogDirectory, binNumber)

			# Set default values
			ts = numpy.nan
			photonFluxUpperLimit = numpy.nan

			# parse the log file in search of the results
			if os.path.isfile(logfile):
				# print logfile

				# Loop through each line in the file
				for line in fileinput.input([logfile]):

					# Catch and resubmit failed likelihood fits
					if 'Return Code: 102' in line:

						failedJobs = failedJobs + 1
						numberOfReturnCode102 = numberOfReturnCode102 + 1

					if 'Return Code: 0' in line:
						numberOfReturnCode0 = numberOfReturnCode0 + 1

					if 'Return Code: 1' in line and 'Return Code: 102' not in line:
						numberOfReturnCode1 = numberOfReturnCode1 + 1

					if 'Return Code: 2' in line:
						numberOfReturnCode2 = numberOfReturnCode2 + 1

					# Extract the TS value
					if 'TS:' in line:
						LineContents = line.split()	
						ts = float(LineContents[1])

						# Make sure we don't have negative values
						if ts < 0: ts = 0

					# Extract the photon flux upper limit
					if 'Photon Flux Upper Limit:' in line:
						LineContents = line.split()	
						photonFluxUpperLimit = LineContents[4]

			else:
				
				print 'file not found: %s' % logfile
				failedJobs = failedJobs + 1

			# Store the ts and photon flux values
			# binNumbersRead.append(binNumber)
			# TSs.append(ts)
			# photonFluxUpperLimits.append(photonFluxUpperLimit)

			LikelihoodResults[binNumber] = ts
			UpperLimitResults[binNumber] = photonFluxUpperLimit

		print "\nDone."

		print "\nNumber of return code 0: %s" % numberOfReturnCode0
		print "Number of return code 1: %s" % numberOfReturnCode1
		print "Number of return code 2: %s" % numberOfReturnCode2
		print "Number of return code 102: %s\n" % numberOfReturnCode102

		# print len(binNumbersRead)
		# print len(TSs)

		# Put the results in a dictionary for easy retrieval later
		# LikelihoodResults = {key:value for key, value in zip(binNumbersRead,TSs)}
		# UpperLimitResults = {key:value for key, value in zip(binNumbersRead,photonFluxUpperLimits)}

		# Make a matrix to store the ts values.  Saving any missing values as NaNs
		TSMap = numpy.zeros(shape=(xsize, ysize))
		UpperLimitMap = numpy.zeros(shape=(xsize, ysize))

		binNumber = 0
		xStep = 0
		yStep = 0
		badBinFound = False

		for x in range(xsize):
			for y in range(ysize):
	
				if binNumber in LikelihoodResults:

					if numpy.isnan(LikelihoodResults[binNumber]):
						badBinFound = True

					else:

						badBinFound = False
						TSMap[x][y] = LikelihoodResults[binNumber]
						UpperLimitMap[x][y] = UpperLimitResults[binNumber]			
						pass

				else:

					badBinFound = True

				if badBinFound == True:

					print "Bin not found: %s" % binNumber
					print "%s/dtsmap_bin%s.log" % (LogDirectory, binNumber)
					TSMap[x][y] = numpy.nan
					UpperLimitMap[x][y] = numpy.nan		

					if resubmit == True:

						# Get the ra and dec step
						raStep = RaDecPairs[binNumber][0]
						decStep = RaDecPairs[binNumber][1]

						# Define the log file
						logfile = "%s/dtsmap_bin%s.log" % (LogDirectory, binNumber)

						# Construct the command		
						command = """LikelihoodAnalysis.py %s %s %s %s %s irfs=%s fixedModel=%s nuke=%s removeWeakSources=%s fixIndex=%s""" % (uniqueSourceName, raStep, decStep, tmin, tmax, irfs, fixedModel, nuke, removeWeakSources, fixIndex)

						# Specify where to find the python script
						if ScriptDirectory != None:
							command = ScriptDirectory + "/" + command 

						# Make the job name
						jobName = "dtsmap_b%s" % binNumber	

						# Construct the process call
						process = 'bsub -oo ' + logfile + ' -J ' + jobName + ' -W 2880 -R rhel60 -g ' + JobsDirectory + ' "' + command + '"'
					
						# Display the command
						if verbose == True:
							print process

						# Start the process
						if test == False:
							if batch == True:
								print process
								subprocess.call(process, shell=True)
							else:
								print command
								os.system(command)


				# if TSMap[x][y] < -1:
				# 	print "Bin with TS < -1: %s" % binNumber
				# 	print "%s/dtsmap_bin%s.log" % (LogDirectory, binNumber)					
				# 	TSMap[x][y] = numpy.nan
				# 	UpperLimitMap[x][y] = numpy.nan										
				# 	pass

				binNumber = binNumber + 1
				pass

			pass


		# Finding the maximum TS
		MaxTS = numpy.nanmax(TSMap)
		MaxBin = numpy.nanargmax(TSMap)


		# Check to see if a maximum couldn't be found.
		if numpy.isnan(MaxBin) == True:
			print "\nAnalysis Complete."
			print "Maximum TS: %s" % 'None'
			print "Coordinates: RA = %s, Dec = %s" % ('NA', 'NA')
			print "*** Unable to locate maximum TS ***"

		else:

			MaxRa = RaDecPairs[MaxBin][0]
			MaxDec = RaDecPairs[MaxBin][1]


		# Rotate and flip the TS map matrix in order to it to match ra and dec ordering conventions
		TSMap = numpy.rot90(TSMap)
		TSMap = numpy.flipud(TSMap)
		# TSMap = numpy.fliplr(TSMap)

		# Rotate and flip the upper limit matrix in order to it to match ra and dec ordering conventions
		UpperLimitMap = numpy.rot90(UpperLimitMap)
		UpperLimitMap = numpy.flipud(UpperLimitMap)
		# UpperLimitMap = numpy.fliplr(UpperLimitMap)


		# Calculate the media upper limit within the ROI
		if ROI != None:

			# Save the upper limits results
			output = open(upperLimitResultsFile, 'w')
			output.write("%s Median Upper Limit Results:\n" % sourceName)

			print "Median Upper Limits:"

			# Convert float values to a list
			if type(ROI) is not list:
				ROI = [ROI]

			# Loop throuh each ROI size and calculate the median upper limit
			for radius in ROI:

				# Start a list to contain all of the upper limits within the specified ROI
				upperlimits = []

				# Loop through each ra dec pair and calculate the angular seperation from the center of the dtsmap
				for binNumber in RaDecPairs.keys():

					# Get the coordinates
					ra_bin = RaDecPairs[binNumber][0]
					dec_bin = RaDecPairs[binNumber][1]

					# Get the angular seperation
					distance = AngularSeperation(ra, dec, ra_bin, dec_bin)

					# Record the upper limit if it's within the ROI
					if distance <= radius:
						upperlimits.append( UpperLimitResults[binNumber] )

				# Turn the upper limits list into a numpy array
				upperlimits = numpy.array(upperlimits)
				upperlimits = upperlimits.astype(float)

				# Print the results
				print 'ROI < %s deg: %s photons cm-2 s-1' % (radius, numpy.median(upperlimits))
				output.write('ROI < %s deg: %s photons cm-2 s-1\n' % (radius, numpy.median(upperlimits)) )

			# Close the file
			output.close()
			print '\nSaving median upper limits results to: %s\n' % upperLimitResultsFile

		# Plot the results
		if plotMaps == True:

			# Plot the ts map
			print '\nCreating ts map image...'
			plotImage(TSMap, wcs, filename=tsMapFigure, region=region, colorbarLabel='TS', ROI=ROI, stretchColorValue=25, maxValue=maxValue)

			# Plot the upper limit map
			print 'Creating upper limits map image...'
			# plotImage(UpperLimitMap, wcs, filename=ulMapFigure, region=None, colorbarLabel=r'log Photon Flux UL (95%) (ph cm$^{-2}$ s$^{-1}$)', ROI=ROI, stretchColorValue=3e-7, maxValue=maxValue)
			plotImage(UpperLimitMap, wcs, filename=ulMapFigure, region=None, colorbarLabel=r'log Photon Flux UL (95%) (ph cm$^{-2}$ s$^{-1}$)', ROI=ROI, maxValue=maxValue)

		# Pickle the results
		if pickleResults == True:
			print "\nPickling results..."
			print 'Saving ts map results to: %s' % tsMapPickle
			pickle.dump( TSMap, open( tsMapPickle, "wb" ) )
			print 'Saving upper limit map results to: %s\n' % ulMapPickle		
			pickle.dump( UpperLimitMap, open( ulMapPickle, "wb" ) )


		return TSMap, UpperLimitMap


##########################################################################################

def Lightcurve(sourceName, ra, dec, tmin, tmax, dt, emin=100, emax=1e5, tsMin=25, irfs='P8R2_SOURCE_V6', association='none', action='both', getData=True, generateFiles=True, nuke=True, fixedModel=True, statistic='UNBINNED', verbose=True, test=False, maxJobs=200, batch=True, resubmit=True, plotLC=True, pickleResults=False, removeWeakSources=True, logCenter=True, annotations=None, xlog=False, ylog=True, fixIndex=True):

	# Define the home directory
	LikelihoodDirectory = '/nfs/slac/g/ki/ki08/kocevski/Likelihood'

	# Define the script directory
	ScriptDirectory = LikelihoodDirectory + '/Scripts'

	# Define an output directory
	OutputDirectory = "%s/Results/%s" % (LikelihoodDirectory, sourceName)

	# Check if the user has write access to this location.  If not, create a subdirectory in the currnet working directory
	if os.access(OutputDirectory, os.W_OK) == False:
		OutputDirectory = os.getcwd() + '/%s' % sourceName

	# Define the log directory
	LogDirectory = OutputDirectory + '/lightcurve'
	JobsDirectory = OutputDirectory + '/lightcurve'

	# Create the output directory if necessary
	if 'submit' in action:

		print "\nCreating custom output directory:\n%s" % OutputDirectory        
		cmd = "mkdir -p " + OutputDirectory
		os.system(cmd)

		# Create the output directory
		print "\nCreating log/job directory:\n%s" % LogDirectory        
		cmd = "mkdir -p " + LogDirectory
		os.system(cmd)


	# Generate the output plot file names
	tsPlot = OutputDirectory + '/tsPlot.png' 
	photonFluxPlot = OutputDirectory + '/PhotonFluxPlot.png'
	photonIndexPlot = OutputDirectory + '/PhotonIndexPlot.png'

	# Generate the pickle filenames
	lightcurvePickle = OutputDirectory + '/lightcuve.pickle'

	# Generate the temporal bins
	timebins = numpy.arange(tmin, tmax, dt)

	# Calculate the duration 
	duration = tmax - tmin

	# Create the tstart and tstop arrays
	tstarts = timebins[0:-1]
	tstops = timebins[1:]

	# Loop through the timebins and submit a job for each one
	timebinPairs = {}
	binNumbers = []
	binNumber = 0

	for tstart, tstop in zip(tstarts, tstops):

		# Record the coordinate pair for the given bin number for later use
		timebinPairs[binNumber] = [tstart, tstop]

		# Increment the bin number
		binNumber = binNumber + 1

	# Submit the jobs
	if 'submit' in action or 'all' in action:

		# Print the total number of jobs
		print "\nSubmitting a total of %s Jobs\n" % len(tstarts)

		# Keep track of the jobs submitted. This allows the user to specify how many jobs should be running at any given time.
		JobsInQueue = 0

		# Get the bin numbers
		binNumbers = timebinPairs.keys()

		for binNumber in binNumbers:

			# Get the ra and dec step
			tstart = timebinPairs[binNumber][0]
			tstop = timebinPairs[binNumber][1]

			# Check to see if the number of jobs in the queue is greater than the max.  If so, wait for the jobs to leave the queue before submitting more.
			while JobsInQueue >= int(maxJobs):

				print "\nMaximum number of submitted jobs (%s) reached.  Waiting..." % maxJobs
				print "Total number of remaining jobs to submit: %s" % remainingJobs

				# Wait 60 seconds before polling the job statuses
				if test == False:
					time.sleep(60)

				# Get the number of jobs actively in the queue
				command = "bjobs -g %s | wc" % JobsDirectory	
				process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
				lines = process.stdout.readlines()
				
				# Extract the number of jobs that are running
				if 'No' in lines[0].split()[0]:
					JobsInQueue = 0
				else:
					JobsInQueue = int(lines[0].split()[0])
					
			# Setup the job
			logfile = LogDirectory + "/lc_bin%s.log" % binNumber

			# Make the job name
			jobName = "lc_b%s" % binNumber				

			# Generate a unique source name
			uniqueSourceName = '%s_b%s' % (sourceName, binNumber)

			# Construct the command		
			command = """LikelihoodAnalysis.py %s %s %s %s %s irfs=%s fixedModel=%s nuke=%s removeWeakSources=%s fixIndex=%s""" % (uniqueSourceName, ra, dec, tstart, tstop, irfs, fixedModel, nuke, removeWeakSources, fixIndex)

			# Specify where to find the python script
			if ScriptDirectory != None:
				command = ScriptDirectory + "/" + command 

			# Construct the process call
			process = 'bsub -oo ' + logfile + ' -J ' + jobName + ' -W 2880 -R rhel60 -g ' + JobsDirectory + ' "' + command + '"'
		
			# Display the command
			if verbose == True or test == True:
				print process

			# Start the process
			if test == False:
				if batch == True:

					# Wait 2 seconds between job submissions
					time.sleep(2)

					# Execute the command
					subprocess.call(process, shell=True)

				else:

					# Execture the command
					os.system(command)


			# Increment the bin number
			JobsInQueue = JobsInQueue + 1

			# Get the number of remaining jobs
			remainingJobs = (len(tstarts) - binNumber)


		print "\nTotal number of jobs submitted: %s" % binNumber
		print "All jobs submitted."

		return 		



	if 'collect' in action or 'all' in action:

		print '\nCollecting results...'

		failedJobs = 0
		binNumbersRead = []
		TSs = []
		photonFlux = []	
		photonFluxError = []					
		photonFluxUpperLimits = []

		# Create dictionaries to store the results
		timebinResults = []
		tsResults = []
		photonFluxUpperLimitResults = []
		photonFluxResults = []
		photonFluxErrorResults = []
		photonIndexResults = []
		photonIndexErrorResults = []

		numberOfReturnCode102 = 0
		numberOfReturnCode0 = 0
		numberOfReturnCode1 = 0
		numberOfReturnCode2 = 0

		print '\nReading logs from:\n%s\n' % LogDirectory

		# Get the bin numbers
		binNumbers = timebinPairs.keys()

		# Calculate the total number
		totalNumberOfBins = len(tstarts)

		# Loop through each bin numbers
		for binNumber in binNumbers:

			# Display the progress
			sys.stdout.write('\r')
			sys.stdout.write("Progress: %d%%" % ( ( float(binNumber) / float(totalNumberOfBins) ) * 100 ) )
			sys.stdout.flush()

			# Get the ra and dec step for each bin number
			tstart = timebinPairs[binNumber][0]
			tstop = timebinPairs[binNumber][1]

			# Generate a unique source name
			uniqueSourceName = '%s_b%s' % (sourceName, binNumber)

			# Define the log file
			logfile = "%s/lc_bin%s.log" % (LogDirectory, binNumber)

			# Set default values
			ts = numpy.nan
			photonFlux = numpy.nan
			photonFluxError = numpy.nan
			photonFluxUpperLimit = numpy.nan
			photonIndex = numpy.nan
			photonIndexError = numpy.nan

			# parse the log file in search of the results
			if os.path.isfile(logfile):
				# print logfile

				# Loop through each line in the file
				for line in fileinput.input([logfile]):

					# Catch and resubmit failed likelihood fits
					if 'Return Code: 102' in line:

						failedJobs = failedJobs + 1
						numberOfReturnCode102 = numberOfReturnCode102 + 1

					if 'Return Code: 0' in line:
						numberOfReturnCode0 = numberOfReturnCode0 + 1

					if 'Return Code: 1' in line and 'Return Code: 102' not in line:
						numberOfReturnCode1 = numberOfReturnCode1 + 1

					if 'Return Code: 2' in line:
						numberOfReturnCode2 = numberOfReturnCode2 + 1

					# Extract the TS value
					if 'TS:' in line:
						LineContents = line.split()	
						ts = float(LineContents[1])

						# Make sure we don't have negative values
						if ts < 0: ts = 0

					# Extract the photon flux upper limit
					if 'Photon Flux Upper Limit:' in line:
						LineContents = line.split()	
						photonFluxUpperLimit = float(LineContents[4])

					# Extract the photon flux upper limit
					if 'Photon Flux:' in line:
						LineContents = line.split()	
						photonFlux = float(LineContents[2])
						photonFluxError = float(LineContents[4])

					# Extract the photon flux upper limit
					if 'Photon Index:' in line:
						LineContents = line.split()	
						photonIndex = float(LineContents[2])
						photonIndexError = float(LineContents[4])


			else:
				print 'file not found: %s' % logfile
				failedJobs = failedJobs + 1

			# Store the results
			tsResults.append(ts)
			photonFluxUpperLimitResults.append(photonFluxUpperLimit)

			photonFluxResults.append(photonFlux)
			photonFluxErrorResults.append(photonFluxError)

			photonIndexResults.append(photonIndex)
			photonIndexErrorResults.append(photonIndexError)

			# Store the bin centers
			if logCenter == True:
				timebin_center = (tstart*tstop)**0.5 
			else:
				timebin_center = tstart + (tstart+tstop)*0.5

			timebinResults.append(timebin_center)

			# Determine which timebins have detections 
			detections = numpy.where(tsResults >= tsMin)[0]
			nondetections = numpy.where(tsResults >= tsMin)[0]


		# Make sure everything is an array
		timebinResults = numpy.array(timebinResults)
		tsResults = numpy.array(tsResults)
		photonFluxResults = numpy.array(photonFluxResults)
		photonFluxErrorResults = numpy.array(photonFluxErrorResults)
		photonIndexResults = numpy.array(photonIndexResults)
		photonIndexErrorResults = numpy.array(photonIndexErrorResults)
		photonFluxUpperLimitResults = numpy.array(photonFluxUpperLimitResults)

		print "\nDone."

		print "\nNumber of return code 0: %s" % numberOfReturnCode0
		print "Number of return code 1: %s" % numberOfReturnCode1
		print "Number of return code 2: %s" % numberOfReturnCode2
		print "Number of return code 102: %s\n" % numberOfReturnCode102

		print "\nLightcurve Results:"
		print 'Time, TS, PhotonFlux +/-, PhotonIndex +/-, UpperLimit'
		for tstart, tstop, ts, photonFlux, photonFluxError, photonIndex, photonIndexError, photonFluxUpperLimit in zip(tstarts, tstops, tsResults, photonFluxResults, photonFluxErrorResults, photonIndexResults, photonIndexErrorResults, photonFluxUpperLimitResults):
			print tstart, tstop, ts, photonFlux, photonFluxError, photonIndex, photonIndexError, photonFluxUpperLimit




		# Plot the results
		if plotLC == True:

			# Make the plot look pretty
			plot.rc('font',**{'size':12,'weight':'medium'})
			plot.rcParams['xtick.major.width'] = 1.0
			plot.rcParams['xtick.minor.width'] = 1.0
			plot.rcParams['ytick.major.width'] = 1.0
			plot.rcParams['ytick.minor.width'] = 1.0
			plot.rc('axes', linewidth=1.0)
			# plot.rcParams['figure.figsize'] = 10, 6.39
			plot.rcParams['legend.fontsize'] = 12
			plot.rc('font', family='sans-serif') 
			plot.rc('font', serif='Helvetica') 	
			# plot.rcParams['figure.figsize'] = 10, 6.39
			plot.rcParams['figure.figsize'] = 9, 5.75

			# Create the photon flux figure
			fig, ax = plot.subplots()

			# Setup the y-axis
			# majorLocator   = MultipleLocator(2)
			# majorFormatter = FormatStrFormatter('%d')
			# minorLocator   = MultipleLocator(0.5)
			# ax.yaxis.set_major_locator(majorLocator)
			# ax.yaxis.set_major_formatter(majorFormatter)
			# ax.yaxis.set_minor_locator(minorLocator)

			# # Setup the x-axis
			# majorLocator   = MultipleLocator(int(duration)/10)
			# majorFormatter = FormatStrFormatter('%d')
			# minorLocator   = MultipleLocator(duration/40)
			# ax.xaxis.set_major_locator(majorLocator)
			# ax.xaxis.set_major_formatter(majorFormatter)
			# ax.xaxis.set_minor_locator(minorLocator)

			# Setup the labels
			plot.xlabel('Time (sec)')
			plot.ylabel(r'Photon Flux cm$^{-2}$ s$^{-1}$')

			if xlog == True:
				plot.xscale('log')

			if ylog == True:
				plot.yscale('log')

			# Plot the flux values
			plot.errorbar(timebinResults[detections], photonFluxResults[detections], yerr=photonFluxErrorResults[detections], fmt='o', color='#3e4d8b', ecolor='#3e4d8b', markeredgecolor='black', label='Photon Flux')

			# Plot the upper limits
			plot.scatter(timebinResults[nondetections], photonFluxUpperLimitResults[nondetections], marker='v', color='#3e4d8b', edgecolor='black', label='Photon Flux')

			# Save the plot
			print '\nSaving photon flux plot to:\n%s' % photonFluxPlot
			plot.savefig(photonFluxPlot, dpi=72)


			# Create the photon flux figure
			fig, ax = plot.subplots()

			# Setup the y-axis
			# majorLocator   = MultipleLocator(2)
			# majorFormatter = FormatStrFormatter('%d')
			# minorLocator   = MultipleLocator(0.5)
			# ax.yaxis.set_major_locator(majorLocator)
			# ax.yaxis.set_major_formatter(majorFormatter)
			# ax.yaxis.set_minor_locator(minorLocator)

			# # Setup the x-axis
			# majorLocator   = MultipleLocator(int(duration)/10)
			# majorFormatter = FormatStrFormatter('%d')
			# minorLocator   = MultipleLocator(duration/40)
			# ax.xaxis.set_major_locator(majorLocator)
			# ax.xaxis.set_major_formatter(majorFormatter)
			# ax.xaxis.set_minor_locator(minorLocator)

			# Setup the labels
			plot.xlabel('Time (sec)')
			plot.ylabel('TS')

			# Plot the ts values
			plot.scatter(timebinResults, tsResults, marker='o', color='#3e4d8b', edgecolor='black', label='TS')

			# Save the plot
			print '\nSaving photon TS plot to:\n%s' % tsPlot
			plot.savefig(tsPlot, dpi=72)


			# Create the photon index figure
			fig, ax = plot.subplots()

			# Setup the y-axis
			# majorLocator   = MultipleLocator(2)
			# majorFormatter = FormatStrFormatter('%d')
			# minorLocator   = MultipleLocator(0.5)
			# ax.yaxis.set_major_locator(majorLocator)
			# ax.yaxis.set_major_formatter(majorFormatter)
			# ax.yaxis.set_minor_locator(minorLocator)

			# # Setup the x-axis
			# majorLocator   = MultipleLocator(int(duration)/10)
			# majorFormatter = FormatStrFormatter('%d')
			# minorLocator   = MultipleLocator(duration/40)
			# ax.xaxis.set_major_locator(majorLocator)
			# ax.xaxis.set_major_formatter(majorFormatter)
			# ax.xaxis.set_minor_locator(minorLocator)

			# Setup the labels
			plot.xlabel('Time (sec)')
			plot.ylabel('Photon Index')

			# Plot the flux values
			plot.errorbar(timebinResults[detections], photonIndexResults[detections], yerr=photonIndexErrorResults[detections], fmt='o', color='#3e4d8b', ecolor='#3e4d8b', markeredgecolor='black', label='Photon Index')

			# Save the plot
			print '\nSaving photon index plot to:\n%s' % photonIndexPlot
			plot.savefig(photonIndexPlot, dpi=72)


		# Pickle the results
		if pickleResults == True:
			print "\nPickling results..."
			print 'Saving light curve results to: %s' % lightcurvePickle
			pickle.dump( [tstarts, tstops, tsResults, photonFluxResults, photonFluxErrorResults, photonIndexResults, photonIndexErrorResults, photonFluxUpperLimitResults], open( lightcurvePickle, "wb" ) )




##########################################################################################


if __name__ == '__main__':


	if len(sys.argv) > 1:

		# Get the required keywords
		sourceName = sys.argv[1]
		ra = sys.argv[2]
		dec = sys.argv[3]
		tmin = sys.argv[4]
		tmax = sys.argv[5]

		# Extact the keywords
		kwargs = {}
		for keyword in sys.argv:
			if '=' in keyword:
				key, value = keyword.split('=', 1)
				kwargs[key] = value

		# Set the default values
		getData = True
		generateFiles = True
		justGetData = False			
		performLikelihoodFit = True
		makeSummaryMap = False
		makeModelMap = False			
		makeLightCurve = False
		makeApertureLightCurve = False
		skipDiffuseResponse = False
		statistic = 'UNBINNED'

		# emin and emax
		if 'emin' in kwargs:
			emin = kwargs['emin']
		else:
			emin = 100

		if 'emax' in kwargs:
			emax = kwargs['emax']
		else:
			emax = 1e5
	
		if 'dra' in kwargs:
			dra = kwargs['dra']
		else:
			dra = 7

		if 'ddec' in kwargs:
			ddec = kwargs['ddec']
		else:
			ddec = 7

		if 'binsize' in kwargs:
			binsize = kwargs['binsize']
		else:
			binsize = 0.15

		if 'association' in kwargs:
			association = kwargs['association']
		else:
			association = 'none'

		if 'getData' in kwargs:
			getData = BOOL(kwargs['getData'])
		else:
			getData = True

		if 'generateFiles' in kwargs:
			generateFiles = BOOL(kwargs['generateFiles'])
		else:
			generateFiles = True

		if 'justGetData' in kwargs:
			justGetData = BOOL(kwargs['justGetData'])
		else:
			justGetData = False

		if 'statistic' in kwargs:
			statistic = kwargs['statistic']
		else:
			statistic = 'UNBINNED'

		if 'fixedModel' in kwargs:
			fixedModel = BOOL(kwargs['fixedModel'])
		else:
			fixedModel = True

		if 'performLikelihoodFit' in kwargs:
			performLikelihoodFit = BOOL(kwargs['performLikelihoodFit'])
		else:
			performLikelihoodFit = True

		if 'maketsmap' in kwargs:
			maketsmap = BOOL(kwargs['maketsmap'])
		else:
			maketsmap = False

		if 'makedtsmap' in kwargs:
			makedtsmap = BOOL(kwargs['makedtsmap'])
		else:
			makedtsmap = False

		if 'makeModelMap' in kwargs:
			makeModelMap = BOOL(kwargs['makeModelMap'])
		else:
			makeModelMap = False

		if 'makeSummaryMap' in kwargs:
			makeSummaryMap = BOOL(kwargs['makeSummaryMap'])
		else:
			makeSummaryMap = False

		if 'cleanup' in kwargs:
			cleanup = BOOL(kwargs['cleanup'])
		else:
			cleanup = False

		if 'cleanupAll' in kwargs:
			cleanupAll = BOOL(kwargs['cleanupAll'])
		else:
			cleanupAll = False

		if 'nuke' in kwargs:
			nuke = BOOL(kwargs['nuke'])
		else:
			nuke = False

		if 'skipDiffuseResponse' in kwargs:
			skipDiffuseResponse = BOOL(kwargs['skipDiffuseResponse'])
		else:
			skipDiffuseResponse = False

		if 'resultsPrefix' in kwargs:
			resultsPrefix = kwargs['resultsPrefix']
		else:
			resultsPrefix = 'Results'

		if 'removeWeakSources' in kwargs:
			removeWeakSources  = BOOL(kwargs['removeWeakSources'])
		else:
			removeWeakSources  = False

		if 'maxJobs' in kwargs:
			maxJobs = kwargs['maxJobs']
		else:
			maxJobs = 1000

		if 'irfs' in kwargs:
			irfs = kwargs['irfs']
		else:
			irfs = 'P8R2_SOURCE_V6'

		if 'fixIndex' in kwargs:
			fixIndex  = BOOL(kwargs['fixIndex'])
		else:
			fixIndex  = True



		# Set some defaults if this is a tsmap request
		if maketsmap == True:
			performLikelihoodFit = False
			# fixedModel= True

		if makedtsmap == True:

			# Create a distributed ts map
			dtsmap(sourceName, ra, dec, tmin, tmax, dra=dra, ddec=ddec, binsize=binsize, action='submit', irfs=irfs, verbose=True, fixedModel=True, removeWeakSources=True, maxJobs=maxJobs, test=False, fixIndex=fixIndex)

		else:

			# Run a point source likelihood analysis
			SourceAnalysis(sourceName, ra, dec, tmin, tmax, emin=emin, emax=emax, dra=dra, ddec=ddec, tsmapBinSize=binsize, association=association, irfs=irfs, getData=getData, generateFiles=generateFiles, justGetData=justGetData, fixedModel=fixedModel, statistic=statistic, performLikelihoodFit=performLikelihoodFit, maketsmap=maketsmap, makeSummaryMap=makeSummaryMap, makeModelMap=makeModelMap, cleanup=cleanup, cleanupAll=cleanupAll, skipDiffuseResponse=skipDiffuseResponse, resultsPrefix=resultsPrefix, nuke=nuke, removeWeakSources =removeWeakSources, fixIndex=fixIndex)


	else:	

		print "usage: LikelihoodAnalysis, sourceName, ra, dec, tmin, tmax, emin=emin, emax=emax, association=association, getData=getData, generateFiles=generateFiles, justGetData=justGetData, fixedModel=fixedModel, statistic=statistic, performLikelihoodFit=performLikelihoodFit, maketsmap=tsmap, makeSummaryMap=makeSummaryMap, makeModelMap=makeModelMap, makeLightCurve=makeLightCurve, makeApertureLightCurve=makeApertureLightCurve, cleanup=cleanup, skipDiffuseResponse=skipDiffuseResponse)"
	 	sys.exit()


		
