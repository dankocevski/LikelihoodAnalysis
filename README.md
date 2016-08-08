# LikelihoodAnalysis

This module contains a collection of routines to perform likelihood analysis of gamma-ray data collected by the Large Area Telescope (LAT) on NASA's Fermi Gamma-ray Space Telescope.

An overview of the detection, flux determination and spectral modeling of Fermi LAT sources through the use of the maximum likelihood optimization technique is outlined in further detail here:

http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/likelihood_tutorial.html

### H3 Enviroment Setup

source /nfs/slac/g/ki/ki08/kocevski/Likelihood/Scripts/SetupScienceTools_10-01-01.csh

### H3 Usage Examples

Determine the detection significance of a new point source and calculate it's flux (or flux upper limit):
LikelihoodAnalysis.SourceAnalysis('CandidateSource', RA, Dec, METStart, METStop)

Search for a new point source over a 10 deg x 10 deg grid with 0.1 deg binning:
LikelihoodAnalysis.tsmap('CandidateSource', RA, Dec, METStart, METStop, dra=10, ddec=10, binsize=0.15)

Create a light curve for a new point source from METStart to METStop with a binning of 1 day (86400 sec):
LikelihoodAnalysis.Lightcurve('CandidateSource', RA, Dec, METStart, METStop, 86400) 
