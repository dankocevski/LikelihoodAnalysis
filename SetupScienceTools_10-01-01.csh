echo ""
echo "Setting Up Fermi Science Tools Enviroment..."
echo ""
#source /afs/slac.stanford.edu/g/ki/software/etc/csh.login.xray
#setenv PFILES "`/u/gl/jchiang/bin/strip_pfiles.py`"
setenv IDL_STARTUP /u/gl/jchiang/.idl_startup

# Define the science tools build to be used
set ST_VER = 10-01-01 
set VARIANT = Optimized
set BUILD_VERSION = redhat6-x86_64-64bit-gcc44

# Setup science tools
setenv GLAST_EXT /afs/slac/g/glast/ground/GLAST_EXT/${BUILD_VERSION}
setenv BINDIR ${BUILD_VERSION}-${VARIANT}
setenv INST_DIR /nfs/farm/g/glast/u35/ReleaseManagerBuild/${BUILD_VERSION}/${VARIANT}/ScienceTools/${ST_VER}
source ${INST_DIR}/bin/${BINDIR}/_setup.csh
setenv PYTHONPATH `pwd`/scripts:/u2/jchiang/local/lib/python2.7/site-packages:${PYTHONPATH}

# Set the pfiles
#setenv PFILES "`/u/gl/jchiang/bin/strip_pfiles.py`"
setenv PFILES /nfs/slac/g/ki/ki08/kocevski/Likelihood/pfiles

# Setup CALDB
setenv CALDB ${INST_DIR}/irfs/caldb/CALDB
setenv CALDBCONFIG ${CALDB}/software/tools/caldb.config
setenv CALDBALIAS ${CALDB}/software/tools/alias_config.fits

# Setup the diffuse models
setenv GALPROP_MODEL /nfs/slac/g/ki/ki08/kocevski/FAVA/diffuseModels/template_4years_P8_V2_scaled_trim.fits
setenv ISOTROPIC_MODEL /nfs/slac/g/ki/ki08/kocevski/FAVA/diffuseModes/iso_P8R2_SOURCE_V6_v06.txt

# Setup custom IRFs
setenv CUSTOM_IRF_DIR /afs/slac/g/glast/groups/canda/standard/custom_irfs
setenv CUSTOM_IRF_NAMES P8R2_TRANSIENT_SFR_V6,P8R2_TRANSIENT_ALLTKR_R100_V6,P8R2_TRANSIENT_TKRONLY_R100_V6,P8R2_TRANSIENT_TKRONLY_R020_V6,P8R2_TRANSIENT_TKRONLY_R010_V6,P8R2_TRANSIENT_R010_V6,P8R2_TRANSIENT_R020_V6,P8R2_TRANSIENT_R100_V6,P8R2_SOURCE_V6,P8R2_CLEAN_V6,P8R2_ULTRACLEAN_V6

# Add custom libraries to the library path
setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:/nfs/slac/g/ki/ki08/kocevski/LATBA/lib/geos-3.3.3_rhel6-64/lib

# Add custom python libraries to the python path
setenv PYTHONPATH {$PYTHONPATH}:/nfs/slac/g/ki/ki08/kocevski/Likelihood/Scripts

echo "ST_VER = $ST_VER"
echo "VARIANT = $VARIANT"
echo "BUILD_VERSION = $BUILD_VERSION"
echo "PFILES = $PFILES"
echo "GALPROP_MODEL = $GALPROP_MODEL"
echo "ISOTROPIC_MODEL = $ISOTROPIC_MODEL"
echo ""
echo "Done."
echo ""
