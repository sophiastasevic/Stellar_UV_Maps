# Stellar_UV_Maps
---------------------------------------------------------------------------------------------
ALL FLUX IS SAVED IN UNITS OF LOG(Lo/pc^2) BUT PLOT IN LOG(Go) --> 1 Go = 3.98 Lo/pc^2 IF CONVERTING FLUX TO SAVE IN Go
OB_plotting: 

Data needed for calculating FUV/EUV/MS lifetime:

    -'ra' 
    -'dec'
    -'mass' (M_Sun)
    -'spt' (spectral type)

Text file output from SpType_to_Mass.py will need to be edited so that the first row is: 
'#ra   dec   teff   mass   spt'

If using SaveOutput(), need to edit output so that the first row is:
'#ra   dec   FUV   EUV   MS_lifetime'

FUV and EUV in L_Sun
MS lifetime in Myr
    
---------------------------------------------------------------------------------------------
SpType_to_Mass:

Data needed for mass estimates:

    -Spectroscopic information in at least one band AND:
    -Spectral Type OR Effective Temperature OR Spectroscopic information in a different band
    -Distance to cluster OR individual parallax informaton
    -At least a median cluster extinction
    -MIST isochrone for the age of the cluster with relevant spectroscopic information
  
Outputs text file with the following columns: RA, Dec, Effective Temp (log), Mass (M_Sun), Spectral Type

Use SpType_logT.txt for clusters closer to 3.5 Myrs old, and SpType_logT.txt for clusters
closer to 13 Myrs old

---------------------------------------------------------------------------------------------
SpType_logT:

Uses SpT to Teff conversion from Wright et al. 2015 for B5 and earlier, rest from Currie et al. 2010 
