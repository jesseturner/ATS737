Here’s a description:
 
orbit_parms.exe reads TLE files, or TBUS files, or takes manual input and translates the data for one satellite into classic orbital parameters and additional parameters. You’ve already compiled this for Linux. This version works on DOS/Windows machines.

epoch_elem.exe is like orbit_parms.exe, except it processes all the satellites in a TLE file and outputs a tab-delimited table that can be imported into Excel for further analysis.

synoptic_elem.exe is like epoch_elem.exe except that it translates all of the orbital elements and other parameters to a single time that you choose. The output is also in a tab-delimited table. This makes comparison of the orbits easier.

Get_elem.exe, ephemeris.exe, and subtrack.exe are used in series.

get_elem.exe, like orbit_parms.exe, calculates orbital elements and other parameters necessary to create an ephemeris and writes them to a file (ORB_ELEM.DAT) in NAMELIST format so that other programs (like ephemeris.exe) can use them without repeatedly entering the TLE file name

ephemeris.exe reads ORB_ELEM.DAT and calculates an ephemeris for a time period that you specify. The output (EPHEMERIS.DAT) is an ephemeris file to be used by subtrack.exe or other programs.

subtrack.exe plots the ephemeris on a world map. It produces a file named subtrack.tif, which is in Tagged Image File Format (TIFF). Among the executables is tiff_to_gif.exe, which translates the TIFF to a GIF.(Unfortunately I can’t supply F90 code for tiff_to_gif.) Note that subtrack.exe requires that the file Continental_Outlines.DAT be in the same working directory as subtrack.exe.

