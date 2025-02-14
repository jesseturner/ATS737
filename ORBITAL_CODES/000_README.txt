README

get_elem.f90: gets orbital elements for a satellite from TLE files, TBUS bulletins (if they exist anymore) or manually

ephemeris.f90: uses the orbital elements produced by get_elem to calculate an ephemeris file for the satellite

NOTE:  get_elem, ephemeris, and subtrack are meant to be used in series to (1) get the necessary orbital elements (and associated parameters) from the TLE file and place them in a file for use by (2) ephemeris, which calculates an ephemeris table (position of the satellite over a selected period of time), and (3) subtrack, which plots the ephemeris on a global map.  

NOTE:  Routines subtrack (plots ephemeris on a global map) and space_view (plots ephemeris in the right ascension-declination coordinate system) routines do not currently work on smiller2 (linux box) but may work on a windows machine.  --> Use a different software (Python, IDL?) to plot the ephemeris on a map.

orbit_parms.f90: gets orbital elements, like get_elem, but then calculates other useful parameters.

NOTE:  Orbit_parms (short for orbit parameters) is used by itself to read TLE files and print out the contents (less drag coefficients) and parameters derived from the TLE data for the chosen satellite.

smiller@smiller2>>pgf90 synoptic_elem.f90 Time_Module.f90 Orbits_Module.f90

smiller@smiller2>>pgf90 get_elem.f90 Time_Module.f90

smiller@smiller2>>pgf90 Track.f90 Time_Module.f90

