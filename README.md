# radarpack
module for calculating geometry parameters and magnetic aspect angle of radar targets monitored by any radar

Installation
------------
After cloning this repository:

    pip install .

Examples
--------

    >>> from radarpack import radarbeam
    >>> import numpy as np
    >>> 
    >>> pfisr = radarbeam.radarspecs(location = 'PFISR')
    >>> # Finding a target at:
    >>> target_range = 300 # km
    >>> target_el = 70 * np.pi/180
    >>> target_az = -30 * np.pi/180
    >>> target_xyz = pfisr.elaz2xyz(target_range, target_el, target_az)
    >>> target_lat, target_lon, target_h = radarbeam.xyz2llh(*target_xyz)
    >>> 
    >>> print(f"PFISR site altitude = {pfisr.h0} km." )
    >>> print(f"PFISR site longitude = {pfisr.lon0*180/np.pi} deg." )
    >>> print(f"PFISR site latitude = {pfisr.lat0*180/np.pi} deg." )
    >>> print("PFISR ECEF xyz=",pfisr.xyz0) # print PISR ECEF coordinates 
    >>> 
    >>> print(f"target range = {target_range} km." )
    >>> print(f"target elevation = {target_el*180/np.pi} deg." )
    >>> print(f"target azimuth = {target_az*180/np.pi} deg." )
    >>> 
    >>> print("target ECEF xyz=",target_xyz)
    >>> print(f"target altitude = {target_h} km." )
    >>> print(f"target longitude = {target_lon*180/np.pi} deg." )
    >>> print(f"target latitude = {target_lat*180/np.pi} deg." )
    
    PFISR site altitude = 0.213 km.
    PFISR site longitude = -147.47104 deg.
    PFISR site latitude = 65.12992 deg.
    PFISR ECEF xyz= [-2267.91913085 -1446.43590632  5764.00992906]
    target range = 300 km.
    target elevation = 70.0 deg.
    target azimuth = -29.999999999999996 deg.
    target ECEF xyz= [-2327.4947121  -1423.58305252  6057.14555841]
    target altitude = 282.9096602325635 km.
    target longitude = -148.54848033844206 deg.
    target latitude = 65.88933220011259 deg.
