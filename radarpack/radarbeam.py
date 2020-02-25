#
#  radarbeam.py
#
#       module for calculating geometry parameters and magnetic aspect
#       angle of radar targets monitored by any radar
#
#       use aspect_elaz or aspect_txty to calculate aspect angles of targets
#       specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08 as jrobeam.py
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#  history
#  - Aug29,2013 by P. Reyes
#    -Generate a module that accepts the lon,lat,h coordinates for the location
#    of any radar.
#    -flattening has been changed from 1/298.257        to 1./298.257223563
#    using the WGS84 reference in:
#    http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
#    - A new routine called enu2xyz to move a point from xr,yr,zr to some
#      direction east, north, up

from __future__ import print_function
from __future__ import division

def llh2xyz(latg,lon,h):
    # returns geocentric xyz coordinates (ECEF) in km of a target with
    # latitude       latg (rad) --- geodetic
    # longitude      lon (rad)
    # height         h (km above local ellipsoid)
    n=a_WGS / np.sqrt(1.-flatness*(2.-flatness) * np.sin(latg)**2.)
    # cartesian geocentric coordinates wrt Greenwich
    x=(n+h)*np.cos(latg)*np.cos(lon)
    y=(n+h)*np.cos(latg)*np.sin(lon)
    z=(n*(1.-eccentricity**2.)+h)*np.sin(latg)
    return x,y,z

def xyz2llh(x,y,z):
    # returns longitude 'lon', geodetic latitude 'lat', and height 'h'
    # of position (x,y,z) defined in geocentric coordinate system (ECEF)

    # on Oct23,2013 by P. Reyes, adding the .all() in order to support
    # arrays

    p=np.sqrt(x**2.+y**2.)
    lon=np.arctan2(y,x)
    lat=np.arctan2(z,p)
    latp=lat.copy()
    for i in range(10):
        n=a_WGS/np.sqrt(1.-flatness*(2-flatness)*np.sin(latp)**2.)
        h=p/np.cos(latp)-n
        lat=np.arctan(z/(p*(1.-n*eccentricity**2./(n+h))))
        if (abs(lat-latp)<3.*eps).all():
            n=a_WGS/np.sqrt(1.-flatness*(2.-flatness)*np.sin(lat)**2.)
            h=p/np.cos(lat)-n
            break
        latp=lat.copy()
    return lat,lon,h

def enu2xyz(xr,yr,zr,east,north,up):
    # moves a point from xr,yr,zr to x,y,z by moving into the direction
    # specified by east,north,up (enu) coordinates in km
    latg,lon,h = xyz2llh(xr,yr,zr)
    A = np.array([
        [-np.sin(lon), -np.sin(latg)*np.cos(lon), np.cos(latg)*np.cos(lon)],
        [ np.cos(lon), -np.sin(latg)*np.sin(lon), np.cos(latg)*np.sin(lon)],
        [           0,  np.cos(latg)            , np.sin(latg)]])
    x,y,z = np.dot(A,np.array([east,north,up]))+np.array([xr,yr,zr])
    return x,y,z

def blowup(arr,rep):            # enlarges square array arr by a factor rep
    dimarr = int(np.sqrt(np.size(arr)))
    dim = dimarr * rep
    arrl = np.zeros([dim,dim])
    for m in range(dimarr):
        for n in range(dimarr):
            arrl[m*rep:(m+1)*rep,n*rep:(n+1)*rep] = arr[m,n]
    return arrl

# --------------------------------------------------------------
import numpy as np
from pyigrf.igrf import igrf
igrf0 = igrf() # loading the latest coefficients

eps=np.finfo(float).eps         # float resolution
deg=np.pi/180.                  # to express angles in degree values

# WGS84 constants
# reference:
# http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf

a_WGS=6378.137         # equatorial radius WGS 84 (semi-major axis) in km
#flatness=1/298.257
flatness = 1./298.257223563  # flatenning
b_WGS=a_WGS*(1.-flatness)    # WGS polar radius (semi-minor axis) in km
eccentricity=np.sqrt(a_WGS**2-b_WGS**2)/a_WGS

# ------------ radar specifications -------------------------
class Antenna_Beam:
    def __init__(self, Bname, antenna_specs, wide_beam_pattern, beam_pattern,
            Txdir, Tydir, el, az,
            TxComb, RxComb, dir_mode, dB_threshold,Rx_order):
        self.Bname = Bname
        self.Txdir = Txdir
        self.Tydir = Tydir
        self.el = el
        self.az = az
        self.wide_beam_pattern = wide_beam_pattern
        self.wide_beam_N_pts = antenna_specs.wide_beam_N_pts
        #N = self.wide_beam_N_pts
        self.beam_pattern = beam_pattern
        self.narrow_beam_N_pts = antenna_specs.narrow_beam_N_pts
        self.ThetaX = antenna_specs.ThetaX
        self.ThetaY = antenna_specs.ThetaY
        self.rotation_xy = antenna_specs.rotation_xy
        self.ThetaXprime = antenna_specs.ThetaXprime
        self.ThetaYprime  = antenna_specs.ThetaYprime
        self.wide_ThetaX = antenna_specs.wide_ThetaX
        self.wide_ThetaY = antenna_specs.wide_ThetaY
        self.TxComb = TxComb
        self.RxComb = RxComb
        self.dir_mode = dir_mode
        self.dB_threshold = dB_threshold
        self.Rx_order = Rx_order

    def get_basic_info(self):
        return {'Bname':self.Bname,'el':self.el,'az':self.az,
                'Txdir':self.Txdir, 'Tydir':self.Tydir,
                'TxComb':self.TxComb,'RxComb':self.RxComb,
                'dir_mode':self.dir_mode,'dB_threshold':self.dB_threshold,
                'Rx_order':self.Rx_order
                }

    def print_summary(self):
        print ( "%8s"%self.Bname +
                ": (Theta_x,Theta_y)=(%11.8f,%11.8f)"%(self.Txdir,self.Tydir) +
                ", (el,az)=(%6.3f,%8.3f)"%(self.el, self.az))

class radarspecs:
    """Will contain radar coordinates and coordinate conversions
    saved locations:
    JRO : lat: -11.947917 , lon: -76.872306, h0: 0.463 km
    JRO_GE : as zoom in with GoogleEarth to the center of the antenna.
    IRIS@ROI
    ALTAIR
    IRIS@URBANA
    """
    def __init__(self,lat0=None,lon0=None,h0=None,location=None,fload='',dec=0.,ha=0.,
            rotation_xy=0.):
        self.ph_arr_d = {}
        self.Uphs = {}
        self.ph_exp_d = {}
        self.comb4q_arr = {}
        self.wide_comb4q_arr = {}
        self.beam_pattern = {}
        self.rotation_xy = rotation_xy
        self.next_order = 0
        if lat0 != None and lon0 != None and h0 != None:
            self.lat0 = lat0 * deg
            self.lon0 = lon0* deg
            self.h0 = h0 # local height in km above reference ellipsoid
            if location == None:
                self.location = "Unknown"
            else:
                self.location = location
            self.dec = dec
            self.ha = ha
        elif location!=None:
            self.location = location.upper()
            if self.location == "JRO":
                # geodetic, the usual map or GPS latitude
                self.lat0 = -11.947917 * deg
                self.lon0 = -76.872306 * deg
                self.h0 = 0.463 # local height in km above reference ellipsoid
                self.dec=-12.88*deg
                self.ha=-(4.+37./60.)*deg # on-axis direction at JRO
            elif self.location == "JRO_GE":
                # gusing google earth to the center of the Antenna
                # -11.9514944444 = -(11.+57./60.+5.38/3600.) # 11deg57'5.38"S
                self.lat0 = -11.9514944444 * deg
                # -76.8743916667#-(76.+52./60.+27.81/3600.) #  76deg52'27.81"W
                self.lon0 = -76.8743916667 * deg
                self.h0 = 0.463 # local height in km above reference ellipsoid
                self.dec=-12.88*deg
                self.ha=-(4.+37./60.)*deg # on-axis direction at JRO
            elif self.location == "IRIS@ROI":
                #  9.39794444444 = (9.+23./60.+52.6/3600.) # 9deg23'52.60"N
                self.lat0 = 9.39794444444 * deg
                #  167.469166667 = (167.+28./60.+9./3600.) # 167deg28'9.00"E
                self.lon0 = 167.469166667 * deg
                self.h0 = 0.012
                self.dec= 8.*deg
                self.ha= 0. #
            elif self.location == "ALTAIR":
                #  9.39794444444 = (9.+23./60.+43.5/3600.) # 9deg23'43.50"N
                self.lat0 = 9.39541666667 * deg
                #  167.469166667 = (167.+28./60.+45.6/3600.) # 167deg28'45.60"E
                self.lon0 = 167.479333333 * deg
                self.h0 = 0.012
                self.dec= 30.*deg
                self.ha= 0. #
            elif self.location == "IRIS@URBANA":
                # 40.16683888888889 = (40.+10./60.+0.62/3600.) #40deg10'0.62"N
                self.lat0 = 40.16683888888889 * deg
                #-88.1586 = -(88.+9./60.+30.96/3600.) #88deg9'30.96"W
                self.lon0 = -88.1586 * deg
                self.h0 = 0.221
                self.dec= 30.*deg
                self.ha= 0. #
            elif self.location == "SANYA":
                self.lat0 = 18.34 * deg
                self.lon0 = 109.62 * deg
                self.h0 = 0.035
                self.dec= - 30.*deg
                self.ha= 0. #
        elif fload!='':
            fp = np.load(fload, encoding='bytes',allow_pickle=True)
            self.comb4q_arr = fp['comb4q_arr'][0]
            self.wide_comb4q_arr = fp['wide_comb4q_arr'][0]
            self.h0 = fp['h0'][0]
            self.aspect_ranges = fp['aspect_ranges']
            self.aspects = fp['aspects']
            self.B_xyz = fp['B_xyz']
            self.month = fp['month'][0]
            self.yyyymmdd = fp['yyyymmdd'][0]
            self.year = fp['year'][0]
            self.ph_exp_d = fp['ph_exp_d'][0]
            self.day = fp['day'][0]
            self.fyear = fp['fyear'][0]
            self.month_name = fp['month_name'][0]
            self.location = fp['location'][0]
            self.ph_arr_d = fp['ph_arr_d'][0]
            self.ThetaX = fp['ThetaX']
            self.ThetaY = fp['ThetaY']
            self.rotation_xy = fp['rotation_xy']
            self.ThetaXprime = fp['ThetaXprime']
            self.ThetaYprime = fp['ThetaYprime']
            self.wide_ThetaY = fp['wide_ThetaY']
            self.wide_ThetaX = fp['wide_ThetaX']
            self.lat0 = fp['lat0'][0]
            self.wide_beam_N_pts = fp['wide_beam_N_pts'][0]
            self.narrow_beam_N_pts = fp['narrow_beam_N_pts'][0]
            self.doy = fp['doy'][0]
            self.Uphs = fp['Uphs'][0]
            self.lon0 = fp['lon0'][0]
            if 'declination' in fp:
                self.dec = fp['declination'][0]
            else:
                self.dec = -12.88*deg
            if 'hour_angle' in fp:
                self.ha = fp['hour_angle'][0]
            else:
                self.ha = -(4.+37./60.)*deg # on-axis direction at JRO
            wide_beam_patterns = fp['wide_beam_patterns'][0]
            beam_patterns = fp['beam_patterns'][0]
            beams_info = fp['beams_info'][0]
            for Bname in beams_info:
                beam = beams_info[Bname]
                wide_beam_pattern = wide_beam_patterns[Bname]
                beam_pattern = beam_patterns[Bname]
                BnameAntenna = Antenna_Beam(Bname,self, wide_beam_pattern,
                        beam_pattern,
                        beam[b'Txdir'], beam[b'Tydir'], beam[b'el'],beam[b'az'],
                        beam[b'TxComb'], beam[b'RxComb'], beam[b'dir_mode'],
                        beam[b'dB_threshold'],beam[b'Rx_order'])
                self.beam_pattern.update({Bname:BnameAntenna})
            fp.close()
        else:
            # gusing google earth to the center of the Antenna
            # -11.9514944444 = -(11.+57./60.+5.38/3600.) # 11deg57'5.38"S
            self.lat0 = -11.9514944444 * deg
            # -76.8743916667#-(76.+52./60.+27.81/3600.) #  76deg52'27.81"W
            self.lon0 = -76.8743916667 * deg
            self.h0 = 0.463 # local height in km above reference ellipsoid
            self.dec=-12.88*deg
            self.ha=-(4.+37./60.)*deg # on-axis direction at JRO

            self.location = 'JRO'

        x0,y0,z0 = llh2xyz(self.lat0,self.lon0,self.h0)
        self.xyz0 = np.array([x0,y0,z0])
        xy0 = np.array([x0,y0])
        p0 = np.sqrt(np.dot(xy0,xy0))

        # unit vectors from jro
        self.east0 = np.array([-y0,x0,0])/p0
        # zenith and north directions wrt local ellipsoid
        self.zenith0 = np.array([np.cos(self.lat0) * np.cos(self.lon0),
                            np.cos(self.lat0) * np.sin(self.lon0),
                            np.sin(self.lat0)])
        self.north0 = np.cross(self.zenith0,self.east0)

        # orthonormal basis vectors including the jro on-axis direction
        self.uo = np.array([np.cos(self.dec) * np.cos(self.ha/4. + self.lon0), # on axis
                    np.cos(self.dec) * np.sin(self.ha/4. + self.lon0), np.sin(self.dec)])
        self.ux = np.cross(self.zenith0,self.uo)
        # along the building to the right
        self.ux = self.ux / np.sqrt(np.dot(self.ux,self.ux))
        # away from the building into the valley
        self.uy = np.cross(self.uo,self.ux)
        #rotation of xy plane:
        self.uDn = self.ux.copy()
        self.uUp = self.uy.copy()

    def get_beam_pattern(self, Bname):
        return self.beam_pattern[Bname]

    def get_beam_names(self):
        return list(self.beam_pattern.keys())

    def get_beam_names_ordered(self):
        oindex = []
        unordered = self.get_beam_names()
        for bname in unordered:
            oindex += [self.beam_pattern[bname].Rx_order]
        ordered = []
        for i in oindex:
            ordered +=[unordered[i]]
        return ordered

    def print_beams_summary(self,ordered=False):
        if ordered:
            bnames = self.get_beam_names_ordered()
        else:
            bnames = self.get_beam_names()
        for Bname in bnames:
            self.beam_pattern[Bname].print_summary()

    def print_html_stearing_phases(self):

        """
        Returns html code for the stearing phases table.
        """
        out = "<html><body>\n"
        out += '<table border="1" width="400" cellspacing="0" style="text-align: center;">\n'
        for i in range(4):
            if i==3:
                comm_style = "border-bottom: 5px double;"
            else:
                comm_style = ""
            out += "<tr>\n"
            for j in range(4):
                up = self.ph_arr_d['N']['up'][::12,::12][i,j]
                dn = self.ph_arr_d['N']['dn'][::12,::12][i,j]
                if j<3:
                    css_style = comm_style
                else:
                    css_style = comm_style + "border-right: 5px double;"
                out += '  <td style="%s">%g<br>%g</td>\n'%(css_style,up,dn)
            for j in range(4):
                up = self.ph_arr_d['E']['up'][::12,::12][i,j]
                dn = self.ph_arr_d['E']['dn'][::12,::12][i,j]
                out += '  <td style="%s">%g<br>%g</td>\n'%(comm_style,up,dn)
            out += "</tr>\n"
        for i in range(4):
            comm_style = ""
            out += "<tr>\n"
            for j in range(4):
                up = self.ph_arr_d['W']['up'][::12,::12][i,j]
                dn = self.ph_arr_d['W']['dn'][::12,::12][i,j]
                if j<3:
                    css_style = comm_style
                else:
                    css_style = comm_style + "border-right: 5px double;"
                out += '  <td style="%s">%g<br>%g</td>\n'%(css_style,up,dn)
            for j in range(4):
                up = self.ph_arr_d['S']['up'][::12,::12][i,j]
                dn = self.ph_arr_d['S']['dn'][::12,::12][i,j]
                out += '  <td style="%s">%g<br>%g</td>\n'%(comm_style,up,dn)
            out += "</tr>\n"
        out += "</table>\n"
        out += "<br>Ues:\n"
        out += '<table border="1" width="400" cellspacing="0" style="text-align: center;">\n'
        out += "<tr>\n"
        for i,antenna in enumerate(['W','N','E','S']):
            comm_style = "text-align: center;"
            pol = ["D","U"][(i+1)%2]
            out += '  <td>%s<sub>%s</sub></td>\n'%(antenna,pol)
            pol = ["D","U"][(i)%2]
            out += '  <td>%s<sub>%s</sub></td>\n'%(antenna,pol)
        out += "</tr>\n"
        out += "<tr>\n"
        for i,antenna in enumerate(['W','N','E','S']):
            comm_style = "text-align: center;"
            pol = ["dn","up"][(i+1)%2]
            out += '  <td>%.2f</sub></td>\n'%(
                                            self.Uphs[antenna][pol])
            pol = ["dn","up"][(i)%2]
            out += '  <td>%.2f</sub></td>\n'%(
                                            self.Uphs[antenna][pol])
        out += "</tr>\n"

        out += "</table></body></html>\n"
        return out

    def save_beams_info(self, expname,folder=''):
        fname = (folder+'/'+'Beams_'+expname+
                '_%.4d.%.2d.%.2d.npz'%(self.year,self.month,self.day))
        print('saving... '+fname)
        beams_info = {}
        beam_patterns = {}
        wide_beam_patterns = {}
        for Bname in self.get_beam_names():
            beams_info.update({Bname:self.beam_pattern[Bname].get_basic_info()})
            beam_patterns.update({Bname:self.beam_pattern[Bname].beam_pattern})
            wide_beam_patterns.update({Bname:self.beam_pattern[Bname].wide_beam_pattern})
        np.savez_compressed(fname,
                beams_info = [beams_info],
                beam_patterns = [beam_patterns],
                wide_beam_patterns = [wide_beam_patterns],
                ThetaX = self.ThetaX,
                ThetaY = self.ThetaY,
                rotation_xy = self.rotation_xy,
                ThetaXprime = self.ThetaXprime,
                ThetaYprime = self.ThetaYprime,
                wide_ThetaX = self.wide_ThetaX,
                wide_ThetaY = self.wide_ThetaY,
                wide_beam_N_pts = [self.wide_beam_N_pts],
                narrow_beam_N_pts = [self.narrow_beam_N_pts],
                declination = [self.dec],
                hour_angle = [self.ha],
                ph_arr_d = [self.ph_arr_d],
                Uphs = [self.Uphs], lon0 = [self.lon0], lat0 = [self.lat0],
                h0 = [self.h0], location=[self.location], day = [self.day],
                doy = [self.doy], fyear = [self.fyear], month = [self.month],
                month_name = [self.month_name], year = [self.year],
                comb4q_arr= [self.comb4q_arr],
                wide_comb4q_arr= [self.wide_comb4q_arr],
                yyyymmdd = [self.yyyymmdd],
                aspect_ranges = self.aspect_ranges,
                ph_exp_d = [self.ph_exp_d],
                aspects = self.aspects,
                radar_east = self.east0,
                radar_north = self.north0,
                radar_zenith = self.zenith0,
                radar_xyz = self.xyz0,
                radar_ux = self.ux,
                radar_uy = self.uy,
                radar_uo = self.uo,
                radar_ux_prime = self.ux_prime,
                radar_uy_prime = self.uy_prime,
                B_xyz = self.B_xyz,
                ordered_bnames = self.get_beam_names_ordered())

    def set_beam_pattern(self, Bname, TxComb, RxComb, dir_mode='max',
            dB_threshold=-6, order=None):
        """
        Inputs:
            Bname: string containing the name to give to this pattern
            TxComb: string containing the antenna name, combination of
                JRO quarters corresponding to a transmitting antenna,
                created usin the self.create_4quarters routine.
            RxComb: string containing the antenna name, combination of
                JRO quarters corresponding to a transmitting antenna,
                created usin the self.create_4quarters routine.
            dir_mode: A string specifying the method of calculating
                the direction of the beam. The default is 'max' which
                finds the maximum value of the beam pattern
                'weighted_avg' is also available. This is the weighted average
                and it calculates the first moment of the beam pattern
                over the values db_threshold to 0 dB.
            dB_threshold: Used for 'weighted_avg' mode in dir_mode.
            order: The order of the receiving channel starting at 0.
                If order == None, then the order will be given as an
                increasing sequence each time this function is called.
        """
        wide_logarg = abs(self.wide_comb4q_arr[TxComb]['FFTa'])**2. * \
                abs(self.wide_comb4q_arr[RxComb]['FFTa'])**2.
        wide_logarg[wide_logarg<=0] = np.nan  # to avoid errors when calculating log10
        wide_beam_pa = 10. * np.log10(wide_logarg)
        wide_beam_pa -= np.nanmax(wide_beam_pa)  # normalizing for max value = 0 dB
        Txdir, Tydir = self.calc_Beam_Direction(wide_beam_pa, method=dir_mode,
                dB_threshold = dB_threshold)
        el, az = self.dir_cos2el_az(Txdir, Tydir)

        logarg = abs(self.comb4q_arr[TxComb]['FFTa'])**2. * \
                abs(self.comb4q_arr[RxComb]['FFTa'])**2.
        logarg[logarg<=0] = np.nan  # to avoid errors when calculating log10
        beam_pa = 10. * np.log10(logarg)
        beam_pa -= np.nanmax(beam_pa)  # normalizing for max value = 0 dB

        if type(order)==type(None):
            current_order = self.next_order
            self.next_order += 1
        else:
            current_order = order
        BnameAntenna = Antenna_Beam(Bname,self,wide_beam_pa,beam_pa,
                Txdir,Tydir,el,az,
                TxComb, RxComb, dir_mode, dB_threshold, current_order)
        self.beam_pattern.update({Bname:BnameAntenna})

        return BnameAntenna

    def calc_Beam_Direction(self, Bpattern, method='max', dB_threshold=-6):
        """Bpattern is the beam pattern with respect to the antenna
        where the coordinates of the 2D array are the direction
        cosines self.ThetaX, and self.ThetaY.
        The values are in dBs and has been normalized with the
        maximum being 0dB.
        If method = 'weighted_avg', then the input dB_threshold is used
        """
        b = Bpattern
        if method == 'max':
            [i,j]=np.unravel_index(np.nanargmax(b), b.shape)
            ThetaX_dir, ThetaY_dir = self.wide_ThetaX[i,j], self.wide_ThetaY[i,j]
        elif method == 'weighted_avg':
            b[np.isnan(b)] = dB_threshold - 1
            th_select = b >= dB_threshold
            lin_b = 10 ** (b/10.)
            b_avg = np.trapz(lin_b[th_select])
            ThetaX_dir = np.trapz((lin_b * self.wide_ThetaX)[th_select]) / b_avg
            ThetaY_dir = np.trapz((lin_b * self.wide_ThetaY)[th_select]) / b_avg
        return ThetaX_dir, ThetaY_dir

    def create_4quarters(self, CombName, N, E, W, S):
        quarters ={}
        ons = np.ones([48,48]); offs = ons * 0 # amplitude distributions
        col0 = np.zeros([48,1]); row0 = np.zeros([1,97])
        for q in 'NEWS':
            if locals()[q] == 'ons':
                quarters.update({q:ons})
            elif locals()[q] == 'offs':
                quarters.update({q:offs})
            else:
                name,pol = locals()[q].split('.')
                quarters.update({q:self.ph_exp_d[name][pol]})

        qN = quarters['N']
        qE = quarters['E']
        qW = quarters['W']
        qS = quarters['S']

        N_pts = self.wide_beam_N_pts
        N2_pts = self.narrow_beam_N_pts

        temp_wide_FFTa = np.fft.fftshift(np.fft.fft2(
            np.asarray(np.bmat('qN,col0,qE;row0;qW,col0,qS')), s=(N_pts,N_pts)))

        #temp_FFTa = temp_wide_FFTa[N_pts/2.-N_pts/16.:N_pts/2.+N_pts/16.,
        #                                N_pts/2.-N_pts/16.:N_pts/2.+N_pts/16.]
        from scipy import interpolate
        f_real = interpolate.interp2d(self.wide_ThetaX[0,:], self.wide_ThetaY[:,0],
                temp_wide_FFTa.real, kind='linear')
        f_imag = interpolate.interp2d(self.wide_ThetaX[0,:], self.wide_ThetaY[:,0],
                temp_wide_FFTa.imag, kind='linear')

        temp_FFTa = np.empty_like(self.ThetaX, dtype="complex128")
        for row,TX_row in enumerate(self.ThetaX):
            for col,TX in enumerate(TX_row):
                TY = self.ThetaY[row,col]
                temp_FFTa.real[row,col] = f_real(TX,TY)
                temp_FFTa.imag[row,col] = f_imag(TX,TY)


        if self.wide_comb4q_arr.has_key(CombName):
            # This will overwrite the quarter combination
            self.wide_comb4q_arr[CombName]['FFTa'] = temp_wide_FFTa
            self.wide_comb4q_arr[CombName]['N'] = N
            self.wide_comb4q_arr[CombName]['E'] = E
            self.wide_comb4q_arr[CombName]['W'] = W
            self.wide_comb4q_arr[CombName]['S'] = S
        else:
            # This will create a quarter combination
            self.wide_comb4q_arr.update({CombName:
                {'FFTa':temp_wide_FFTa,
                    'N':N,'E':E,'W':W,'S':S}
                })

        if self.comb4q_arr.has_key(CombName):
            # This will overwrite the quarter combination
            self.comb4q_arr[CombName]['FFTa'] = temp_FFTa
            self.comb4q_arr[CombName]['N'] = N
            self.comb4q_arr[CombName]['E'] = E
            self.comb4q_arr[CombName]['W'] = W
            self.comb4q_arr[CombName]['S'] = S
        else:
            # This will create a quarter combination
            self.comb4q_arr.update({CombName:{'FFTa':temp_FFTa,'N':N,'E':E,'W':W,'S':S} })

    def __update_exp_array(self,Aname, Apol):
        val = np.exp(2j * np.pi *
                (self.ph_arr_d[Aname][Apol] + self.Uphs[Aname][Apol])/4.)
        if self.ph_exp_d.has_key(Aname):
            self.ph_exp_d[Aname].update({Apol:val})
        else:
            self.ph_exp_d.update({Aname:{Apol:val}})

    def copy_phase_arr(self, Aname, Apol, newAname, newApol, Uph):
        """
        Copy the Aname, Apol to newAname, newApol
        """
        if self.ph_arr_d.has_key(newAname):
            self.ph_arr_d[newAname].update({newApol:self.ph_arr_d[Aname][Apol]})
            self.Uphs[newAname].update({newApol:Uph})
        else:
            self.ph_arr_d.update({newAname:{newApol:self.ph_arr_d[Aname][Apol]}})
            self.Uphs.update({newAname:{newApol:Uph}})
        self.__update_exp_array(newAname, newApol)

    def create_beam_aspect_angle_grid(self, N,N2, yyyymmdd, ranges_km, xp_overs_fact=1.,
            yp_overs_fact=1., verbose=False ):
        """
        Create grid of direction cosines ThetaX, ThetaY
        yyyymmdd : year, month and day separated by dots, e.g 2014.10.12
                   for October 12, 2014
        xp_overs_fact: oversample factor in the x prime axis(rotated by rotation_xy)
        yp_overs_fact: oversample factor in the y prime axis(rotated by rotation_xy)
        """
        from calendar import isleap, timegm
        from time import gmtime, strftime
        self.yyyymmdd = yyyymmdd
        ts = gmtime(timegm(np.array((yyyymmdd+'.0.0.0').split('.'),dtype=int)))
        self.year = ts.tm_year
        self.month = ts.tm_mon
        self.month_name = strftime('%b',ts)
        self.day = ts.tm_mday
        self.doy = ts.tm_yday
        yeardays = 366. if isleap(self.year) else 365.
        self.fyear = self.year + self.doy/yeardays
        self.wide_beam_N_pts = N
        self.narrow_beam_N_pts = N2
        [self.wide_ThetaY,self.wide_ThetaX] = np.mgrid[-N/2.:N/2.,-N/2.:N/2.]
        self.wide_ThetaY = - self.wide_ThetaY / (N / 2.)
        self.wide_ThetaX =   self.wide_ThetaX / (N / 2.)

        [self.ThetaYprime,self.ThetaXprime] = np.mgrid[
                -N2/2.*yp_overs_fact:N2/2.*yp_overs_fact,
                -N2/2.*xp_overs_fact : N2/2.*xp_overs_fact]
        self.ThetaYprime = - self.ThetaYprime / (N / 2. * yp_overs_fact)
        self.ThetaXprime =   self.ThetaXprime / (N / 2. * xp_overs_fact)

        self.ThetaX = np.empty_like(self.ThetaXprime)
        self.ThetaY = np.empty_like(self.ThetaYprime)

        if self.rotation_xy == "opt_rot":
            rtarget = 300. * self.uo # from Antenna
            xyz = self.xyz0 + rtarget # from center of Earth
            r,lat,lon,aspect,B = self.aspect_angle(self.fyear,xyz)
            mag_B = np.sqrt(np.dot(B,B))
            mag_uy = np.sqrt(np.dot(self.uy,self.uy))
            opt_rot = np.arccos(np.dot(self.uy,B)/mag_B/mag_uy)
            self.rotation_xy = opt_rot
        cosrot = np.cos(self.rotation_xy)
        sinrot = np.sin(self.rotation_xy)
        self.ux_prime = cosrot * self.uDn + sinrot * self.uUp
        self.uy_prime = -sinrot * self.uDn + cosrot * self.uUp

        for row,Xp_row in enumerate(self.ThetaXprime):
            for col,Xp in enumerate(Xp_row):
                Yp = self.ThetaYprime[row,col]
                self.ThetaX[row,col] = cosrot * Xp - sinrot * Yp
                self.ThetaY[row,col] = sinrot * Xp + cosrot * Yp

        self.aspect_ranges = ranges_km
        self.aspects = np.empty(((len(ranges_km),)+self.ThetaX.shape))
        self.B_xyz = np.empty(((3,len(ranges_km),)+self.ThetaX.shape))
        for i, current_range in enumerate(ranges_km):
            if verbose:
                print("calculating B at %.2f. range %d out of %d."%(
                        current_range,i+1,len(ranges_km)))
            for row,Xp_row in enumerate(self.ThetaX):
                for col,Xp in enumerate(Xp_row):
                    [r,lon,lat,dec,ha,aspect,B] = self.aspect_txty(self.fyear,
                            current_range, self.ThetaX[row,col], self.ThetaY[row,col])
                    self.aspects[i,row,col] = aspect * 180. / np.pi
                    self.B_xyz[:,i,row,col] = B

    def add_phase_arr(self, Aname, Apol, rep, phs, Uph=0.,add_to_cables=0.):
        """
        Aname: Name of the Antenna e.g. 'N','S','W','E'
        Apol: Polarizarion of the Antenna, e.g 'up', 'dn'
        rep: Number of dipoles per unit block, e.g. 12
        phs: nxn phases, e.g. 4x4 phases
        """
        phs = np.array(phs) + add_to_cables
        phs[phs>5] -= 4
        phs[phs<2] += 4
        if self.ph_arr_d.has_key(Aname):
            self.ph_arr_d[Aname].update({Apol:blowup(phs,rep)})
            self.Uphs[Aname].update({Apol:Uph})
        else:
            self.ph_arr_d.update({Aname:{Apol:blowup(phs,rep)}})
            self.Uphs.update({Aname:{Apol:Uph}})
        self.__update_exp_array(Aname, Apol)

    def locations(self):
        return ["JRO","JRO_GE","IRIS@ROI","ALTAIR","IRIS@URBANA"]

    def dec_ha2el_az(self,dec,ha):
        # returns elevation and azimuth angles of a radar beam
        # with respect to local tangent plane.
        # the beam is specified by:
        #		declination dec (deg)
        #		hour angle  ha  (min)
        # with respect to radar location at longitude lon0 and height h0
        # above reference ellipsiod at geodetic latitude lat0

        lat = dec * deg                 # on celestial sphere
        lon = 2. * np.pi * (ha/(24.*60.))
        lon = lon + self.lon0                # on celestial sphere
        vec = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon),
                        np.sin(lat)])
        hor = vec - np.dot(vec, self.zenith0) * self.zenith0
        hor = hor / np.sqrt(np.dot(hor, hor))
        el = np.arccos(np.dot(hor, vec)) / deg
        north = np.dot(hor, self.north0)
        east = np.dot(hor, self.east0)
        az = np.arctan2(east, north)/deg
        return el,az

    def xyz2dec_ha(self,vec):
        """    declination and hour angle in target direction used to
        describe radar beam direction at JRO, corresponding to latitude and
        relative longitude of the beam-spot on the celestial sphere, corresponds
        to rr->\infty, in which case:"""
        vec = vec/np.sqrt(np.dot(vec,vec))
        p = np.sqrt(vec[0]**2.+vec[1]**2.)
        dec = np.arctan2(vec[2],p)/deg                           # in degrees
        ha = (np.arctan2(vec[1],vec[0]) - self.lon0)*(24./(2.*np.pi))*60.    # in minutes
        return dec,ha

    def dir_cos2xyz(self, Theta_x, Theta_y):
        """    Direction Cosines Theta_x and Theta_y describe radar beam
        direction with respect to the antenna of a radar.
        input: direction cosines
        output: unit vector in ECEF"""
        Theta_z = np.sqrt(1. - (Theta_x ** 2. + Theta_y ** 2.))
        return  Theta_x * self.ux + Theta_y * self.uy + Theta_z * self.uo

    def dir_cos2el_az(self, Theta_x, Theta_y):
        """    Direction Cosines Theta_x and Theta_y describe radar beam
        direction with respect to the antenna of a radar.
        input: direction cosines
        output: elevation and azimuth angles"""
        beam_hat = self.dir_cos2xyz(Theta_x, Theta_y)
        vert_comp = np.dot(beam_hat, self.zenith0) * self.zenith0
        horiz_comp = beam_hat - vert_comp
        horiz_hat = horiz_comp / np.sqrt(np.dot(horiz_comp,horiz_comp))
        el = np.arccos(np.dot(horiz_hat, beam_hat)) / deg
        north_comp = np.dot(horiz_hat, self.north0)
        east_comp  = np.dot(horiz_hat, self.east0)
        az = np.arctan2(east_comp, north_comp) / deg
        return el,az

    def aspect_angle(self,year,xyz,radar_xyz=None):
        # returns the magnetic aspect angle (rad) of a target with
        # geocentric vector xyz defined in geocentric coordinates

        r   = np.sqrt(np.dot(xyz,xyz))    # from center of earth
        p   = np.sqrt(xyz[0]**2. + xyz[1]**2.)
        lat = np.arctan2(xyz[2],p)      # geocentric latitude (theta)
        lon = np.arctan2(xyz[1],xyz[0]) # geocentric or geodetic longitude (phi)
        radial = xyz/r;     # radial direction of target from center of Earth
        east   = np.array([-xyz[1],xyz[0],0.])/p  # geocentric or geodetic east
        north  = -np.cross(east,radial)   # geocentric north unit vector
        if type(radar_xyz) == type(None):
            rr = xyz - self.xyz0     # target vector from radar antenna
        else:
            rr = xyz - radar_xyz
        u_rr = rr / np.sqrt(np.dot(rr,rr))   # unit vector from radar to target
        [bX,bY,bZ,bB] = igrf0.igrf_B(year, r - igrf0.a, lon/deg, lat/deg)
        bfield = np.array([bX,bY,bZ])
        B = bX*north + bY*east - bZ*radial
        u_B = B / np.sqrt(np.dot(B,B))
        aspect = np.arccos(np.dot(u_B, u_rr))
        return r,lat,lon,aspect,B

    def aspect_txty(self,year,rr,tx,ty):
        # returns magnetic aspect angle and geocentric coordinates of a target
        # tracked by jro at
        # range rr (km)
        # tx along jro building
        # ty into the building

        tz = np.sqrt(1.-tx**2.-ty**2.)
        #geocentric coordinates of target
        xyz = self.xyz0 + rr*(tx*self.ux + ty*self.uy + tz*self.uo)

        [r,lat,lon,aspect,B] = self.aspect_angle(year,xyz)
        [dec,ha] = self.xyz2dec_ha(xyz - self.xyz0)
        return r,lon,lat,dec,ha,aspect,B

    def aspect_elaz(self,year,rr,el,az):
        # returns magnetic aspect angle and geocentric coordinates of a target
        # tracked by jro at
        # range       rr (km)
        # elevation   el (rad above local tangent plane to ellipsoid)
        # azimuth     az (rad east of local north)

        tx = np.cos(el) * np.sin(az)		# direction cosines wrt east and north
        ty = np.cos(el) * np.cos(az)
        tz = np.sin(el)
        #geocentric coordinates of target :
        xyz = self.xyz0 + rr*(tx * self.east0 + ty*self.north0+tz*self.zenith0)

        [r,lat,lon,aspect,B] = self.aspect_angle(year,xyz)
        [dec,ha] = self.xyz2dec_ha(xyz - self.xyz0)
        return r,lon,lat,dec,ha,aspect,B
