#     - Observer: acts as our space camera, sees objects from any perspective (instead of just Earth's)

# import astropy
import astropy
from astropy import wcs
import astropy.table as tb
from astropy.wcs import WCS
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_body_barycentric, SkyCoord
# import others
import math
import numpy as np
import transformation as tr
import matplotlib.pyplot as pl

# *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** #

class Observer:
    '''moveable camera in ~space~'''
    def __init__(self, parallax_factor, z_add):
        '''
        parallax_factor: value between 0 and 1, AKA "earth shrink."
                  1 ~ same orbit size as earth: squiggly TNOs
                  0 ~ no orbiting: "straight" TNOs
        z_add: this value lifts our camera above the solar system
                  units: AU
                  note: the Observer moves along the ecliptic pole,
                  which is why adding to z impacts y. Hence, y_add
                  (see below)''' 
        self.p_factor = parallax_factor
        self.z_add = z_add + (self.p_factor*math.sin(23.5)) #cos(23.5) #remove p_factor
        self.y_add = z_add*(-1/2) - (self.p_factor*math.cos(23.5)) #sin(23.5) #remove p_factor
        self.inter = 70 # *** should match with Orbit's interval value!! again, room for improvement ***

    def getPositions(self, start, span):
        '''prepares a list of Observer positions over given time period'''
        interval = self.inter * span
        #prepare the timestamps which we give to astropy
        astrotimes = []
        for j in range(0,int(interval)):
            astrotimes.append(astropy.time.Time(((2016.0001+(2016-start))+((span/(interval)*j))), format='decimalyear', scale='tdb'))
        #get the positions of earth at these times
        earths = get_body_barycentric('earth', Time(astrotimes))
        #revise earth's positions to get our desired Observer positions
        obs_pos = np.array([self.p_factor*earths.x.value, 
                                       (self.p_factor*(earths.y.value) + self.y_add),
                                       (self.p_factor*(earths.z.value) + self.z_add)]).T
        return obs_pos
    
    def getVectors(self, obj_pos, start, span): #obj_pos: Orbit(#).getPath(start,span)[4]
        '''calculates vectors between the Observer and object(s) (obj_pos - obs_pos) 
        and returns the ra&dec values to visualize what the Observer sees'''
        interval = self.inter * span
        # retrieve Observer positionos
        obs_pos = self.getPositions(start,span)
        # subtract Observer positions from object positions
        vectors = [(obj_pos[i] - obs_pos) for i in range(len(obj_pos))]
        
        # convert from cartesian to ra and dec
        allras, alldecs = [], []
        for i in range(len(obj_pos)):    
            all_radec = SkyCoord(x=vectors[i][:,0], y=vectors[i][:,1], z=vectors[i][:,2], unit='AU', representation_type='cartesian')
            all_radec.representation_type = 'spherical'
            all_radec = np.array([all_radec.ra.value, all_radec.dec.value])
            # correct ra outliers (ex: some are recorded as -20 instead of 340)
            ras = np.where(all_radec[0] < 250, all_radec[0], all_radec[0] - 360)
            decs = all_radec[1]
            allras.append(ras)
            alldecs.append(decs)
        return allras, alldecs
    
    def wcsPoints(self, obj_pos, start, span, wcs):
        '''list of pixel coordinates of seen object(s) (for animation)'''
        rasdecs = self.getVectors(obj_pos, start, span)
        allrapix, alldecpix = [], []
        for i in range(len(obj_pos)):
            rapix, decpix = wcs.all_world2pix(rasdecs[0][i], rasdecs[1][i], 0)
            allrapix.append(rapix)
            alldecpix.append(decpix)
        return allrapix, alldecpix
    
    def plot(self, obj_pos, start, span, col):
        '''plots objects as a picture seen through the Observer's perspective
           col should be a list of colors. If color-coding, translate from colormap to color value beforehand'''
        rapoints = self.getVectors(obj_pos, start, span)[0]
        decpoints = self.getVectors(obj_pos, start, span)[1]
        for i in range(len(obj_pos)):
            pl.plot([px for px in rapoints[i]], [py for py in decpoints[i]], '-', color=col[i], alpha=1, markersize=4)
        #make sure axes don't stretch, and flip RA for astro purposes
        pl.gca().set_aspect('equal')
        
    def wcsPlot(self, obj_pos, start, span, col, axis, wcs):
        '''same as plot but for pixels / wcs'''
        rapoints = self.getVectors(obj_pos, start, span)[0]
        decpoints = self.getVectors(obj_pos, start, span)[1]
        for i in range(len(obj_pos)):
            rapix, decpix = wcs.all_world2pix(rapoints[i], decpoints[i], 0)
            axis.plot(rapix, decpix, '-', color=col[i], alpha=1, lw=1)
    
    def wcsPlotHead(self, obj_pos, start, span, col, axis, wcs, alpha_):
        '''plots objects on top of an image, using corresponding wcs
           plus: larger dot in front for "current position"'''
        rapoints = self.getVectors(obj_pos, start, span)[0]
        decpoints = self.getVectors(obj_pos, start, span)[1]
        for i in range(len(obj_pos)):
            rapix, decpix = wcs.all_world2pix(rapoints[i], decpoints[i], 0)
            axis.plot(rapix, decpix, '-', color=col[i], alpha=alpha_, lw=1)
            axis.plot(rapix[-1], decpix[-1], 'o', color=col[i],alpha=alpha_, markersize=7)
        
    def wcsStretch(self, points, xdata, ydata, obj_pos, start, span, wcs, frame):
        '''to visualize the shrinking of Earth's orbit (and unsquiggling of TNO Orbits), 
        this function appends point data to a pl.plot function defined in animation function'''
        self.p_factor = 1 - (0.01*frame)
        self.z_add = 0
        self.y_add = 0
        rapoints = self.getVectors(obj_pos, start, span)[0]
        decpoints = self.getVectors(obj_pos, start, span)[1]  
        for i in range(len(obj_pos)):
            rapix, decpix = wcs.all_world2pix(rapoints[i], decpoints[i], 0)     
            xdata[i].clear()      
            ydata[i].clear()
            xdata[i].append(rapix)  
            ydata[i].append(decpix)
            points[i].set_data(xdata[i], ydata[i])
        return points