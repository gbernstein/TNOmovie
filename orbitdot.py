# Contents: Dot class, Orbit class, getPlanetPos function
#     - Dot: catalogued detections of a TNO, visualized as flashes or dots on a plot
#     - Orbit: the path of a given TNO (squiggles from Earth's POV, orbits when using Observer)
#     - getPlanetPos: uses astropy to retrieve Solar System planet coordinates

# import astropy
import astropy
from astropy import wcs
import astropy.table as tb
from astropy.wcs import WCS
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_body_barycentric, SkyCoord
# import others
import numpy as np
import transformation as tr
import matplotlib.pyplot as pl

# relevant FITS files
# - table with orbital elements, IDs, # of detections, distance of each TNO
t = tb.Table.read('dataa/orbits_y6_preview.fits')
tnoinfo = t['ORBITID','OLDID','NDETECT','a','e','i','Omega','omega','T_p','DISTANCE']
tnoinfo['DISTANCE'].name = 'd'
tnoinfo['T_p'] = tnoinfo['T_p'] + 2000
# - table with all recorded TNO detections
tnodetect = tb.Table.read('dataa/orbits_y6_preview.fits', 2)
tnodetect['RA_DET'].name = 'RA'
tnodetect['DEC_DET'].name = 'DEC'

# *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** #

class Dot:
    '''TNO detections to be drawn on top of Orbits'''
    def __init__(self, tno, number):
        self.tno = tno # the TNO to which this Dot belongs
        self.number = number # the index in list of this TNO's Dots
        orbitid = tnoinfo[self.tno]['ORBITID']
        ra = tnodetect[(tnodetect['ORBITID']== orbitid)][self.number]['RA']
        dec = tnodetect[(tnodetect['ORBITID']== orbitid)][self.number]['DEC']
        self.coords = (ra,dec)
        self.tdb0 = 2000 + tnodetect[(tnodetect['ORBITID']== orbitid)][self.number]['TDB'] # Julian time of corresponding detection
        self.lastdrawn = None

    def plot(self, col):
        '''basic plot'''
        pl.plot(self.coords[0], self.coords[1], 'o', color=str(col))
        pl.gca().invert_xaxis()
        
    def wcsplot(self, wcs, ax, col, marker, alph):
        '''plot using pixel coordinates given a wcs'''
        rapix, decpix = wcs.all_world2pix(self.coords[0], self.coords[1], 0)
        pixcoords = (rapix, decpix)

    def draw(self, tdb, col):
        '''draws itself at the right time in an animation'''
        # tdb: the time of the frame i'm drawing right now (comes from a 'clock' value in the orbit class)
        age = (tdb - self.tdb0)
        # increasing either of these makes fade and shrink effect happen more quickly
        fadetime = 4
        shrinktime = 45
        if age < 0: # detection has not occured yet
            return
        if age >= 0: # detection has occured
            # first, remove any previous drawings of this Dot (so it doesn't draw on top of itself)
            if self.lastdrawn != None:
                self.lastdrawn[0].set_xdata([])
                self.lastdrawn[0].set_ydata([])
            alpha = max(0.3, (1 - (age * fadetime)))
            marker = max(5, (9 - (age * shrinktime)))
            #plot the Dot with the appropriate size and opacity according to its age - making a nice "Pop!"
            self.lastdrawn = pl.plot(self.coords[0], self.coords[1], 'o', color=str(col), alpha=alpha, markersize=marker)

    def wcsdraw(self, wcs, ax, tdb, col):
        '''same as draw, but for pixel coordinates / wcs'''
        # tdb: the time of the frame i'm drawing right now
        age = (tdb - self.tdb0)
        # increasing either of these makes fade and shrink effect happen more quickly
        fadetime = 4
        shrinktime = 45
        if age < 0: #detection has not occured yet
            return
        if age >= 0: #detection has occured
            #first, remove any previous drawings of this Dot (so it doesn't draw on top of itself)
            if self.lastdrawn != None:
                self.lastdrawn[0].set_xdata([])
                self.lastdrawn[0].set_ydata([])
            alpha = max(0.3, (1 - (age * fadetime)))
            marker = max(6, (10 - (age * shrinktime)))
            # plot the Dot with the appropriate size and opacity according to its age
            rapix, decpix = wcs.all_world2pix(self.coords[0], self.coords[1], 0)
            pixcoords = (rapix, decpix)
            self.lastdrawn = ax.plot(pixcoords[0], pixcoords[1], 'o', color=str(col), alpha=alpha, markersize=marker)

            
            
class Orbit:
    '''the TNO path which can draw itself'''
    def __init__(self, ID):
        self.ID = ID # ID = row index, different from 'ORBITID' in table
        a = tnoinfo[self.ID]['a']
        e = tnoinfo[self.ID]['e']
        i = tnoinfo[self.ID]['i']
        O = tnoinfo[self.ID]['Omega']
        o = tnoinfo[self.ID]['omega']
        self.elements = [a,e,i,O,o]
        self.T_p = tnoinfo[self.ID]['T_p']
        self.orbitid = tnoinfo[self.ID]['ORBITID']
        self.total = tnodetect[(tnodetect['ORBITID']== self.orbitid)]['RA'].size # total # of detections
        self.dots = [Dot(self.ID, k) for k in range(self.total-1)]
        self.firsttime = self.dots[0].tdb0
        self.distance = tnoinfo[self.ID]['d']
        # for color coding based on distance
        max_d = tnoinfo['d'][(tnoinfo['d'] > 92)][0]
        min_d = tnoinfo['d'][(tnoinfo['d'] < 28.7)][0]
        self.color_val = (self.distance - min_d)/(55 - min_d) # determines the color of this Orbit given a colormap
        
    def plotDots(self, **argv): #argv can be color, ex: orb.plotDots(col='k')
        '''plots all detection dots of the TNO, used for plots (not animations)'''
        for d in self.dots:
            d.plot(**argv)
            
    def wcsPlotDots(self, wcs, ax, col, marker, alph):
        '''same as plotDots, but for pixel coordinates / wcs'''
        for d in self.dots:
            d.wcsplot(wcs, ax, col, marker, alph)
            
    def drawDots(self, tdb, col):
        '''time-dependent, flashes detection dots, used in animations'''
        for d in self.dots:
            d.draw(tdb,col)
            
    def wcsDrawDots(self, wcs, ax, tdb, col):
        '''time-dependent, flashes detection dots, used in animations (pixel / wcs)'''
        for d in self.dots:
            d.wcsdraw(wcs,ax,tdb,col)
            
    def wcsDrawOneDot(self, wcs, ax, tdb, col, num):
        '''for cherry-picking dot drawing
           for instance, if we want to start drawing from the 2nd or 3rd dot instead of first'''
        self.dots[num].wcsdraw(wcs, ax, tdb, col)
    
    def getPath(self, tdbstart, tdbspan):
        '''get (ra,dec) coordinates of TNO path from its orbital elements (the meat and potatoes!!!)'''
        # ** biggest issue: having to change this value often if animating on several different scales ** #
        # interval determines how many coordinate points to grab when drawing over a span of time
        #   - for up-close orbit drawing, interval should be larger 
        #     (~730*tdbspan = update TNO position ~2 times a day)
        #   - for far away orbit drawing, like all-survey animations, interval should be smaller. 
        #     (122*tdbspan = update position ~once every 3 days)
        #   - for Observer drawing, where the camera leaves earth and zooms out, interval should be even smaller!
        # ** this should be an argument in a function instead so it's easier to keep track of ** #
        interval = 300 * tdbspan
        t_ps, astrotimes, clocktimes, displaytimes, earth_pos, tno_pos, clocktimesval = [], [], [], [], [], [], []
        varyelement = np.zeros([int(interval)+1,6])
        
        for j in range(0,int(interval)):
            # generate several timelines:
            # - t_ps: varying time of perihelion, only used for Pedro's transformation code
            # - astrotimes: astropy-formatted times for retrieving earth's location
            # - clocktimes: the timeline used for final animation and Dot drawing (ex '2016.0001')
            t_ps.append(((self.T_p + (2016-tdbstart))-((tdbspan/(interval))*j)))
            astrotimes.append(astropy.time.Time(((2016.0001+(2016-tdbstart))+((tdbspan/(interval)*j))), format='decimalyear', scale='tdb'))
            clocktimes.append(astropy.time.Time(((2016.0001-(2016-tdbstart))+((tdbspan/(interval)*j))), format='decimalyear', scale='tdb'))
            displaytimes.append(clocktimes[j].utc.iso[0:10]) # clocktimes but in calendar format (ex '2016-01-01')
            
            # get earth's cartesian coordinates
            earths = get_body_barycentric('earth', Time(astrotimes[j]))
            earth_pos.append(np.array((earths.x.value, earths.y.value, earths.z.value)))
            # get tno's cartesian coordinates
            b = np.insert(self.elements, 5, t_ps[j])
            varyelement[j] = b        
        tno_pos = tr.keplerian_to_cartesian(varyelement[0:int(interval)], epoch = 2016)[:,0:3] # also given to Observer class
        
        # subtract earth and TNO coordinates to get vector coordinates (i.e. corrected TNO coordinates)
        vectors = tno_pos - earth_pos
        
        # convert from cartesian to ra and dec
        all_radec = SkyCoord(x=vectors[:,0], y=vectors[:,1], z=vectors[:,2], unit='AU', representation_type='cartesian')
        all_radec.representation_type = 'spherical'
        all_radec = np.array([all_radec.ra.value, all_radec.dec.value])
        # correct ra outliers (ex: some are recorded as -20 instead of 340)
        ras = np.where(all_radec[0] < 250, all_radec[0], all_radec[0] - 360)
        decs = all_radec[1]
        
        # optional modification: only retrieve and visualize data after the TNOs first detection
        if self.firsttime < (tdbstart+tdbspan):
            indexstart = list(map(lambda i: i> self.firsttime, clocktimesval)).index(True) 
            seenpath = [ras[indexstart:], decs[indexstart:]]
        else:
            seenpath = [0,0]
        
        return ras, decs, clocktimes, displaytimes, tno_pos, vectors, seenpath, indexstart
    
    def wcsPath(self, wcs, start, span):
        '''converts getPath values into pixel coordinates given a wcs'''
        info = self.getPath(start,span)
        seenpath = info[6]
        rapix, decpix = wcs.all_world2pix(info[0], info[1], 0)
        pixcoords = (rapix, decpix)
        seenrapix, seendecpix = wcs.all_world2pix(seenpath[0], seenpath[1], 0)
        seenpixcoords = (seenrapix, seendecpix)
        return pixcoords, seenpixcoords # choose when animating: wait to draw until after first detection or not
    
    def extremes(self, tdbstart, tdbspan):
        '''over the span given, this returns the extreme values of RA and DEC (axis limit purposes)'''
        values = self.getPath(tdbstart,tdbspan)
        minra = min(values[0])
        mindec = min(values[1])
        maxra = max(values[0])
        maxdec = max(values[1])
        return minra, mindec, maxra, maxdec
    
    def plot(self, tdbstart, tdbspan, col):
        '''plots an entire orbit (image, not animation)'''
        rapoints = self.getPath(tdbstart, tdbspan)[0]
        decpoints = self.getPath(tdbstart, tdbspan)[1]
        pl.plot(rapoints, decpoints, '-', color=col, alpha=1, markersize=4)
        # invert the ra values, for astronomer purposes
        pl.gca().invert_xaxis()
        
    def wcsPlot(self, axis, wcs, tdbstart, tdbspan, col):
        '''plots an entire orbit (image, not animation; pixels / wcs)'''
        rapoints = self.wcsPath(wcs,tdbstart,tdbspan)[0][0]
        decpoints = self.wcsPath(wcs,tdbstart,tdbspan)[0][1]
        axis.plot(rapoints, decpoints, '-', color=col, alpha=1, lw=1)
    
    def draw(self, points, xdata, ydata, tdbstart, tdbspan, frame):
        ''' appends point data to a pl.plot function defined in animation function
            - ** very slow ... room for improvement!
            - ** I've been using an animation template for drawing Orbits,
            - ** it would be good to make self-drawing Orbits (like how the Dots are self-drawing)'''
        x = info[0][frame]
        y = info[1][frame]
        xdata.append(x)  
        ydata.append(y)
        points.set_data(xdata, ydata)
        return points

    
    
def getPlanetPos(planets, start, span):
    '''function which grabs planet positions from astropy'''
    #coordinate intervals with that in the Observer (+ Orbit) ** room for improvement!! **
    inter = 300
    interval = span*inter  
    #prepare the times and get positions at those times
    astrotimes, planets_pos = [], []
    for j in range(0,int(interval)):
        astrotimes.append(astropy.time.Time(((2016.0001-(2016-start))+
                ((span/(interval)*j))), format='decimalyear', scale='tdb'))        
    for i in range(len(planets)):
        get_planet = get_body_barycentric(planets[i], Time(astrotimes))
        planet_pos = np.array([get_planet.x.value, 
                               get_planet.y.value, 
                               get_planet.z.value]).T
        planets_pos.append(planet_pos)    
    return planets_pos