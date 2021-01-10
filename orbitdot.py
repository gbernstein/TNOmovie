# Contents: Orbit class, Dot class, Timeline class
#     - Dot: the recorded detections of a TNO
#     - Orbit: the path of a given TNO
#     - Timeline: an add-on to an animation, fills up over time

import transformation as tr
import matplotlib.pyplot as pl
from astropy.io import fits
import astropy.table as tb
from astropy import wcs
from astropy.wcs import WCS
import astropy
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_body_barycentric_posvel, get_body, get_body_barycentric, SkyCoord
import numpy as np



#relevant FITS files
alltnoi = tb.Table.read('dataa/y4_tno_table.fits')
big_tno = tb.Table.read('dataa/orbits_culled_ast.fits', 2)

#include calculated RA and DEC columns
big_tno['RA'] = big_tno['DETECT'][:,0]*180/np.pi
big_tno['DEC'] = big_tno['DETECT'][:,1]*180/np.pi
alltnod = big_tno[(big_tno['OBJECTID'] > 0)]

#sort them into 2 main lists,
# - tnoinfo: list with the orbital elements of each TNO
# - tnodetect: list with the times and positions of each detection of each TNO
tnoinfo = []
for i in range(0,316):
    if -1 < i and i < 10:
        tnoinfo.append(alltnoi[(alltnoi['DES']=="DES000%s" % (i))])
    if i == 10:
        tnoinfo.append(alltnoi[(alltnoi['DES']=="DES0010")])
    if 10 < i and i < 100:
        tnoinfo.append(alltnoi[(alltnoi['DES']=="DES00%s" % (i))])
    if i == 100:
        tnoinfo.append(alltnoi[(alltnoi['DES']=="DES0100")])
    if 100< i and i < 1000:
        tnoinfo.append(alltnoi[(alltnoi['DES']=="DES0%s" % (i))])
tnodetect = []
for i in range(0,316):
    tnodetect.append(alltnod[(alltnod['ORBITID']== i)])

class Dot:
    '''the detection flashes drawn on top of orbits'''
    def __init__(self, tno, number):
        self.tno = tno
        self.number = number #the index in list of a single TNO's dots
        ra = tnodetect[self.tno]['RA'][self.number]
        dec = tnodetect[self.tno]['DEC'][self.number]
        self.coords = (ra,dec)
        self.tdb0 = 2000 + tnodetect[self.tno]['TDB'][self.number] #Julian time of corresponding detection
        self.lastdrawn = None

    def plot(self, col):
        '''basic plot'''
        pl.plot(self.coords[0], self.coords[1], 'o', color=str(col))
        pl.gca().invert_xaxis()
        
    def wcsplot(self, wcs, col):
        '''plot on top of an image, given its wcs'''
        rapix, decpix = wcs.all_world2pix(self.coords[0], self.coords[1], 0)
        pixcoords = (rapix, decpix)
        pl.plot(pixcoords[0], pixcoords[1], 'o', color=str(col), markersize=1.5)
        
    def draw(self, tdb, col):
        '''draws itself at the right time in an animation'''
        #tdb: the time of the frame i'm drawing right now (comes from a 'clock' value in the orbit class)
        age = (tdb - self.tdb0)
        #increasing either of these makes fade and shrink effect happen more quickly
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
            marker = max(5, (9 - (age * shrinktime)))
            #plot the Dot with the appropriate size and opacity according to its age - making a nice "Pop!"
            self.lastdrawn = pl.plot(self.coords[0], self.coords[1], 'o', color=str(col), alpha=alpha, markersize=marker)

    def wcsdraw(self, wcs, ax, tdb, col):
        '''same as draw, but for drawing atop an image'''
        #tdb: the time of the frame i'm drawing right now
        age = (tdb - self.tdb0)
        #increasing either of these makes fade and shrink effect happen more quickly
        fadetime = 10 #was 4
        shrinktime = 55 #was 45
        if age < 0: #detection has not occured yet
            return
        if age >= 0: #detection has occured
            #first, remove any previous drawings of this Dot (so it doesn't draw on top of itself)
            if self.lastdrawn != None:
                self.lastdrawn[0].set_xdata([])
                self.lastdrawn[0].set_ydata([])
            alpha = max(0.2, (1 - (age * fadetime)))
            marker = max(1, (4 - (age * fadetime)))
            #plot the Dot with the appropriate size and opacity according to its age
            rapix, decpix = wcs.all_world2pix(self.coords[0], self.coords[1], 0)
            pixcoords = (rapix, decpix)
            self.lastdrawn = ax.plot(pixcoords[0], pixcoords[1], 'o', color=str(col), alpha=alpha, markersize=marker)

            
            
class Orbit:
    '''the TNO path which can draw itself'''
    def __init__(self, ID):
        self.ID = ID
        a = tnoinfo[self.ID]['a'][0]
        e = tnoinfo[self.ID]['e'][0]
        i = tnoinfo[self.ID]['i'][0]
        O = tnoinfo[self.ID]['Omega'][0]
        o = tnoinfo[self.ID]['omega'][0]
        self.elements = [a,e,i,O,o]
        self.T_p = tnoinfo[self.ID]['T_p']
        self.total = tnodetect[self.ID]['RA'].size
        self.dots = [Dot(self.ID, k) for k in range(self.total-1)]
        self.firsttime = self.dots[0].tdb0
        self.distance = alltnoi[self.ID]['d']
        #TNO 237 has the farthest distance (~90 AU) for reference
        max_d = alltnoi['d'][(alltnoi['d'] > 92)][0]
        #TNO 122 has the smallest distance (~35 AU) for reference
        min_d = alltnoi['d'][(alltnoi['d'] < 28)][0]
        self.color_val = (self.distance - min_d)/(55 - min_d) #determines the color of the drawn orbit, given a colormap
        
    def plotDots(self,**argv): #argv can be color, ex: orb.plotDots(col='k')
        '''plots all detection dots of the TNO, used for plots (pictures)'''
        for d in self.dots:
            d.plot(**argv)
            
    def wcsPlotDots(self, wcs, **argv):
        for d in self.dots:
            d.wcsplot(wcs, **argv)
            
    def drawDots(self,tdb,col): #ex: orb.drawDots(tdb, col='k')
        '''time-dependent, flashes detection dots, used in animations'''
        for d in self.dots:
            d.draw(tdb,col)
            
    def wcsDrawDots(self,wcs,ax,tdb,col): #ex: orb.drawDots(tdb, col='k')
        '''time-dependent, flashes detection dots, used in animations'''
        for d in self.dots:
            d.wcsdraw(wcs,ax,tdb,col)
    
    def getPath(self, tdbstart, tdbspan): #tdbspan = year span (integer)
        '''get (ra,dec) coordinates of the TNO path from the orbital elements (the meat and potatoes of my code)'''
        
        interval = 122 * tdbspan ## !!!!!! 
        #for up-close orbit drawing, interval should be larger (~730*tdbspan = update TNO position ~2 times a day)
        #for far away orbit drawing, like all-survey animations, interval should be smaller. (122*tdbspan = update position ~once every 3 days)
        t_ps, astrotimes, clocktimes, displaytimes, earth_pos, tno_pos, clocktimesval = [], [], [], [], [], [], []
        varyelement = np.zeros([interval+1,6])
        
        for j in range(0,interval):
            #generate several timelines:
            # - t_ps: varying time of perihelion, only used for Pedro's transformation code
            # - astrotimes: astropy-formatted times for retrieving earth's location
            # - clocktimes: the timeline used for final animation and Dot drawing (ex '2016.0001')
            t_ps.append(((self.T_p[0] + (2016-tdbstart))-((tdbspan/(interval))*j)))
            astrotimes.append(astropy.time.Time(((2016.0001+(2016-tdbstart))+((tdbspan/(interval)*j))), format='decimalyear', scale='tdb'))
            clocktimes.append(astropy.time.Time(((2016.0001-(2016-tdbstart))+((tdbspan/(interval)*j))), format='decimalyear', scale='tdb'))
            clocktimesval.append((astropy.time.Time(((2016.0001-(2016-tdbstart))+((tdbspan/(interval)*j))), format='decimalyear', scale='tdb')).value) #is this necessary?
            displaytimes.append(clocktimes[j].utc.iso[0:10]) #clocktimes but in calendar format (ex '2016-01-01')
            
            #get earth's cartesian coordinates
            earths = get_body_barycentric('earth', Time(astrotimes[j]))
            earth_pos.append(np.array((earths.x.value, earths.y.value, earths.z.value)))
            
            #get tno's cartesian coordinates
            b = np.insert(self.elements, 5, t_ps[j])
            varyelement[j] = b        
        tno_pos = tr.keplerian_to_cartesian(varyelement[0:interval], epoch = 2016)[:,0:3] 
    
        #subtract earth and TNO coordinates to get vector coordinates (i.e. corrected TNO coordinates)
        vectors = tno_pos - earth_pos
        
        #convert from cartesian to ra and dec
        all_radec = SkyCoord(x=vectors[:,0], y=vectors[:,1], z=vectors[:,2], unit='AU', representation_type='cartesian')
        all_radec.representation_type = 'spherical'
        all_radec = np.array([all_radec.ra.value, all_radec.dec.value])
        ras = np.where(all_radec[0] < 250, all_radec[0], all_radec[0] - 360) #correcting outlier values
        decs = all_radec[1]
        
        #optional modification of list of ras and decs: 
        #     only plot/draw points after the first detection
        if self.firsttime < (tdbstart+tdbspan):
            indexstart = list(map(lambda i: i> self.firsttime, clocktimesval)).index(True) 
            seenpath = [ras[indexstart:], decs[indexstart:]]
        
        return ras, decs, clocktimes, displaytimes, tno_pos, clocktimesval, seenpath, indexstart #tno_pos is given to Observer class
    
    def wcsPath(self, wcs, start, span):
        '''converts getPath values into pixel coordinates, given a wcs'''
        info = self.getPath(start,span)
        seenpath = info[6]
        rapix, decpix = wcs.all_world2pix(info[0], info[1], 0)
        pixcoords = (rapix, decpix)
        seenrapix, seendecpix = wcs.all_world2pix(seenpath[0], seenpath[1], 0)
        seenpixcoords = (seenrapix, seendecpix)
        return pixcoords, seenpixcoords
    
    def extremes(self, tdbstart, tdbspan):
        '''over the span given, this returns the extreme values of RA and DEC (axis limit purposes)'''
        values = self.getPath(tdbstart,tdbspan)
        minra = min(values[0])
        mindec = min(values[1])
        maxra = max(values[0])
        maxdec = max(values[1])
        return minra, mindec, maxra, maxdec
    
    def plot(self, tdbstart, tdbspan, col):
        '''plots an entire orbit as a picture'''
        rapoints = self.getPath(tdbstart, tdbspan)[0]
        decpoints = self.getPath(tdbstart, tdbspan)[1]
        pl.plot(rapoints, decpoints, '-', color=str(col), alpha=1, markersize=4)
        #invert the ra values, for astronomer purposes
        pl.gca().invert_xaxis()
        
    def wcsPlot(self, axis, wcs, tdbstart, tdbspan, colormap): ##added AXIS
        rapoints = self.wcsPath(wcs,tdbstart,tdbspan)[0][0]
        decpoints = self.wcsPath(wcs,tdbstart,tdbspan)[0][1]
        axis.plot(rapoints, decpoints, '-', color=colormap(self.color_val), alpha=1, markersize=4)
    
    def plotSeen(self, tdbstart, tdbspan, col):
        '''only plot points which come after the TNO's first detection'''
        info = self.getPath(tdbstart,tdbspan)
        seenpath = info[6]
        pl.plot(seenpath[0], seenpath[1], '-', color=str(col), alpha=1, markersize=4)
        pl.gca().invert_xaxis()
    
    def draw(self, points, xdata, ydata, tdbstart, tdbspan, frame):
        ''' appends point data to a pl.plot function defined in animation function
            - room for improvement, very slow when info is calculated in this function as opposed to in the animation cell'''
        #info = self.getPath(tdbstart, tdbspan)
        x = info[0][frame]
        y = info[1][frame]
        xdata.append(x)  
        ydata.append(y)
        points.set_data(xdata, ydata)
        return points
    

    
class Timeline:
    def __init__(self,start,end):
        self.start = start #values in years (integer)
        self.end = end

    def plot(self, main_axes, length_start, height_start):
        '''situate the timeline on a larger plot or image
             - length_start: decimal value, representing the length of the larger axes at which the timeline begins
             - height_start: decimal value, representing the height of the larger axes at which the timeline begins'''
        with pl.rc_context({'axes.edgecolor':'white', 'xtick.color':'white'}):
        # Temporary rc parameters in effect
            ax1 = pl.axes([length_start,height_start,0.4,0.01], yticklabels=[], xlim=[self.start, self.end])
            ax1.tick_params(axis = "y", which = "both", bottom = False, top = False, left=False)
            ax1.spines['top'].set_visible(False)
            ax1.patch.set_facecolor('none')
            pl.rcParams['axes.xmargin'] = 0
            pl.rcParams['axes.ymargin'] = 0
            #these parameters give turn this timeline subplot into more of a spine
            #note: color is automatically white, I am using this on top of a black sky image

    def time(self,frame,totalframe):
        '''incrementally adds points to the timeline to fill it up'''
        inter = (self.end - self.start) / totalframe #the interval which seperates each point on the line
        
        #to fill up smoothly, there should be > 37 total frames in the animation
        if totalframe < 37:
            marker = 10
        else:
            marker = 6
        pl.plot(self.start+(frame*inter), 0.5, color='white', markersize=marker, marker = "s")