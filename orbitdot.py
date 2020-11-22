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
        self.number = number
        ra = tnodetect[self.tno]['RA'][self.number]
        dec = tnodetect[self.tno]['DEC'][self.number]
        self.coords = (ra,dec)
        self.tdb0 = 2000 + tnodetect[self.tno]['TDB'][self.number] #Julian time of corresponding detection
        self.lastdrawn = None

    def plot(self, col):
        '''basic plot picture'''
        pl.plot(self.coords[0], self.coords[1], 'o', color=str(col))
        
    def draw(self, tdb, col):
        '''draws itself at the right time in an animation'''
        #tdb: the time of the frame i'm drawing right now
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
            #plot the Dot with the appropriate size and opacity according to its age
            self.lastdrawn = pl.plot(self.coords[0], self.coords[1], 'o', color=str(col), alpha=alpha, markersize=marker)

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
        
    def plotDots(self,**argv): #argv can be color, ex: orb.plotDots(col='k')
        '''plots all detection dots of the TNO, used for plots (pictures)'''
        for d in self.dots:
            d.plot(**argv)
    
    def drawDots(self,tdb,col): #ex: orb.drawDots(tdb, col='k')
        '''time-dependent, flashes detection dots, used in animations'''
        for d in self.dots:
            d.draw(tdb,col)
    
    def getPath(self, tdbstart, tdbspan): #tdbspan = year span (integer)
        '''get (ra,dec) coordinates of the TNO path (the meat and potatoes)'''
        interval = 730 * tdbspan #retrieve and update the tno's position ~2 times a day
        t_ps, astrotimes, clocktimes, displaytimes, earth_pos, tno_pos = [], [], [], [], [], []
        varyelement = np.zeros([interval+1,6])
        
        for j in range(0,interval):
            #generate several timelines:
            # - t_ps: varying time of perihelion, only used for Pedro's transformation code
            # - astrotimes: astropy-formatted times for retrieving earth's location
            # - clocktimes: the timeline used for final animation and Dot drawing (ex '2016.0001')
            t_ps.append(((self.T_p[0] + (2016-tdbstart))-((tdbspan/(interval))*j)))
            astrotimes.append(astropy.time.Time(((2016.0001+(2016-tdbstart))+((tdbspan/(interval)*j))), format='decimalyear', scale='tdb'))
            clocktimes.append(astropy.time.Time(((2016.0001-(2016-tdbstart))+((tdbspan/(interval)*j))), format='decimalyear', scale='tdb'))
            displaytimes.append(clocktimes[j].utc.iso[0:10]) #clocktimes but in calendar format (ex '2016-01-01')
            
            #get earth's cartesian coordinates
            earths = get_body_barycentric('earth', Time(astrotimes[j]))
            earth_pos.append(np.array((earths.x.value, earths.y.value, earths.z.value)))
            #numpy faster than append, earth_pos = np.array(earth_pos)
            
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
            
        return ras, decs, clocktimes, displaytimes, tno_pos, vectors #tno_pos is given to observer
    
    def wcsPath(self, wcs, start, span):
        '''converts getPath values into pixel coordinates, given a wcs'''
        rapix, decpix = wcs.all_world2pix(self.getPath(start,span)[0], self.getPath(start,span)[1],0)
        pixcoords = (rapix, decpix)
        return pixcoords
    
    def extremes(self, tdbstart, tdbspan):
        '''over the span given, this returns the extreme values of RA and DEC (axis limit purposes)'''
        values = self.getPath(tdbstart,tdbspan)
        minra = min(values[0])
        mindec = min(values[1])
        maxra = max(values[0])
        maxdec = max(values[1])
        return minra, mindec, maxra, maxdec
    
    def plot(self, tdbstart, tdbspan, col):
        '''plots an entire orbit, as a picture'''
        rapoints = self.getPath(tdbstart, tdbspan)[0]
        decpoints = self.getPath(tdbstart, tdbspan)[1]
        pl.plot(rapoints, decpoints, '-', color=str(col), alpha=1, markersize=4)
    
    def draw(self, points, xdata, ydata, tdbstart, tdbspan, frame):
        ''' appends point data to a pl.plot function defined in animation function'''
        info = self.getPath(tdbstart, tdbspan)
        x = info[0][frame]
        y = info[1][frame]
        xdata.append(x)  
        ydata.append(y)
        points.set_data(xdata, ydata)
        return points