#     - Timeline: add-on for animations, fills up with every frame

# import matplotlib
import matplotlib.pyplot as pl

# *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** # *** #

class Timeline:
    '''self-filling timeline'''
    def __init__(self,start,end):
        # start and end dates in years (integer)
        self.start = start
        self.end = end 

    def plot(self, main_axes, length_start, height_start):
        '''situate the timeline on a larger plot or image
             - length_start: decimal value, representing the length of the larger axes at which the timeline begins
             - height_start: decimal value, representing the height of the larger axes at which the timeline begins'''
        with pl.rc_context({'axes.edgecolor':'white', 'xtick.color':'white'}):
        # temporary rc parameters in effect
            ax1 = pl.axes([length_start,height_start,0.4,0.01], yticklabels=[], xlim=[self.start, self.end])
            ax1.tick_params(axis = "y", which = "both", bottom = False, top = False, left=False)
            ax1.spines['top'].set_visible(False)
            ax1.patch.set_facecolor('none')
            pl.rcParams['axes.xmargin'] = 0
            pl.rcParams['axes.ymargin'] = 0
            # these parameters turn this timeline subplot into more of a spine
            # note: spine color is automatically white since I am using this on top of ~black sky images
    
    def time(self, frame, totalframe):
        '''incrementally adds points to the timeline to fill it up'''
        # the interval which seperates each point on the line; totalframe = how many frames you want it to take to fill
        inter = (self.end - self.start) / totalframe
        # to fill up smoothly, there should be > 37 total frames in the animation
        if totalframe < 37:
            marker = 10
        else:
            marker = 6
        pl.plot(self.start+(frame*inter), 0.5, color='white', markersize=marker, marker = "s")