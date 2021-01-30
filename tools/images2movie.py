"""
filename: 
     images2movie.py
 
PURPOSE:
     combine image files to generate movie file

Routines:
     images2movie: generate movie file

Written by:
     Doosoo Yoon
     University of Amsterdam, 
     Anton Pannekoek Institute of Astronomy
   
History:
     Written, 11 July 2019
"""
import numpy as np
import matplotlib.image as image
from matplotlib import pyplot as plt
from matplotlib import animation

def images2movie(filelist,dpi=100,interpolation='spline16',fps=24,movie_dpi=100, \
        out='movie_sample.mp4',silent=False,**kwargs):
    """
    combine image files to generate movie file
    args:
        filelist: list of file names
    keywords:
        dpi:           pixel/inch of the image (integer, default=100)
        interpolation: interpolation method in imshow (default='spline16', 
                        other options can be found in imshow document)
        fps:           frame per second (float, default=24)
        movie_dpi:     pixel/inch for the movie (integer, default=100)
        out:           file name of the output (default='movie_sample.mp4')
        silent:        if True, print out the process (boolean, default=False)
    """

    # if ffmpeg is not in the global environment, the full directory should be set.
    #plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\astrodoo\\ffmpeg-20190711-af9dc02-win64-static\\bin\\ffmpeg.exe'

    img = image.imread(filelist[0])
    height, width, depth = img.shape

    # What size does the figure need to be in inches to fit the image?
    figsize = width / float(dpi), height / float(dpi)
    
    # Create a figure of the right size with one axes that takes up the full figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0, 0, 1, 1])

    # Hide spines, ticks, etc.
    ax.axis('off')

    def init():
        res = ax.imshow(img,interpolation=interpolation,animated=True)
        return res,

    def animate(flist):
        img = image.imread(flist)       
        res = ax.imshow(img,interpolation=interpolation,animated=True)
        if (not silent):
            print('Insert %s to the movie (last snapshot:%s)'%(flist,filelist[-1]))
        return res,
    
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=filelist, blit=True, **kwargs)
    #anim = animation.FuncAnimation(fig, animate, frames=filelist, interval=25, blit=True, repeat=False)
    FFwriter = animation.FFMpegWriter(fps=fps)
    anim.save(out, writer = FFwriter, dpi=movie_dpi)

    print("movie file is generated: %s"%out)
    
    plt.close()
