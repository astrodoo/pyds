"""
filename: 
     readcol.py
 
PURPOSE:
     read the columns with formats from ascii data

     ### Note that in current verion of python (2.7), nan cannot 
     be managed in integer type (only work for float format). 
     Instead, it will produce weird number (i.e. -91212141412412)

Routines:
     readcol: read the columns

Written by:
     Doosoo Yoon
     Shanghai Astronomical Observatory
   
History:
     Written, 23 November 2017
"""
import numpy as np
import sys

def readcol(file,format=[],nskip=0):
    """
    read the columns with formats from ascii data
    args:
        file: file name
    keywords:
        format: array of formats. The choices of the format are 'i','f','s','x',
               coressponding to np.int32, np.float64, string (S25), and skipping, respectively.
        nskip: the number of lines to be skipped from the top 
    return:
        (vars)  -> number of vars are same as the number of 'i','f','s' in fortmats 
    """

    # read the file
    target = open(file,'r')
    text = target.readlines()
    target.close()
    
    # skip the lines
    text = text[nskip:]

    nlines = np.str(len(text))
    ncols = len(text[0].split())

    nformat = len(format)

    # number of variables to be read
    nvars = 0
    for iformat in format:
        if ((iformat=='i') | (iformat=='f') | (iformat=='s')): nvars += 1

    types = ["" for i in range(nvars)]
    cols = np.asarray([False for i in range(ncols)])

    if (nformat > ncols):
        sys.exit('number of formats should not be larger than that of columns')

    # assign the type of the variable
    j = 0
    for i,iformat in enumerate(format):
        if (iformat == 'i'):
            types[j] = 'np.int32'
            cols[i] = True
        elif (iformat == 'f'):
            types[j] = 'np.float64'
            cols[i] = True
        elif (iformat == 's'):
            types[j] = "'S25'"    # for longer string array
            cols[i] = True
        elif (iformat =='x'):
            continue
        else:
            sys.exit('the allowed formats are i/f/s/x for integer,float,string and skip, respectively.')

        j += 1

    # col1,col2,col3, ...
    colarr = ['col'+str(i) for i in range(nvars)]
    
    # initialize the arrays
    for i,icolarr in enumerate(colarr):
        exec(icolarr+'=np.zeros('+nlines+',dtype='+types[i]+')')

    for i,iline in enumerate(text):
        colstrs_tot = np.asarray((iline.split()))
        colstrs = np.asarray(colstrs_tot[cols])

        for j,jcolarr in enumerate(colarr):
            
            if (types[j] == 'np.float64'):
                exec(jcolarr+"[i]=np.float64(colstrs[j].lower().replace('d','e'))")
            elif (types[j] == 'np.int32'):
    # In current version of python(2.7), nan cannot be managed in integer type (only work for float format). 
    # Instead, it will produce weird number (i.e. -91212141412412)
                exec(jcolarr+"[i]=np.int32(np.float64(colstrs[j].lower().replace('d','e')))")
            else:
                exec(jcolarr+'[i]=np.str(colstrs[j])')
                
    colarrstr = ', '.join(colarr)
    return eval(colarrstr)
