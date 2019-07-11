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
        format: array of formats. The choices of the format are 'i','f','f64','s','x',
               coressponding to np.int32, np.float32, np.float64, string (S25), and skipping, respectively.
        nskip: the number of lines to be skipped from the top 
    return:
        (vars)  -> number of vars are same as the number of 'i','f','f64','s' in fortmats 
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
        if ((iformat=='i') | (iformat=='f') | (iformat=='f64') | (iformat=='s')): nvars += 1

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
            types[j] = 'np.float32'
            cols[i] = True
        elif (iformat == 'f64'):
            types[j] = 'np.float64'
            cols[i] = True
        elif (iformat == 's'):
            types[j] = "'S25'"    # for longer string array (no longer to be used)
            cols[i] = True
        elif (iformat =='x'):
            continue
        else:
            sys.exit('the allowed formats are i/f/f64/s/x for integer,float,string and skip, respectively.')

        j += 1
    
    # separate data elements 
    alldata = np.asarray([[data for data in line.split()] for line in text])
    selectdata = np.squeeze(alldata[:,cols])
        
    # col1,col2,col3, ...
    colarr = ['col'+str(i) for i in range(nvars)]
        
    for i,icolarr in enumerate(colarr):
        if (types[i]=='np.float32'):
            exec(icolarr+"= np.asarray([np.float32(num.lower().replace('d','e')) for num in selectdata[:,i]])")
        elif (types[i]=='np.float64'):
            exec(icolarr+"= np.asarray([np.float64(num.lower().replace('d','e')) for num in selectdata[:,i]])")
    # In current version of python(2.7), nan cannot be managed in integer type (only work for float format). 
    # Instead, it will produce weird number (i.e. -91212141412412)
        elif (types[i]=='np.int32'):
            exec(icolarr+"= np.asarray([np.int32(np.float64(num.lower().replace('d','e'))) \
                for num in selectdata[:,i]])")
        else:
            exec(icolarr+"=selectdata[:,i]")
                
    colarrstr = ', '.join(colarr)
    return eval(colarrstr)
