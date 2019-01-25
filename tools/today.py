"""
filename: 
     today.py
 
PURPOSE:
     generate the string that indicates the today with the format of "yymmdd_"

Written by:
     Doosoo Yoon 
     University of Amsterdam
   
History:
     Written, 1 October 2018
"""
import datetime

class today:
    """ 
    generate the string that indicates the today with the format of "yymmdd_"

    ex) from pyds.tools import today
        print today.datestr
    """

    datestr = datetime.datetime.now().strftime("%y%m%d") + '_'
