'''
Created on May 18, 2014
This one reads dpi
@author: Keisuke
'''
from pylab import * 

def sf(figName):
    First_part="C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\02\\19 Rhodopsin bacteria\\"
    Figure_name=str(figName)
    Last_part=".png"
    savefig(First_part+Figure_name+Last_part,dpi=450)
