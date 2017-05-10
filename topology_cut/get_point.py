import rhinoscriptsyntax as rs
import math
import Rhino

def GetPoint ():
    pl = rs.ObjectsByType(1,False,1)
    for p in pl:
        print rs.PointCoordinates(p)
    
GetPoint()
    