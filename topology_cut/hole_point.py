import rhinoscriptsyntax as rs
import math
import Rhino

def GetHolepoint ():
    id = rs.GetObject("select hole curve",4,True,True)
    n = rs.GetInteger("How may hole point")
    if not id or n <=0:
        print "illegal selection or interger"
        return 
    domain = rs.CurveDomain(id)
    t = 0
    while t<=1.0:
        point = rs.EvaluateCurve(id,domain[0]+t*(domain[1]-domain[0]))
        print point[0],' ',point[1],' ',point[2]
        t+=1.0/n
        
GetHolepoint()