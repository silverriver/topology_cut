import rhinoscriptsyntax as rs
import math
import Rhino

def GetNormal(dip_dir, dip,dx):
    dip_dir = math.radians(dip_dir)
    dip = math.radians(dip)
    dx = math.radians(dx)
    return [math.sin(dip)*math.cos(dip_dir-dx),-math.sin(dip)*math.sin(dip_dir-dx),math.cos(dip)]
def DrawDisc (filename,dx):
    result = []
    try:
        fp = open(filename)
    except:
        print "can't open "+filename
        return
    curve_group = rs.AddGroup("curve")
    surf_group = rs.AddGroup("surface")
    for line in fp:
        result.append(map(float,[x for x in line.split(' ') if len(x)!=0]))
    n = int(result[0][0])
    for i in range(n):
        i += 1
        if len(result[i])!=9:
            break
        print [result[i][3],result[i][4],dx]
        print GetNormal(result[i][3],result[i][4],dx)
        id = rs.AddCircle(rs.PlaneFromNormal(result[i][0:3],GetNormal(result[i][3],result[i][4],dx)),result[i][5])
        rs.AddObjectToGroup(id,curve_group)
        rs.AddObjectToGroup(rs.AddPlanarSrf(id),surf_group)
    
dx = rs.GetReal("Enter x-axis direction:")
if not dx:
    dx = 0.0
DrawDisc ("contri_fracture_xyzabr.dat",dx)