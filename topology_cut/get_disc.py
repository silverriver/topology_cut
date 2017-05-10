# Create a circle from a center point and a circumference.
import rhinoscriptsyntax as rs
import math
import Rhino

def GetFracture ():
    dx = rs.GetReal("Enter x-axis direction:")
    if not dx:
        dx = 0.0
    print "x-axis direction: ", dx
    objs = rs.ObjectsByType(8,False)
    cir_objs = rs.ObjectsByType(4)
    
    #delete all circles
    for cur_obj in cir_objs:
        if rs.IsCircle(cur_obj):
            rs.DeleteObject(cur_obj)
    occor = []
    radius = []
    center = []
    
    #for all surface
    for obj in objs:
        rs.UnselectAllObjects()
        rs.SelectObject(obj)
        rs.Command ("_Silhouette")
        created_objs = rs.LastCreatedObjects()
        #not a circle
        if len(created_objs) !=1:
            rs.DeleteObjects(created_objs)
            print "unvailded surface"
            continue
        created_objs = created_objs[0]
        #not a circle
        if not rs.IsCircle(created_objs):
            rs.DeleteObject(created_objs)
            print "unvailded surface, not circle"
            continue
        point = rs.CircleCenterPoint(created_objs)
        center.append(point)
        r = rs.CircleRadius(created_objs)
        radius.append(r)
        normal = rs.SurfaceNormal(obj,[0,0])
        occor.append(GetDirDip(normal,dx))
        rs.DeleteObject(created_objs)
    print center
    print occor
    print radius
    path = rs.DocumentPath()
    path_l = path.split("\\")
    path_l[len(path_l)-1] = "fracture.dat"
    file = open("\\".join(path_l),"w")
    file.write(str(len(occor)))
    file.write('\n')
    for i in range (len(occor)):
        file.write("%.15f " % center[i][0])
        file.write("%.15f " % center[i][1])
        file.write("%.15f " % center[i][2])
        file.write ("%.15f " % occor[i][0])
        file.write ("%.15f " % occor[i][1])
        file.write ("%.15f " % radius[i])
        file.write ("0.0001 0.1 30\n")
    file.close ()
def GetDirDip (normal,dx):
    if normal[2]<0:
        normal[0] = -normal[0]
        normal[1] = -normal[1]
        normal[2] = -normal[2]
    x = normal[0]
    y = normal[1]
    z = normal[2]
    d = x*x+y*y
    if d==0:
        return [0,0]   #dip-direction==0 ; dip-angle==0
    d = math.sqrt(d)
    theta = math.acos (x/d)
    theta = theta*180/math.pi
    if (y<0):
        theta = 360-theta
    dip_dir = dx-theta;
    if (dip_dir<0):
        dip_dir +=360
    alpha = z/math.sqrt (x*x+y*y+z*z)
    alpha = math.asin (alpha)
    alpha = alpha*180/math.pi
    dip = 90-alpha
    return [dip_dir,dip]
    
    
GetFracture()
#rs.Command("_Silhouette")