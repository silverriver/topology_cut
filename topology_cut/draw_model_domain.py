import rhinoscriptsyntax as rs
import math
import Rhino

def CrossProduct (v1,v2):
    return [v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]]

def VPlus (v1,v2):
    return map (lambda x,y:x+y, v1,v2)
    
def VMinus (v1,v2):
    return map (lambda x,y:x-y, v1,v2)

def VMultNum (v,n):
    return map (lambda x: x*n,v)

def GetPlane (pl):
    center = [0.0,0.0,0.0]
    for p in pl:
        center = VPlus(center, p)
    center = VMultNum(center,1.0/len(pl))
    normal = [0,0,0]
    for i in range(len(pl)):
        i1 = i+1
        i0 = i-1
        if i1 >= len(pl):
            i1 = 0
        if i0<0:
            i0 = len(pl)-1
        normal = VPlus(normal, CrossProduct(VMinus(pl[i1],pl[i]),VMinus(pl[i0],pl[i])))
    return rs.PlaneFromNormal(center,normal)


def DrawDomain (filename):
    try:
        fp = open(filename)
    except:
        print "can't open "+filename
        return
    result = []
    for line in fp:
        line=line.strip()
        result += map(float,[x for x in line.split(' ') if len(x)!=0])
    pos = 0
    dx = result[pos]
    pos += 1
    
    slope_or = int(result[pos])
    pos += 1
    
    nf = int(result[pos])
    pos += 1
    
    np = int(result[pos])
    pos += 1
    
    ifnp = map(int,result[pos:pos+nf])
    pos += nf
    
    pos += nf
    
    face_type = map(int,result[pos:pos+nf])
    pos += nf
    
    facep = []
    for i in range(nf):
        facep.append (map(int,result[pos:pos+ifnp[i]]))
        pos += ifnp[i]
        
    xyzp = []
    for i in range(np):
        xyzp.append (result[pos:pos+3])
        pos += 3
    
    nsubdomain =int(result[pos])
    pos += 1
    
    pos += nsubdomain
    
    frami_nface = map(int,result[pos:pos+nsubdomain])
    pos+= nsubdomain
    
    frami_facej = []
    for i in range(nsubdomain):
        frami_facej.append (map(int,result[pos:pos+frami_nface[i]]))
        pos += frami_nface[i]
        
    print "Render finished"
    
    line_group = rs.AddGroup("line")
    poly_group = rs.AddGroup("poly")
    
    for face in facep:
        pl = []
        for pn in face:
            pl.append (xyzp[pn-1])
        plane = GetPlane(pl)
        pl.append (xyzp[face[0]-1])
        pl = map(lambda x:rs.PlaneClosestPoint(plane,x),pl)
        id = rs.AddPolyline(pl)
        rs.AddObjectToGroup(id,line_group)
        rs.AddObjectToGroup(rs.AddPlanarSrf(id),poly_group)

DrawDomain ("model_domain.dat")