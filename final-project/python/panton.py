from pyplasm import *
import scipy
from scipy import *

#---------------------------------------------------------

dom1D = INTERVALS(1)(32)
dom2D = PROD([dom1D, dom1D])


def VERTEXTRUDE((V,coords)):
    return CAT(AA(COMP([AA(AR),DISTR]))(DISTL([V,coords])))


def cumsum(iterable):
    # cumulative addition: list(cumsum(range(4))) => [0, 1, 3, 6]
    iterable = iter(iterable)
    s = iterable.next()
    yield s
    for c in iterable:
        s = s + c
        yield s

def larExtrude(model,pattern):
    V,FV = model
    d = len(FV[0])
    offset = len(V)
    m = len(pattern)
    outcells = []
    for cell in FV:
        # create the indices of vertices in the cell "tube"
        tube = [v + k*offset for k in range(m+1) for v in cell]
        # take groups of d+1 elements, via shifting by one
        rangelimit = len(tube)-d
        cellTube = [tube[k:k+d+1] for k in range(rangelimit)]
        outcells += [scipy.reshape(cellTube,newshape=(m,d,d+1)).tolist()]
    outcells = AA(CAT)(TRANS(outcells))
    outcells = [group for k,group in enumerate(outcells) if pattern[k]>0 ]
    coords = list(cumsum([0]+(AA(ABS)(pattern))))
    outVerts = VERTEXTRUDE((V,coords))
    newModel = outVerts, CAT(outcells)
    return newModel

def GRID(args):
    model = ([[]],[[0]])
    for k,steps in enumerate(args):
        model = larExtrude(model,steps*[1])
    V,cells = model
    verts = AA(list)(scipy.array(V) / AA(float)(args))
    return MKPOL([verts, AA(AA(lambda h:h+1))(cells), None])


def nodes(points,s):
	m = s
	k = 2
	n = (m + k + 1)
	l = n - 3
	j = 1
	knots = []
	for i in range(n):
		if (i<3):
			knots += [0]
		if (3<=i<l):
			knots += [j]
			j = j + 1
		if (i>=l):
			knots += [j]
	return knots


def trasla(p,v):
    q = []
    length=len(p)
    for i in range(length):
        q += [ADD([p[i],v])]
    return q


def trasla_fix(points,k1,k2,t,d):
  p = points
  x = 1
  q = k2-k1
  for i in range(11):
    z = k1
    for k in range(q):
        p[x][z][d] = p[x][z][d] + t
        z = z + 1
    x = x + 1
  return p


def trasla_y(points):
    p = []
    for i in range(30):
        p += [trasla(points,[0,i*(0.333333333),0])]
    return p


def trasla_x0(points,v):
    p = points
    for i in range(12):
        p[i][0][0] = p[i][0][0] - v[i]
    return p


def trasla_z12(points,v):
    p = points
    for i in range(12):
        p[i][12][2] = p[i][12][2] + v[i]
        p[i][13][2] = p[i][13][2] + v[i]
    return p


def panton(p):
    v0 = [0,2.3,3.1,3.55,3.8,4,4.1,4,3.8,3.55,3.1,2.3,0]           #vettore per la base
    v1 = [0,0.5,0.6,0.7,0.8,0.85,0.85,0.85,0.8,0.7,0.6,0.5,0]      #vettore per la parte alta
    p0 = trasla_y(p)                   #genero le curve tutte alla stessa y a partire da quella data
    p1 = trasla_x0(p0,v0)              #traslo le x dei primi punti di ogni curva per ottenere la base circolare
    p2 = trasla_fix(p1,1,2,-2,0)       #modifico le x e le z di alcuni punti  
    p3 = trasla_fix(p2,2,3,-0.5,0)     
    p4 = trasla_fix(p3,2,4,0.5,2)      
    p5 = trasla_fix(p4,7,10,-0.25,2)     
    p6 = trasla_z12(p5,v1)             #traslo le z degli ultimi punti per ottenere la parte alta tondeggiante
    p7 = trasla_fix(p6,12,13,-0.15,2)
    p8 = trasla_fix(p7,9,11,-0.3,0)          
    return p8

def gen_nubs(points):
    p = []
    for i in range(13):
        knots = nodes(points[i],14) 
        p += [NUBSPLINE(2)(knots)(points[i])]
    return p


def gen_map(p):
    vertical = p[0]
    for i in range(13):
        v = p[i]
        vertical = STRUCT([vertical,v])
    return vertical




#///////////////////////////////////////////////////////////////////////////////////////////////////////////////

b = [0,0,0,1]

p_in = [[5.05,-2.3,0],[5.55,-2.3,1.05],[5.9,-2.3,2.2],[6,-2.3,3.15],[5.65,-2.3,4.05],[4.95,-2.3,4.525],[3.65,-2.3,4.5],[2.4,-2.3,4.3],[1.85,-2.3,4.55],[1.45,-2.3,5.2],[1.15,-2.3,6.4],[0.6,-2.1,6.9]]
knotspin = nodes(p_in,12) 
cpin = NUBSPLINE(2)(knotspin)(p_in)

p_fin = [[5.05,2.3,0],[5.55,2.3,1.05],[5.9,2.3,2.2],[6,2.3,3.15],[5.65,2.3,4.05],[4.95,2.3,4.525],[3.65,2.3,4.5],[2.4,2.3,4.3],[1.85,2.3,4.55],[1.45,2.3,5.2],[1.15,2.3,6.4],[0.6,2.1,7]]
knotspfin = nodes(p_fin,12) 
cpfin = NUBSPLINE(2)(knotspfin)(p_fin)

p0 = [[4.75,-2,0],[5.25,-2,1.05],[5.6,-2,2.2],[5.7,-2,3.15],[5.35,-2,3.9],[4.95,-2,4.275],[3.65,-2,4.4],[2.4,-2,4.3],[1.85,-2,4.55],[1.45,-2,5.2],[1.15,-2,6.4],[0.9,-2,7.3],[0.55,-2,7.5],[0.45,-2,7.1]]

p_0 = panton(p0)
p_1 = gen_nubs(p_0)

panton_0 = gen_map(p_1)

p1 = [[4.9,-2.5,0],[5.35,-2.5,1.05],[5.65,-2.5,2.2],[5.55,-2.5,3.15],[5.25,-2.5,3.7],[4.95,-2.5,4.075],[3.65,-2.5,4.15],[2.4,-2.5,4.1],[1.85,-2.5,4.4],[1.15,-2.5,5.2],[0.8,-2.5,6.4],[0.6,-2.1,6.9]]
knotsp1 = nodes(p1,12) 
cp1 = NUBSPLINE(2)(knotsp1)(p1)

p2 = [[4.9,2.5,0],[5.35,2.5,1.05],[5.65,2.5,2.2],[5.55,2.5,3.15],[5.25,2.5,3.7],[4.95,2.5,4.075],[3.65,2.5,4.15],[2.4,2.5,4.1],[1.85,2.5,4.4],[1.15,2.5,5.2],[0.8,2.5,6.4],[0.6,2.1,6.9]]
knotsp2 = nodes(p2,12) 
cp2 = NUBSPLINE(2)(knotsp2)(p2)

panton_chair = COLOR(b)(STRUCT([panton_0,cp1,cp2]))

VIEW(panton_chair)