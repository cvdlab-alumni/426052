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


#funzione che trasla alcuni punti della curva in funzione del punto in cui si trova la curva

#points : array di array con i punti corrispondenti alle 15 barre verticali
#k1 e k2 : da quale punto a quale punto modificare ogni curva
#t : parametro che viene incrementato per modificare in modo parametrico ogni coppia di curve
#d : coordinata da modificare 

def trasla_param(points,k1,k2,t,d):
  p = points
  j = t
  x = 2
  q = k2-k1
  for i in range(7):
    z = k1
    for k in range(q):
      if x==14:
        p[x][z][d] = p[x-2][z][d]
        p[x+1][z][d] = p[x-2+1][z][d]
      else:
        p[x][z][d] = p[x][z][d]+ j
        p[x+1][z][d] = p[x+1][z][d]+ j
        p[len(points)-2-x][z][d] = p[len(points)-2-x][z][d]+ j
        p[len(points)-2-x+1][z][d] = p[len(points)-2-x+1][z][d]+ j
      z = z + 1 
    j = j + t
    x = x + 2
  return p


def trasla_y(points):
  p = []
  for i in range(30):
    p += [trasla(points,[0,i*(0.227586207),0])]
  return p

def trasla_z(points,t):
  p = points
  for i in range(28):
    if i>1:
      p[i][0][2] = p[i][0][2] + t
  return p


#/////////////////////////////////////////////////// 

def horizontal(p1,p2):
  cp1 = BEZIER(S1)(p1)
  p11 = trasla(p1,[-0.25,0,0])
  cp11 = BEZIER(S1)(p11)
  p_0 = MAP(BEZIER(S2)([cp1,cp11]))(dom2D)
  cp2 = BEZIER(S1)(p2)
  p22 = trasla(p2,[-0.25,0,0])
  cp22 = BEZIER(S1)(p22)
  p_1 = MAP(BEZIER(S2)([cp2,cp22]))(dom2D)
  p_2 = MAP(BEZIER(S2)([cp1,cp2]))(dom2D)
  p_3 = MAP(BEZIER(S2)([cp11,cp22]))(dom2D)
  p = STRUCT([p_0,p_1,p_2,p_3])
  return p


def repeat_horizontal(f):
  fig = f
  for i in range(15):
    f1 = T([2])([i*(0.454)])(f)
    fig = STRUCT([fig,f1])
  return fig


def vertical(p):
  v = [0, 1.5966490409604734, 2.1694314093789644, 2.5439012244975236, 2.8007660666324847, 2.9706908287467413, 3.0681547809717813, 3.1, 3.0681547809717813, 2.9706908287467413, 2.8007660666324847, 2.5439012244975236, 2.1694314093789644, 1.5966490409604734, 0]
  p0 = trasla_y(p)                    #//genero le curve tutte alla stessa altezza a partire da quella data
  p1 = trasla_x0(p0,v)                #//traslo le x dei primi punti di ogni curva per ottenere la base circolare
  p2 = trasla_param(p1,1,2,0.15,2)    #//modifico le z di alcuni punti in modo parametrico 
  p3 = trasla_param(p2,3,4,-0.05,0)   #//modifico le x di alcuni punti in modo parametrico 
  p4 = trasla_param(p3,6,9,-0.05,0)   #//modifico le x di alcuni punti in modo parametrico 
  p5 = trasla_param(p4,6,9,-0.1,2)    #//modifico le z di alcuni punti in modo parametrico 
  p6 = trasla_z(p5,0.25)              #//modifico le z dei primi punti di ogni curva di una quantita prefissata 
  return p6



#//////////////////////////////////////////////////

#//funzione che genera il vettore traslazione da applicare alle curve per ottenere la base circolare 

#def vett_traslX():   
#  p = [-3.1,-2.6572,-2.2144,-1.7716,-1.3288,-0.886,-0.4432,0,0.4432,0.886,1.3288,1.7716,2.2144,2.6572,3.1]
#  v = []
  #v[0] = 0
#  v += [0]
#  for i in range(14):
#    c = 9.61 - p[i]*p[i]
#    v += [math.sqrt(c)]
#  v += [0]
#  return v


#// funzione che utilizza il vettore generato da vett_traslX e trasla le curve 

def trasla_x0(points,v):
  p = points
  j = 1
  x = 2
  for i in range(13):
    p[x][0][0] = p[x][0][0] - v[j] + 0.12
    p[x+1][0][0] = p[x+1][0][0] - v[j] + 0.12
    j = j + 1
    x = x + 2
  p[27][0][0] = p[27][0][0] - v[j]+0.12
  return p



#/////////////////////////////////////////////////

def gen_nubs(points):
  p = []
  for i in range(30):
    knots = nodes(points[i],11) 
    p += [NUBSPLINE(2)(knots)(points[i])]
  return p


def gen_map(p):
  vertical = p[0]
  for i in range(30):
    v = p[i]
    vertical = STRUCT([vertical,v])
  return vertical


def arc(alpha,r,R):
  domain1 = INTERVALS(alpha)(36)
  domain2 = INTERVALS(R)(1)
  domain3 = INTERVALS(r)(1)
  domain_0 = PROD([domain1,domain2])
  domain_1 = PROD([domain1,domain3])
  def mapping(v):
    a = v[0]
    r = v[1]
    return [r*COS(a), r*SIN(a)]
  model1 = MAP(mapping)(domain_0)
  model2 = MAP(mapping)(domain_1)
  model = DIFF([model1,model2])
  return model



#////////////////////////////////////////////////////////////////////////////////////////////////////////////

black = [0,0,0,1]


v = [0, 1.5966490409604734, 2.1694314093789644, 2.5439012244975236, 2.8007660666324847, 2.9706908287467413, 3.0681547809717813, 3.1, 3.0681547809717813, 2.9706908287467413, 2.8007660666324847, 2.5439012244975236, 2.1694314093789644, 1.5966490409604734, 0]

a0 = [[4.8,-3.3,0],[6,-3.3,1.95],[6.35,-3.3,3.35],[6.05,-3.3,4.25],[4.85,-3.3,4.85],
[3.6,-3.3,4.85],[2.55,-3.3,4.8],[1.7,-3.3,5.2],[1.1,-3.3,6.3],[0.7,-3.3,7],[0.5,-3.3,8]]

a1 = [[4.4,-3.3,0],[5.6,-3.3,1.95],[5.95,-3.3,3.35],[5.7,-3.3,4.1],[4.85,-3.3,4.45],
[3.6,-3.3,4.45],[2.55,-3.3,4.4],[1.25,-3.3,5.05],[0.7,-3.3,6.3],[0.3,-3.3,7],[0.1,-3.3,8]]

p0_0 = vertical(a0)
p0_1 = gen_nubs(p0_0)
p0 = gen_map(p0_1)

p1_0 = vertical(a1)
p1_1 = gen_nubs(p1_0)
p1 = gen_map(p1_1)

p2_0 = T([1,2,3])([0.1,-3.3,8])(CUBOID([0.4,0.227586207,0.05]))
p2 = repeat_horizontal(p2_0)

vertical_0 = STRUCT([p0,p1,p2])

a2 = [[4.9,-3.3,4.45],[4.9,-3.2,4.05],[4.9,3.2,4.05],[4.9,3.3,4.45]]

a3 = [[4.9,-3.3,4.45],[4.9,-3.2,3.85],[4.9,3.2,3.85],[4.9,3.3,4.45]]

p4 = horizontal(a2,a3)

p5_0 = T([1,2,3])([0.1,-3.3,7.772])(R([1,3])(PI/20)(CUBOID([0.1,6.6,0.227])))
p5_1 = T([1,3])([0.69,-0.454])(R([1,3])(PI/40)(p5_0))
p5_2 = T([1,3])([0.49,-0.452])(R([1,3])(PI/60)(p5_1))
p5 = T([3])(0.04)(STRUCT([p5_0,p5_1,p5_2]))

p6_0 = PROD([arc(PI,3.073,3.3),Q(0.4)])
p6 = T([1])(4.65)(R([1,2])(PI/2)(p6_0))

horizontal_0 = STRUCT([p4,p5,p6])

pantosh_chair = COLOR(black)(STRUCT([vertical_0,horizontal_0]))

#VIEW(pantosh_chair)