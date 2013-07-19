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


def figure2d(p0,p1):
	c0 = BEZIER(S1)(p0)
	c1 = BEZIER(S1)(p1)
	c = BEZIER(S2)([c0,c1])
	return c


def vertical_back(f,n):
	fig = f
	for i in range(n):
		f1 = R([1,2])(PI/22.5)(fig)
		fig = STRUCT([fig,f1])
	return fig


def horizontal_back(h):
	d = PI-(2*h*PI/22.5)
	domain1 = INTERVALS(d)(36)
	domain2 = INTERVALS(4.21)(1)
	domain3 = INTERVALS(4.11)(1)
	domain_0 = PROD([domain1, domain2])
	domain_1 = PROD([domain1, domain3])
	def mapping(v):
		a = v[0]
		r = v[1]
		return [r*COS(a), r*SIN(a)]
	m = MAP(mapping)(domain_0)
	m1 = PROD([m,Q(0.286)])
	model1 = R([1,2])(h*(PI/22.5))(m1)
	mm = MAP(mapping)(domain_1)
	mm1 = PROD([mm,Q(0.286)])
	model2 = R([1,2])(h*(PI/22.5))(mm1)
	model = DIFF([model1,model2])
	return model


def hor_back():
	domain1 = INTERVALS(PI/22.5)(36)
	domain2 = INTERVALS(4.21)(1)
	domain3 = INTERVALS(4.11)(1)
	domain_0 = PROD([domain1,domain2])
	domain_1 = PROD([domain1,domain3])
	def mapping(v):
		a = v[0]
		r = v[1]
		return [r*COS(a), r*SIN(a)]
	m = MAP(mapping)(domain_0)
	model1 = PROD([m,Q(0.286)])
	model2 = R([1,2])(21.5*PI/22.5)(model1)
	model0 = STRUCT([model1,model2])
	mm = MAP(mapping)(domain_1)
	model11 = PROD([mm,Q(0.286)])
	model22 = R([1,2])(21.5*PI/22.5)(model11)
	model00 = STRUCT([model11,model22])
	model = DIFF([model0,model00])
	return model


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



########### COLORS ////////////////////////////////////////////////////////////////////////////////////////////////

c_black = [0,0,0,1]
c_gray = [0.1,0.1,0.1,1]
c_pillow = [0.87,0.1921,0.3882,1]


########### BACK ///////////////////////////////////////////////////////////////////////////////////////////////

#### Vertical Back

v0 = T([1])([4.1])(CUBOID([0.1,0.286,11.9]))

half_vertical0 = vertical_back(v0,11)
half_vertical1 = T([3])([11.9])(R([1,3])(PI)(half_vertical0))

v1 = T([2])([-0.286])(v0)
v2 = T([1])([-8.3])(v1)


vertical = STRUCT([half_vertical0,half_vertical1,v1,v2])

#### Horizontal Back

h_0 = horizontal_back(0)
h0 = T([3])([11.614])(h_0)
h1 = T([3])([11.059])(h_0)
h2 = T([3])([10.504])(horizontal_back(2))
h3 = T([3])([9.949])(horizontal_back(3))
h4 = T([3])([9.394])(horizontal_back(4))
h5 = T([3])([8.839])(horizontal_back(5))
h6 = T([3])([8.284])(horizontal_back(6))
h7 = T([3])([7.729])(horizontal_back(7))
h8 = T([3])([7.174])(horizontal_back(8))
h9 = T([3])([6.619])(horizontal_back(9))

horizontal_0 = STRUCT([h0,h1,h2,h3,h4,h5,h6,h7,h8,h9])

h_1 = T([3])([6.064])(horizontal_back(10))

horizontal_1 = STRUCT([h_1, T([3])([-0.555])] * 5)


h_2 = T([3])([10.504])(hor_back())

horizontal_2 = STRUCT([h_2, T([3])([-0.555])] * 13)


horizontal = STRUCT([horizontal_0,horizontal_1,horizontal_2])

####################


back = COLOR(c_black)(STRUCT([vertical,horizontal]))


########### BASE ///////////////////////////////////////////////////////////////////////////////////////////////

b0 = [[-4.1,-0.286,0],[-3.5,0,0]]
b1 = [[-4.1,-0.286,0.2],[-3.5,0,0.2]]
b2 = [[-4.1,0,0.2],[-3.5,0,0.2]]

b_0 = MAP(figure2d(b0,b1))(dom2D)
b_1 = MAP(figure2d(b1,b2))(dom2D)

b3 = [[4.1,-0.286,0],[3.5,0,0]]
b4 = [[4.1,-0.286,0.2],[3.5,0,0.2]]
b5 = [[4.1,0,0.2],[3.5,0,0.2]]

b_2 = MAP(figure2d(b3,b4))(dom2D)
b_3 = MAP(figure2d(b4,b5))(dom2D)

base_0 = STRUCT([b_0,b_1,b_2,b_3])

base_1 = PROD([arc(PI,0,4.1),Q(0.2)])

base = COLOR(c_black)(STRUCT([base_0,base_1]))



######### SITTING ///////////////////////////////////////////////////////////////////////////////////////////////

s0 = [[-4,0.1,0.2],[4,0.1,0.2]]
s1 = [[-4,-0.5,3.6],[4,-0.5,3.6]]

s_0 = MAP(figure2d(s0,s1))(dom2D)

s2 = [[-4,0.1,0.2],[-4,-0.5,3.6]]
s3 = [[-4,0.1,0.2],[-4,0.1,3.6]]

s_1 = MAP(figure2d(s2,s3))(dom2D)
s_2 = T([0])([8])(s_1) 

s_3 = T([1,2])([0.1,0.2])(PROD([arc(PI,0,4),Q(3.4)]))

s4 = [[-4,0.1,3.6],[4,0.1,3.6]]
s5 = [[-4,-0.5,3.6],[4,-0.5,3.6]]

s_4 = MAP(figure2d(s4,s5))(dom2D)

sitting = COLOR(c_gray)(STRUCT([s_0,s_1,s_2,s_3,s_4]))


######### PILLOW //////////////////////////////////////////////////////////////////////////////////

p0 = [[-4,0.1,3.6],[-4,-0.7,3.5],[-4,-0.7,3.844],[-4,0.1,3.844]]
p1 = [[4,0.1,3.6],[4,-0.7,3.5],[4,-0.7,3.844],[4,0.1,3.844]]
p2 = [[-4,0.1,3.6],[-4,0.1,3.844]]

pillow_0 = MAP(figure2d(p1,p0))(dom2D)
pillow_1 = MAP(figure2d(p0,p2))(dom2D)
pillow_2 = T([0])([8])(pillow_1)

pillow_3 = T([2,3])([0.1,3.6])(PROD([arc(PI,0,4),Q(0.244)]))

pillow = COLOR(c_pillow)(STRUCT([pillow_0,pillow_1,pillow_2,pillow_3]))


######### WILLOW CHAIR /////////////////////////////////////////////////////////////////////////////

willow = STRUCT([back,base,sitting,pillow])
willow_chair = S([1])([1.146341463])(willow)

#VIEW(willow_chair)
