from pyplasm import *
import scipy
from scipy import *

#---------------------------------------------------------
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


dom1D = INTERVALS(1)(32)

def mapBezier(points):
    curve = BEZIER(S1)(points)
    return MAP(curve)(dom1D)


def mapHermite(points):
    curve = CUBICHERMITE(S1)(points)
    return MAP(curve)(dom1D)


def trasla (p,v):
	q = []
	length=len(p)
	for i in range(length):
		q += [ADD([p[i],v])]
	return q;

def spessore(p0,s):
	c0 = BEZIER(S1)(p0)
	p00 = trasla(p0,s)
	c00 = BEZIER(S1)(p00)
	c = BEZIER(S2)([c0,c00])
	return c

def forma2d(p0,p1):
	c0 = BEZIER(S1)(p0)
	c1 = BEZIER(S1)(p1)
	c = BEZIER(S2)([c0,c1])
	return c


################################### EXERCISE 2 #####################################################################


########################### y = 0 plan ###################################################################

y0 = mapBezier([[-8.9,0,-0.3],[-6,0,1],[-3.6,0,0.7],[-3.4,0,1]])
y1 = mapBezier([[-3.4,0,1],[2.6,0,3.9],[4.5,0,1.55],[8,0,0.9]])
y2 = mapBezier([[8,0,0.9],[9.1,0,1.4],[9.2,0,1.2],[8.4,0,-1.4]])
y3 = mapBezier([[-8.9,0,-0.3],[-8.9,0,-0.7],[-8.6,0,-1.1],[-8.2,0,-1.2]])
y4 = mapBezier([[-8.2,0,-1.2],[-8.4,0,-1.8]])
y5 = mapBezier([[-8.4,0,-1.8],[-6.6,0,-1.8]])
y6 = mapHermite([[-6.6,0,-1.8],[-3.1,0,-1.8],[0.7,0,8],[0.7,0,-8]])
y7 = mapBezier([[-3.1,0,-1.8],[4.1,0,-1.8]])
y8 = mapBezier([[8.4,0,-1.4],[7.3,0,-1.7]])
y9 = mapHermite([[7.3,0,-1.7],[4.1,0,-1.8],[0,0,8],[0,0,-8]])

profileY0 = STRUCT([y0,y1,y2,y3,y4,y5,y6,y7,y8,y9]) 

########################## x = 0 plan ####################################################################

x0 = mapBezier([[0,0,2.2],[0,1.9,2.2],[0,2.4,2.2],[0,3.2,1]])
x1 = mapBezier([[0,3.2,1],[0,4,0.55],[0,4.1,1.1],[0,3.6,-1.7]])
x2 = mapBezier([[0,3.6,-1.7],[0,-3.6,-1.7]])
x3 = mapBezier([[0,-3.2,1],[0,-4,0.55],[0,-4.1,1.1],[0,-3.6,-1.7]])
x4 = mapBezier([[0,0,2.2],[0,-1.9,2.2],[0,-2.4,2.2],[0,-3.2,1]])

profileX0 = STRUCT([x0,x1,x2,x3,x4])

######################### z = 0 plan #####################################################################

z0 = mapBezier([[-8.9,-2.6,0],[-8.4,-3.6,0],[-5.7,-4.3,0],[-3.5,-4.1,0]])
z1 = mapBezier([[-3.5,-4.1,0],[-0.8,-3.5,0],[0.8,-3.5,0],[3.5,-4.1,0]])
z2 = mapBezier([[3.5,-4.1,0],[5.7,-4.3,0],[8.4,-3.6,0],[8.9,-2.6,0]])
z3 = mapBezier([[-8.9,-2.6,0],[-8.9,2.6,0]])
z4 = mapBezier([[-8.9,2.6,0],[-8.4,3.6,0],[-5.7,4.3,0],[-3.5,4.1,0]])
z5 = mapBezier([[-3.5,4.1,0],[-0.8,3.5,0],[0.8,3.5,0],[3.5,4.1,0]])
z6 = mapBezier([[3.5,4.1,0],[5.7,4.3,0],[8.4,3.6,0],[8.9,2.6,0]])
z7 = mapBezier([[8.9,-2.6,0],[8.9,2.6,0]])

profileZ0 = STRUCT([z0,z1,z2,z3,z4,z5,z6,z7])

exercise2 = STRUCT([profileX0,profileY0,profileZ0])



######################## EXERCISE 3 #############################################################################


dom2D = GRID([20,20])

r0 = [[0,1.3,0],[1.3,0,0],[2.4,0,0],[0,-2.4,0]]
r00 = CUBICHERMITE(S1)(r0)

r1 = [[0,0.8,0],[0.8,0,0],[1.4,0,0],[0,-1.4,0]]
r11 = CUBICHERMITE(S1)(r1)

r01 = CUBICHERMITE(S2)([r11,r00,[0,0,2.2],[0,0,-2.2]])
w01 = MAP(r01)(dom2D)

r01b = CUBICHERMITE(S2)([r00,r11,[0,0,-2.2],[0,0,2.2]])
w01b = MAP(r01b)(dom2D)

wheel0 = STRUCT([w01b,w01])
w1 = R([1,2])(PI/2)(wheel0)
w2 = R([1,2])(PI/2)(w1)
w3 = R([1,2])(PI/2)(w2)

wheel1 = COLOR(BLACK)(STRUCT([wheel0,w1,w2,w3]))

r2 = [[0,0.55,0],[0.55,0,0],[1,0,0],[0,-1,0]]
r22 = CUBICHERMITE(S1)(r2)

r12 = CUBICHERMITE(S2)([r22,r11,[0,0,2],[0,0,-2]])
w12 = MAP(r12)(dom2D)

r12b = CUBICHERMITE(S2)([r11,r22,[0,0,-2],[0,0,2]])
w12b = MAP(r12b)(dom2D)

wheel2 = STRUCT([w12b,w12])
w4 = R([1,2])(PI/2)(wheel2)
w5 = R([1,2])(PI/2)(w4)
w6 = R([1,2])(PI/2)(w5)

wheel3 = COLOR(GRAY)(STRUCT([wheel2,w4,w5,w6]))

wheel4 = COLOR(GRAY)(T([3])(-0.5)(CYLINDER([0.2,1])(36)))

w7 = POLYLINE([[-0.19,0.2],[-0.1,0.6],[0.1,0.6],[0.2,0.2],[-0.19,0.2]])
w8 = SOLIDIFY(w7)

wheel5 = COLOR([0.9,0.9,0.9,1])(T([3])(-0.5)(PROD([w8,Q(0.9)])))

wheel6 = STRUCT([wheel5,R([1,2])(2*PI/5)] * 5)

w9 = [[-0.2,0],[-0.17,0.25],[0.17,0.25],[0.2,0]]
sw9 = BEZIERSTRIPE([w9,0.1,20])

w10 = [[0.2,0],[0.17,-0.25],[-0.17,-0.25],[-0.2,0]]
sw10 = BEZIERSTRIPE([w10,0.1,20])

wheel7 = T([3])(-0.5)(PROD([STRUCT([sw9,sw10]),Q(0.9)]))
wheel8 = STRUCT([wheel1,wheel3,wheel4,wheel6,wheel7])
wheel9 = R([2,3])(PI/2)(wheel8)

wheel_0 = T([1,2,3])([-4.8,-3.5,-1.4])(wheel9)
wheel_1 = T([1,2,3])([5.7,-3.5,-1.4])(wheel9)
wheel_2 = T([2])(7)(wheel_0)
wheel_3 = T([2])(7)(wheel_1)

wheels = STRUCT([wheel_0,wheel_1,wheel_2,wheel_3])

exercise3 = STRUCT([exercise2,wheels])

############################## EXERCISE 4 ############################################################



############################## STEERING WHEEL #############################################

dom2D = GRID([20,20])

s0 = [[0,0.5,0],[0.5,0,0],[0.9,0,0],[0,-0.9,0]]
s00 = CUBICHERMITE(S1)(s0)

s1 = [[0,0.4,0],[0.4,0,0],[0.7,0,0],[0,-0.7,0]]
s11 = CUBICHERMITE(S1)(s1)

s01 = CUBICHERMITE(S2)([s11,s00,[0,0,0.2],[0,0,-0.2]])
surface01 = MAP(s01)(dom2D)

s01b = CUBICHERMITE(S2)([s00,s11,[0,0,-0.2],[0,0,0.2]])
surface01b = MAP(s01b)(dom2D)

sw0 = STRUCT([surface01b,surface01])
sw1 = R([1,2])(PI/2)(sw0)
sw2 = R([1,2])(PI/2)(sw1)
sw3 = R([1,2])(PI/2)(sw2)

steering_wheel0 = COLOR(BLACK)(STRUCT([sw0,sw1,sw2,sw3]))

s4 = [[0.45,0.07],[0.1,0.17],[-0.1,0.17],[-0.45,0.07]]
sw4 = BEZIERSTRIPE([s4,0.1,20])

s5 = [[0.06,-0.45],[0.06,-0.05],[0.13,-0.02],[0.45,-0.02]]
sw5 = BEZIERSTRIPE([s5,0.1,20])

s6 = [[-0.45,-0.02],[-0.13,-0.02],[-0.06,-0.05],[-0.06,-0.45]]
sw6 = BEZIERSTRIPE([s6,0.1,20])

centro = COLOR(BLACK)(CYLINDER([0.13,0.06])(36))

steering_wheel1 = COLOR([0.25,0.25,0.25,1])(PROD([STRUCT([sw4,sw5,sw6]),Q(0.04)]))
steering_wheel2 = STRUCT([steering_wheel1,centro])
steering_wheel3 = R([1,2])(PI/2)(steering_wheel2)
steering_wheel4 = S([1,2,3])([1.4,1.4,1.4])(STRUCT([steering_wheel0,steering_wheel3]))
steering_wheel5 = R([1,3])(-PI/3)(steering_wheel4)

steering_wheel = T([1,2,3])([-1.7,-1.2,0.6])(steering_wheel5)

exercise4 = STRUCT([exercise3,steering_wheel]) 



############################ EXERCISE 5 #####################################################################

g = [0.784,0,0,1]

################## PARTE SUPERIORE #####################################################

top1 = [[8,3.5,0.9],[4.3,3.5,2.2],[1.7,3.5,3.1],[-0.5,1.7,2.1]]
top1c = BEZIER(S1)(top1)

top2 = [[8,3.5,0.9],[2.1,3.5,1.1],[6.4,3.5,2],[-0.5,1.7,2]]
top2c = BEZIER(S1)(top2)

top11 = [[8,-3.5,0.9],[4.3,-3.5,2.2],[1.7,-3.5,3.1],[-0.5,-1.7,2.1]]
top11c = BEZIER(S1)(top11)

top22 = [[8,-3.5,0.9],[2.1,-3.5,1.1],[6.4,-3.5,2],[-0.5,-1.7,2]]
top22c = BEZIER(S1)(top22)

top = COLOR(g)(MAP(BEZIER(S2)([top2c,top1c,top11c,top22c]))(GRID([20,20])))


################ FIANCATA 2D ###########################################################

a0 = [[-8.5,0,-1.8],[-6.6,0,-1.8]]
a1 = [[-8.2,0,-1.2],[-6.6,0,-1.2]]
a01 = MAP(forma2d(a0,a1))(GRID([20,20]))
a10 = MAP(forma2d(a1,a0))(GRID([20,20]))

a2 = [[-8.2,0,-1.2],[-9.5,0,-0.1],[-9.3,0,-0.5],[-6.6,0,0.5]]
a21 = MAP(forma2d(a1,a2))(GRID([20,20]))
a12 = MAP(forma2d(a2,a1))(GRID([20,20]))

a3 = [[-6.6,0,0.5],[-5.2,0,0.75],[-4.6,0,0.85],[-3.1,0,0.9]]
ac3 = BEZIER(S1)(a3)
a4 = [[-6.6,0,-1.8],[-3.1,0,-1.8],[0.7,0,8],[0.7,0,-8]]
ac4 = CUBICHERMITE(S1)(a4)
a34 = MAP(BEZIER(S2)([ac4,ac3]))(GRID([20,20]))
a43 = MAP(BEZIER(S2)([ac3,ac4]))(GRID([20,20]))

a5 = [[-3.1,0,0.9],[-0.6,0,1],[1.7,0,1.1],[4.1,0,1.1]]
a6 = [[-3.1,0,-1.8],[4.1,0,-1.8]]
a56 = MAP(forma2d(a6,a5))(GRID([20,20]))
a65 = MAP(forma2d(a5,a6))(GRID([20,20]))

a7 = [[4.1,0,1.1],[5.1,0,1.1],[6.3,0,1.1],[7.3,0,1]]
ac7 = BEZIER(S1)(a7)
a8 = [[4.1,0,-1.8],[7.3,0,-1.7],[0,0,8],[0,0,-8]]
ac8 = CUBICHERMITE(S1)(a8)
a78 = MAP(BEZIER(S2)([ac8,ac7]))(GRID([20,20]))
a87 = MAP(BEZIER(S2)([ac7,ac8]))(GRID([20,20]))

a9 = [[7.3,0,1],[7.6,0,1],[7.8,0,1],[8,0,0.9]]
a_10 = [[7.3,0,-1.7],[8.4,0,-1.4]]
a910 = MAP(forma2d(a_10,a9))(GRID([20,20]))
a109 = MAP(forma2d(a9,a_10))(GRID([20,20]))

a11 = [[7.9,0,0.9],[9.1,0,1.4],[9.2,0,1.2],[8.4,0,-1.4]]
a_12 = [[8.4,0,-1.4],[7.9,0,0.9]]
a1112 = MAP(forma2d(a_12,a11))(GRID([20,20]))
a1211 = MAP(forma2d(a11,a_12))(GRID([20,20]))

profileY02D_1 = T([2])(-3.5)(STRUCT([a01,a21,a34,a56,a78,a910,a1112]))

profileY02D_2 = T([2])(3.5)(STRUCT([a10,a12,a43,a65,a87,a109,a1211]))

profileY02D = COLOR(g)(STRUCT([profileY02D_1,profileY02D_2]))



########################## PARTE ANTERIORE E POSTERIORE ################################


b0 = [[-8.2,-3.5,-1.2],[-9.5,-3.5,-0.1],[-9.3,-3.5,-0.5],[-6.6,-3.5,0.5]]
b00 = COLOR(g)(MAP(spessore(b0,[0,7,0]))(GRID([20,20])))

b1 = [[-8.2,3.5,-1.2],[-8.7,3.5,-1.8],[-8.4,3.5,-1.8]]
b11 = COLOR(BLACK)(MAP(spessore(b1,[0,-7,0]))(GRID([20,20])))

b2 = [[-6.6,-3.5,0.5],[-5.2,-3.5,0.75],[-4.6,-3.5,0.85],[-3.1,-3.5,0.9]]
b22 = COLOR(g)(MAP(spessore(b2,[0,7,0]))(GRID([20,20])))

profileX02D_ant = STRUCT([b00,b11,b22])

d1 = [[7.9,-3.5,0.9],[9.1,-3.5,1.4],[9.2,-3.5,1.2],[8.4,-3.5,-1.4]]

profileX02D_post = COLOR(BLACK)(MAP(spessore(d1,[0,7,0]))(GRID([20,20])))

profileX02D = STRUCT([profileX02D_ant,profileX02D_post])


####################### PARTE INFERIORE ##############################################

i1 = [[4.1,-3.5,-1.8],[-3.1,-3.5,-1.8]]
i11 = MAP(spessore(i1,[0,1,0]))(GRID([20,20]))

i2 = [[-8.4,-3.5,-1.8],[-6.6,-3.5,-1.8]]
i22 = MAP(spessore(i2,[0,7,0]))(GRID([20,20]))

i3 = [[7.3,3.5,-1.7],[8.4,3.5,-1.4]]
i33 = MAP(spessore(i3,[0,-7,0]))(GRID([20,20]))

i4 = [[-8.4,-2.5,-1.8],[7.3,-2.5,-1.8]]
i44 = MAP(spessore(i4,[0,5,0]))(GRID([20,20]))

i5 = [[4.1,3.5,-1.8],[-3.1,3.5,-1.8]]
i55 = MAP(spessore(i5,[0,-1,0]))(GRID([20,20]))

profileZ02D = COLOR([0.25,0.25,0.25,1])(STRUCT([i11,i22,i33,i44,i55])) 

######################################################################################

scocca = STRUCT([top,profileY02D,profileX02D,profileZ02D])

#######################################################################################


###################### FINESTRINI E PARABREZZA ########################################

# bordo finestrini

f0 = [[5.2,-3.5,1.1],[3.8,-3.2,1.6],[1.3,-2.4,2.8],[-0.5,-1.7,2]]
f1 = [[4.6,-3.5,1.1],[2.7,-3.2,1.6],[1.8,-2.4,2.2],[-0.5,-1.7,1.9]]
f = MAP(forma2d(f0,f1))(GRID([20,20]))

f00 = [[5.2,3.5,1.1],[3.8,3.2,1.6],[1.3,2.4,2.8],[-0.5,1.7,2]]
f11 = [[4.6,3.5,1.1],[2.7,3.2,1.6],[1.8,2.4,2.2],[-0.5,1.7,1.9]]
f1 = MAP(forma2d(f00,f11))(GRID([20,20]))

f2 = [[-0.5,1.8,1.9],[-0.5,-1.8,1.9]]
f22 = [[-0.5,1.8,2.1],[-0.5,-1.8,2.1]]
f_2 = MAP(forma2d(f2,f22))(GRID([20,20]))

f3 = [[-3.1,3.5,0.9],[-0.5,1.6,2]]
f33 = [[-3,3.5,0.9],[-0.3,1.7,1.9]]
f_3 = MAP(forma2d(f3,f33))(GRID([20,20]))

f4 = [[-3.1,-3.5,0.9],[-0.5,-1.6,2]]
f44 = [[-3,-3.5,0.9],[-0.3,-1.7,1.9]]
f_4 = MAP(forma2d(f4,f44))(GRID([20,20]))

f5 = [[-1.4,3.5,0.9],[-0.5,1.7,1.9]]
f55 = [[-1.2,3.5,0.9],[-0.3,1.7,1.9]]
f_5 = MAP(forma2d(f5,f55))(GRID([20,20]))

f6 = [[2.1,-3.5,1.1],[2.3,-1.8,2]]
f66 = [[1.9,-3.3,1.1],[2.1,-1.6,2]]
f_6 = MAP(forma2d(f6,f66))(GRID([20,20]))

f7 = [[2.1,3.5,1.1],[2.3,1.8,2]]
f77 = [[1.9,3.3,1.1],[2.1,1.6,2]]
f_7 = MAP(forma2d(f7,f77))(GRID([20,20]))

f8 = [[-0.5,-1.8,1.9],[-1.4,-3.5,0.9]]
f88 = [[-0.3,-1.6,1.9],[-1.2,-3.3,0.9]]
f_8 = MAP(forma2d(f8,f88))(GRID([20,20]))

bordo_fin = COLOR(BLACK)(STRUCT([f,f1,f_2,f_3,f_4,f_5,f_6,f_7,f_8]))

# parabrabrezza

v0 = [[-3.1,-3.5,0.9],[-0.5,-1.6,2]]
v00 = [[-3.1,3.5,0.9],[-0.5,1.6,2]]
parabrezza = COLOR([0.6863,1,0.6745,0.1])(MAP(forma2d(v0,v00))(GRID([20,20])))



#########################################################################################

ferrari = STRUCT([wheels,steering_wheel,scocca,bordo_fin,parabrezza])

exercise5 = STRUCT([ferrari,exercise2])

VIEW(ferrari)
VIEW(exercise5)




