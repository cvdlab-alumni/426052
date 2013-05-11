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

################################## EXERCISE 2 ############################################################

########################### y = 0 plan ####################################################

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

########################## x = 0 plan #####################################################

x0 = mapBezier([[0,0,2.2],[0,1.9,2.2],[0,2.4,2.2],[0,3.2,1]])
x1 = mapBezier([[0,3.2,1],[0,4,0.55],[0,4.1,1.1],[0,3.6,-1.7]])
x2 = mapBezier([[0,3.6,-1.7],[0,-3.6,-1.7]])
x3 = mapBezier([[0,-3.2,1],[0,-4,0.55],[0,-4.1,1.1],[0,-3.6,-1.7]])
x4 = mapBezier([[0,0,2.2],[0,-1.9,2.2],[0,-2.4,2.2],[0,-3.2,1]])

profileX0 = STRUCT([x0,x1,x2,x3,x4])

######################### z = 0 plan ######################################################

z0 = mapBezier([[-8.9,-2.6,0],[-8.4,-3.6,0],[-5.7,-4.3,0],[-3.5,-4.1,0]])
z1 = mapBezier([[-3.5,-4.1,0],[-0.8,-3.5,0],[0.8,-3.5,0],[3.5,-4.1,0]])
z2 = mapBezier([[3.5,-4.1,0],[5.7,-4.3,0],[8.4,-3.6,0],[8.9,-2.6,0]])
z3 = mapBezier([[-8.9,-2.6,0],[-8.9,2.6,0]])
z4 = mapBezier([[-8.9,2.6,0],[-8.4,3.6,0],[-5.7,4.3,0],[-3.5,4.1,0]])
z5 = mapBezier([[-3.5,4.1,0],[-0.8,3.5,0],[0.8,3.5,0],[3.5,4.1,0]])
z6 = mapBezier([[3.5,4.1,0],[5.7,4.3,0],[8.4,3.6,0],[8.9,2.6,0]])
z7 = mapBezier([[8.9,-2.6,0],[8.9,2.6,0]])

profileZ0 = STRUCT([z0,z1,z2,z3,z4,z5,z6,z7])

###########################################################################################

exercise2 = STRUCT([profileX0,profileY0,profileZ0])






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

s7 = [[-0.06,0],[-0.05,0.06],[0.05,0.06],[0.06,0]]
sw7 = BEZIERSTRIPE([s7,0.1,20])

s8 = [[0.06,0],[0.05,-0.06],[-0.05,-0.06],[-0.06,0]]
sw8 = BEZIERSTRIPE([s8,0.07,20])

steering_wheel1 = COLOR([0.25,0.25,0.25,1])(PROD([STRUCT([sw4,sw5,sw6,sw7,sw8]),Q(0.04)]))
steering_wheel2 = R([1,2])(PI/2)(steering_wheel1)
steering_wheel3 = S([1,2,3])([1.4,1.4,1.4])(STRUCT([steering_wheel0,steering_wheel2]))
steering_wheel4 = R([1,3])(-PI/3)(steering_wheel3)

steering_wheel = T([1,2,3])([-1.3,-1.2,0.8])(steering_wheel4)

#######################################################################################################

exercise4 = STRUCT([exercise2,steering_wheel]) 
