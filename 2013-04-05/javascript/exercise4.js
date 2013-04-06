//Utility functions

//T(dims)(values)(object)

T = function (dims) {
	dims = dims.map(function(dim){
		return dim - 1;
	});
    return function (values) {
      return function (object) {
       return object.clone().translate(dims, values);
      };
    };
  }

//R
R = function (dims) {
	dims = dims.map(function(dim){
		return dim - 1;
	});
    return function (angle) {
      return function (object) {
        return object.clone().rotate(dims, angle);
      };
    };
  }

//S
S = function (dims) {
	dims = dims.map(function(dim){
		return dim - 1;
	});
    return function (values) {
      return function (object) {
        return object.clone().scale(dims, values);
      };
    };
  }


//Per scambiare le coordinate in pyplasm -> S1,S2,S3
//Per scambiare le coordinate in Plasm.js -> S0,S1,S2
S3=S2
S2=S1
S1=S0

//Funzioni da pyplasm a Plasm.js
GRID = SIMPLEXGRID
VIEW = DRAW
NN = REPLICA

//arco di circonferenza: il dominio invece che da 0 a 2*PI, sarà tra 0 e alpha (angolo desiderato)
function arc(alpha,r,R){
	var domain = DOMAIN([[0,alpha], [r,R]])([36,1]);
	var mapping = function(v){
		var a = v[0];
		var r = v[1];
		return [r*COS(a), r*SIN(a)];
	};
	var model = MAP(mapping)(domain);
	return model;
};
DRAW(arc(PI/3,2,3));
GRID = SIMPLEX_GRID

//Exercise4

//Windows nord
frameVertNorth1 = GRID([[0.25],[-0.25,0.1,-1.1,0.1,-1.15,0.1,-1.15,0.1,-1.1,0.1],[-2.5,-0.95,1.3,-1.2,1.3,-1.2,1.3]])

frameHorizNorth1 = GRID([[0.25],[-0.25,5],[-2.5,-0.85,0.1,-1.3,0.1,-1,0.1,-1.3,0.1,-1,0.1,-1.3,0.1]])

frameVertNorth2 = GRID([[0.25],[-6.1,0.05,-0.175,0.05,-0.175,0.05],[-2.5,-0.15,2.3,-0.15,2.3,-0.15,2.3]])

frameHorizNorth2 = GRID([[0.25],[-6.1,0.5],[-2.5,-0.15,0.05,-1.075,0.05,-1.075,0.05,-0.15,0.05,-1.075,0.05,-1.075,0.05,-0.15,0.05,-1.075,0.05,-1.075,0.05]])


glassNorth1 = GRID([[0.25],[-0.25,-0.1,1.1,-0.1,1.15,-0.1,1.15,-0.1,1.1],[-2.5,-0.95,1.3,-1.2,1.3,-1.2,1.3]])

glassNorth2 = GRID([[0.25],[-6.1,-0.05,0.175,-0.05,0.175],[-2.5,-0.15,-0.05,1.075,-0.05,1.075,-0.05,-0.15,-0.05,1.075,-0.05,1.075,-0.05,-0.15,-0.05,1.075,-0.05,1.075]])

glassNorth = COLOR([0,110,110,0.2])(STRUCT([glassNorth1,glassNorth2]))

framesNorth = COLOR([0,0,0])(STRUCT([frameHorizNorth1,frameVertNorth1,frameVertNorth2,frameHorizNorth2]))

windowsNorth = T([1])([11])(STRUCT([framesNorth,glassNorth]))

//Railing South

railingSouth = COLOR([0,0,0])(GRID([[0.05],[1.1,-0.25,5],[-8.3,0.05,-0.25,0.05]]))


//Windows South

frameVertSouth1 = GRID([[0.25],[-1.6,0.1,-1.075,0.1,-1.175,0.1,-1.175,0.1,-1.075,0.1],[-2.5,-0.14,-0.1,1.03,-0.1,1.03,-0.1,-0.14,-0.1,1.03,-0.1,1.03]])

frameHorizSouth1 = GRID([[0.25],[-1.6,5],[-2.5,-0.14,0.1,-1.03,0.1,-1.03,0.1,-0.14,0.1,-1.03,0.1,-1.03,0.1]])

glassSouth1 = GRID([[0.25],[-1.6,-0.1,1.075,-0.1,1.175,-0.1,1.175,-0.1,1.075],[-2.5,-0.14,-0.1,1.03,-0.1,1.03,-0.1,-0.14,-0.1,1.03,-0.1,1.03]])

framesSouth = COLOR([0,0,0])(STRUCT([frameVertSouth1,frameHorizSouth1]))

glassSouth = COLOR([0,110,110,0.2])(glassSouth1)

windowsSouth =  T([2])([6.6])(R([1,2])(PI)(STRUCT([framesSouth,glassSouth,railingSouth])))

//Railing east

railingEast = COLOR([0,0,0])(GRID([[0.05],[-0.25,5.25],[-8.3,0.05,-0.25,0.05]]))

//Windows east

frameVertEast1 = GRID([[0.25],[-5.75,0.1,-1.1,0.1,-1.1,0.1],[-2.5,-1.2,1.2,-1.3,1.2,-1.3,1.2]])

frameHorizEast1 = GRID([[0.25],[-5.75,2.5],[-2.5,-1.1,0.1,-1.2,0.1,-1.1,0.1,-1.2,0.1,-1.1,0.1,-1.2,0.1]])

glassEast1 = GRID([[0.25],[-5.75,-0.1,1.1,-0.1,1.1],[-2.5,-1.1,-0.1,1.2,-0.1,-1.1,-0.1,1.2,-0.1,-1.1,-0.1,1.2]])

glassEast = COLOR([0,110,110,0.2])(glassEast1)

framesEast = COLOR([0,0,0])(STRUCT([frameHorizEast1,frameVertEast1]))

windowsEast = R([1,2])(3*PI/2)(STRUCT([framesEast,glassEast,railingEast]))


//Windows west

frameVertWest1 = GRID([[0.25],[-2.7,0.05,-0.4,0.05],[-2,-0.05,0.4]])

frameHorizWest1 = GRID([[0.25],[-2.7,0.5],[-2,0.05,-0.4,0.05]])

glassWest1 = GRID([[0.25],[-2.7,-0.05,0.4],[-2,-0.05,0.4]])

frameVertWest2 = GRID([[0.25],[-1.7,0.1,-1.1,0.1,-1.1,0.1],[-4,-0.1,0.8]])

frameHorizWest2 = GRID([[0.25],[-1.7,2.5],[-4,0.1,-0.8,0.1]])

glassWest2 = GRID([[0.25],[-1.8,1.1,-0.1,1.1],[-4,-0.1,0.8]])

frameVertWest3 = GRID([[0.25],[-1,0.05,-0.1,0.05,-0.3,0.05,-0.1,0.05],[-6.5,-0.05,0.9]])

frameHorizWest3 = GRID([[0.25],[-1,0.2,-0.3,0.2],[-6.5,0.05,0.9,0.05]])

glassWest3 = GRID([[0.25],[-1,-0.05,0.1,-0.05,-0.3,-0.05,0.1],[-6.5,-0.05,0.9]])

framesWest = COLOR([0,0,0])(STRUCT([frameVertWest1,frameHorizWest1,frameVertWest2,frameHorizWest2,frameVertWest3,frameHorizWest3]))

glassWest  = COLOR([0,110,110,0.2])(STRUCT([glassWest1,glassWest2,glassWest3]))

windowsWest = T([1,2])([11,6.6])(R([1,2])(PI/2)(STRUCT([framesWest,glassWest])))

building4 = STRUCT([windowsWest,windowsEast,windowsSouth,windowsNorth])