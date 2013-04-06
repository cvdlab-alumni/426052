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

//Exercise3

//########## NORTH ##########################################################################################

north0 = GRID([[0.25],[0.25,-6.35,0.25],[-2.50,7.5]])
north1 = GRID([[0.25],[-0.25,5],[-2.5,0.85,-1.5,-0.15,0.85,-1.5,-0.15,0.85]])
north2 = GRID([[0.25],[-0.25,5],[-2.5,-0.85,-1.5,0.15,-0.85,-1.5,0.15,-0.85,-1.5,0.15]])
north3 = GRID([[0.25],[-5.25,1.1],[-2.5,7.5]])
north4 = GRID([[0.25],[-6.35,0.25],[-2.5,0.15,-2.3,0.15,-2.3,0.15,-2.3,0.15]])
north5 = GRID([[0.25],[6.85],[-10,0.8]])

north = T([1])([11])(STRUCT([north0,north1,north2,north3,north4,north5]))

//######### SOUTH #############################################################################################

south0 = GRID([[0.25],[6.85],[-2.5,0.14,-2.36,0.14]])
south1 = GRID([[0.25],[6.85],[-7.5,0.8,-1.7,0.8]])
south2 = GRID([[0.25],[-1.6,5],[-2.5,-0.14,-1.13,0.1,-1.13,-0.14,-1.13,0.1]])
south3 = GRID([[0.25],[-2.775,0.1,-1.175,0.1,-1.175,0.1,-1.175,0.1],[-2.5,-0.14,2.36,-0.14,2.36]])
south4 = GRID([[0.25],[0.25,-1.1,0.25],[-2.5,2.5]])
south5 = GRID([[0.25],[1.6],[-5,3.3]])
south6 = GRID([[0.25],[-6.6,0.25],[-2.5,-0.14,2.36,-0.14,2.36]])

south = T([2])([6.6])(R([1,2])(PI)(STRUCT([south0,south1,south2,south3,south4,south5,south6])))

//######## EAST ################################################################################################

east0 = GRID([[0.25],[-5.5,0.25],[-2.5,-1.1,1.4,-1.1,1.4,-1.1,1.4]])
east1 = GRID([[0.25],[-8.25,3],[-2.5,-1.1,1.4,-1.1,1.4,-1.1,1.4]])
east2 = GRID([[0.25],[-5.5,5.75],[-2.5,1.1,-1.4,1.1,-1.4,1.1]])
east3 = GRID([[0.25],[5.5],[-2.5,5.8]])
east4 = GRID([[0.25],[11.25],[-10,0.8]])

east = R([1,2])(3*PI/2)(STRUCT([east0,east1,east2,east3,east4]))

//###### WEST ##################################################################################################

west0 = GRID([[0.25],[-4.2,7.05],[10.8]])
west1 = GRID([[0.25],[4.2],[-7.5,3.3]])
west2 = GRID([[0.25],[1],[-2.5,5]])
west3 = GRID([[0.25],[-1,0.7],[-2.5,4]])
west4 = GRID([[0.25],[-1.2,0.3],[-6.5,1]])
west5 = GRID([[0.25],[-1.7,2.5],[-5,2.5]])
west6 = GRID([[0.25],[-2.9,0.1],[-4,1]])
west7 = GRID([[0.25],[-1.7,2.5],[-2.5,1.5]])
west8 = GRID([[0.25],[-1.1,1.3,-0.5,1.3],[2.5]])
west9 = GRID([[0.25],[-1.1,-1.3,0.5],[2]])

west = T([1,2])([11,6.6])(R([1,2])(PI/2)(STRUCT([west9,west8,west7,west6,west5,west4,west3,west2,west1,west0])))

//########## BUILDING MODEL #####################################################################################################

building3 = STRUCT([north,south,east,west])



