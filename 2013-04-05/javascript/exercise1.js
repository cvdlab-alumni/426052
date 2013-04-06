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
function (dims) {
	dims = dims.map(function(dim){
		return dim - 1;
	});
    return function (values) {
      return function (object) {
        return object.clone().scale(dims, values);
      };
    };
  }

S3=S2
S2=S1
S1=S0

function (array) {
      return array[index];
    }

GRID = SIMPLEX_GRID
VIEW = DRAW
NN = REPLICA

function (n) {
    return function (value) { 
      var result = [];
      var i;
      for (i = 0; i < n; i++) {
        result = result.concat(value);
      }
      return result;
    };
  }

//Funzione che definisce un arco di circonferenza
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


//Exercise1

//Ground floor pillars
pillars01 = GRID([[0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25],[-6.6,0.25],[-0.04,2.46]]);
pillars02 = GRID([[-0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25],[-5.25,0.25],[-0.04,2.46]]);
column_base = arc(2*PI,0,0.125);
column03 = EXTRUDE([2.36])(column_base);
columns0 = STRUCT(CONS(AA(T([1,2,3]))([[0.125,0.125,0.04],[2.875,0.125,0.04],
[5.625,0.125,0.04],[8.375,0.125,0.04],[11.125,0.125,0.04],[0.125,5.375,0.04]]))(column03));
pillars0 = STRUCT([pillars01,pillars02,columns0]);

//First floor pillars
pillars11 = GRID([[0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25],[-6.6,0.25],[-0.14,-2.36,-0.14,2.36]]);
pillars12 = GRID([[0.25,-2.5,0.25,-2.5,0.25,-2.5,-0.25,-2.5,0.25],[-5.25,0.25],[-0.14,-2.36,-0.14,2.36]]);
pillars13 = GRID([[0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25],[0.25],[-0.14,-2.36,-0.14,2.36]]);
column13 = EXTRUDE([2.36])(column_base);
column1 = T([1,2,3])([8.375,5.375,2.64])(column13);
pillars1 = STRUCT([pillars11,pillars12,pillars13,column1]);

//Second floor pillars
pillars2 = GRID([[0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25],[0.25,-5,0.25,-1.1,0.25],[-5.14,2.36]]);

//Third floor pillars
pillars312 = GRID([[0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25],[-5.25,0.25,-1.1,0.25],[-7.64,2.36]]);
pillars33 = GRID([[0.25,-2.5,-0.25,-2.5,0.25,-2.5,0.25,-2.5,0.25],[0.25],[-7.64,2.36]]);
pillars3 = STRUCT([pillars312,pillars33]);

//building model
building1 = STRUCT([pillars0,pillars1,pillars2,pillars3]);