T = function (dims){
	var dims = dims.map(function(dim){
		return dim-1;
	})
	return function(values){
		return function(object) {
			return object.clone().translate(dims, values);
		};
	};
}

var pointfloor = [[0,0],[0,6.85],[11.25,6.85],[11.25,0]]
var cellsfloor = [[0,1,3],[3,2,1]];
var floor = SIMPLICIAL_COMPLEX(pointfloor)(cellsfloor);
var floor3D = EXTRUDE([0.14])(floor);
DRAW(floor);

floor00 = EXTRUDE([0.04])(floor3D);

var pointsfloor0 = [[0,5.25],[0,6.85],[1.5,6.85],[8.65,6.85],[8.825,6.85],[8.825,4],[8.65,4],[2.75,2.15],[8.65,2.15],
			[2.75,1.775],[1.5,1.775],[1.5,5.25]];
var cellsfloor0 = [[0,1,2],[0,11,2],[2,3,6],[2,6,7],[6,8,7],[7,9,10],[7,10,11],[2,10,7],[3,4,6],[4,6,5]];
var floor0 = SIMPLICIAL_COMPLEX(pointsfloor0)(cellsfloor0);
var floor01 = T([3])([0.04])(EXTRUDE([0.1])(floor0));

floor1 = T([3])([2.5])(floor3D)

floor2 = T([3])([5])(floor3D)

floor3 = T([3])([7.5])(floor3D)

floor4 = T([3])([10])(floor3D)

DRAW(STRUCT([floor00,floor01,floor1,floor2,floor3,floor4]))