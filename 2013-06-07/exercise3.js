
///////// EXERCISE 1 ////////////////////////////////////////////////////////////////////////////////////////////////


var c_mountains = [0.9568,0.64313,0.3764,1];

var dom = DOMAIN([[0,10],[0,7]])([24,12]);
v_x = [];
v_y = [];
v_z = [];

var x = function (v) {
	a = v[0];
	v_x.push(a);
	return a;
	};

var y = function (v) {
	b = v[1];
	v_y.push(b);
	return b;
	};

var z = function (v) {
	if (SIN(v[1])+COS(v[0])>0.05 && v[0]<3 && v[0]>0.01 && v[1]>0.01 && v[1]<6.99){
		k = 1 +Math.random();
		v_z.push(k);
		return k;}
	if (SIN(v[1])+COS(v[0])<0.2 && v[0]>9 && v[0]<9.99 && v[1]>0.01 && v[1]<6.99){
		k = 1.9 +Math.random();
		v_z.push(k);
		return k;}
	if (SIN(v[1])+COS(v[0])<0.2 && v[0]>6 && v[0]<9.99 && v[1]>0.01 && v[1]<6.99){
		k = 1.3 +Math.random();
		v_z.push(k);
		return k;}
	else{
		k = 1;
		v_z.push(k);
		return k;}
	};

var base0 = CUBOID([10,7,0.01]);
var base1 = CUBOID([10,0.01,1]);
var base2 = T([1])([6.99])(base1);
var base3 = CUBOID([0.01,7,1]);
var base4 = T([0])([9.99])(base3);

var base = STRUCT([base0,base1,base2,base3,base4]);

var mappings = [x,y,z];

var mountains = COLOR(c_mountains)(STRUCT([MAP(mappings)(dom),base]));



///////// EXERCISE 2 //////////////////////////////////////////////////////////////////////////////////////////////////


var c_lake = [0,0.8,0.6,0.7];

function lake(p0,p1){
	dom2D = DOMAIN([[0,1],[0,1]])([30,50]);
	c0 = BEZIER(S0)(p0);
	c1 = BEZIER(S0)(p1);
	c_m = MAP(BEZIER(S1)([c0,c1]))(dom2D);
	c = COLOR(c_lake)(c_m)
	return c;
}

var l0 = [[5.1,7,1.01],[5.2,6.7,1.01],[5.4,6.4,1.01],[5.6,6.5,1.01],[5.8,6.2,1.01],[5.9,5.5,1.01]];
var l1 = [[7.5,7,1.01],[7.5,5.5,1.01]];

var lake1 = lake(l0,l1);

var l2 = [[5,0,1.01],[4.9,1.5,1.01],[3.4,2,1.01],[2,2.3,1.01],[2.7,3.7,1.01],[2.6,4.2,1.01],[0,5.2,1.01]];
var l3 = [[0,0,1.01],[0,5,1.01]];

var lake2 = lake(l2,l3);

var lakes = STRUCT([lake1,lake2]);



/////// EXERCISE 3 ///////////////////////////////////////////////////////////////////////////////////////////////


var c_tree = [0.133,0.545,0.133,1];
var c_trunk = [0.5882,0.2941,0,1];

function cylinder(r,h){
	d = DISK(r)(36)
	c = EXTRUDE([h])(d)
	return c
}

f = function(c,px_1,px_2,py_1,py_2){
	s = "";
	for (i=0;i<v_x.length;i++){
	if(v_x[i]>=px_1 && v_x[i]<=px_2 && v_y[i]>=py_1 && v_y[i]<=py_2)
		s = STRUCT([T([0,1,2])([v_x[i],v_y[i],v_z[i]])(c),s])
	}
	return s;
}

function cone(p){
	domain = PROD1x1([INTERVALS(1)(20),INTERVALS(1)(6)]);
	apex = [0,0,0.35];
	profile = CUBIC_HERMITE(S0)(p);
	c0 = MAP(CONICAL_SURFACE(apex)(profile))(domain);
	c1 = R([0,1])([PI])(c0);
	c = COLOR(c_tree)(STRUCT([c0,c1]));
	return c;
}


var t0 = [[-0.1,0,0],[0.1,0,0],[0,0.4,0],[0,-0.4,0]];
var t1 = COLOR(c_trunk)(cylinder(0.025,0.3));

var tree = STRUCT([T([2])([0.12])(cone(t0)),t1]);

var forest1 = f(tree,0,2,0.5,3);
var forest2 = f(tree,7,9.9,3,6);
var forest3 = f(tree,9,9.99,0.2,3);

var forests = STRUCT([forest1,forest2,forest3]);

var model = STRUCT([mountains,lakes,forests]);

//DRAW(model);
