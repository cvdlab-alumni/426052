
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

var model = mountains;

//DRAW(model);