var dom1D = INTERVALS(1)(32)
var dom2D = PROD1x1([dom1D,dom1D])

function nodes (points) { 
  var m = points.length;
  var k = 2;
  var n = (m + k + 1);
  var l = n - 3;
  var j = 1;
  var knots = [];
  for (var i = 0; i < 3; i++) {
    knots[i] = 0;
  };
  for (var i = 3; i < l; i++, j++) {
    knots[i] = j;
  };
  for (var i = l; i < n; i++) {
    knots[i] = j;
  };
  return knots;
};

function trasla(punti,s){
  punti = punti.map(function(p){return [p[0]+s[0],p[1]+s[1],p[2]+s[2]]});
  return punti;
}


//funzione che trasla alcuni punti della curva in funzione del punto in cui si trova la curva

//points : array di array con i punti corrispondenti alle 15 barre verticali
//k1 e k2 : da quale punto a quale punto modificare ogni curva
//t : parametro che viene incrementato per modificare in modo parametrico ogni coppia di curve
//d : coordinata da modificare 

function trasla_param(points,k1,k2,t,d){ 
  var p = points;
  var j = t;
  for(i=2; i<16; i++){
    for(k=k1;k<k2;k++){
      if (i==14){
        p[i][k][d] = p[i-2][k][d];
        p[i+1][k][d] = p[i-2+1][k][d];
      }
      else {
      p[i][k][d] = p[i][k][d]+ j;
      p[i+1][k][d] = p[i+1][k][d]+ j;
      p[points.length-2-i][k][d] = p[points.length-2-i][k][d]+ j;
      p[points.length-2-i+1][k][d] = p[points.length-2-i+1][k][d]+ j;
      }
    }
    j = j + t;
    i++;
  }
  return p;
};


function trasla_z(points,t){  //trasla di una quantità fissa 
  var p = points;
  for(i=2; i<28; i++){
    p[i][0][2] = p[i][0][2] + t;
  }
  return p;
};

function trasla_y(points){ //trasla di una quantità fissa
  var p = [];
  for(i=0; i<30; i++){
    p[i] = trasla(points,[0,i*(0.227586207),0]);
  }
  return p;
};

/////////////////////////////////////////////////// 

function horizontal(p1,p2){
  var cp1 = BEZIER(S0)(p1);
  var p11 = trasla(p1,[-0.25,0,0]);
  var cp11 = BEZIER(S0)(p11);
  var p_0 = MAP(BEZIER(S1)([cp1,cp11]))(dom2D);
  var cp2 = BEZIER(S0)(p2);
  var p22 = trasla(p2,[-0.25,0,0]);
  var cp22 = BEZIER(S0)(p22);
  var p_1 = MAP(BEZIER(S1)([cp2,cp22]))(dom2D);
  var p_2 = MAP(BEZIER(S1)([cp1,cp2]))(dom2D);
  var p_3 = MAP(BEZIER(S1)([cp11,cp22]))(dom2D);
  var p = STRUCT([p_0,p_1,p_2,p_3]);
  return p;
};


function repeat_horizontal(f){
  var fig = f;
      for(i=1; i<15; i++){
        f1 = T([1])([i*(0.454)])(f);
        fig = STRUCT([fig,f1]);
    }
  return fig;
};

function vertical(p){
  var v = vett_traslX();                    //calcolo del vettore per la base circolare
  var p0 = trasla_y(p);                     //genero le curve tutte alla stessa altezza a partire da quella data
  var p1 = trasla_x0(p0,v);                 //traslo le x dei primi punti di ogni curva per ottenere la base circolare
  var p2 = trasla_param(p1,1,2,0.15,2);     //modifico le z di alcuni punti in modo parametrico 
  var p3 = trasla_param(p2,3,4,-0.05,0);    //modifico le x di alcuni punti in modo parametrico 
  var p4 = trasla_param(p3,6,9,-0.05,0);    //modifico le x di alcuni punti in modo parametrico 
  var p5 = trasla_param(p4,6,9,-0.1,2);     //modifico le z di alcuni punti in modo parametrico 
  var p6 = trasla_z(p5,0.25);               //modifico le z dei primi punti di ogni curva di una quantità prefissata 
  return p6;
};

function depth(p0,p1){
  var a0 = vertical(p0);
  var a1 = vertical(p1);
  var p = [];
  for(i=0; i<30; i++){
    knots0 = nodes(a0[i]); 
    knots1 = nodes(a1[i]);
    cp0 = NUBS(S0)(2)(knots0)(a0[i]);
    cp1 = NUBS(S0)(2)(knots1)(a1[i]);
    p[i] = BEZIER(S1)([cp0,cp1]);
  }
  return p;
};

//////////////////////////////////////////////////

//funzione che genera il vettore traslazione da applicare alle curve per ottenere la base circolare 

function vett_traslX(){   
 var p = [-3.1,-2.6572,-2.2144,-1.7716,-1.3288,-0.886,-0.4432,0,0.4432,0.886,1.3288,1.7716,2.2144,2.6572,3.1];
 var v = [];
 v[0] = 0;
 v[14] = 0;
 for (i=1;i<p.length-1;i++){
  v[i] = Math.sqrt(9.61 - (p[i]*p[i]));
  }
 return v;
};

// funzione che utilizza il vettore generato da vett_traslX e trasla le curve 

function trasla_x0(points,v){
  var p = points;
  var j = 1;
  for(i=2; i<27; i++){
    p[i][0][0] = p[i][0][0] - v[j] + 0.12;
    p[i+1][0][0] = p[i+1][0][0] - v[j] + 0.12;
    j++;
    i++;
  }
  return p;
};


/////////////////////////////////////////////////

function width(points){
  var p = [];
  var j = 0;
  for(i=0; i<30; i++){
    knots0 = nodes(points[i]); 
    knots1 = nodes(points[i+1]);
    cp0 = NUBS(S0)(2)(knots0)(points[i]);
    cp1 = NUBS(S0)(2)(knots1)(points[i+1]);
    p[j] = BEZIER(S1)([cp0,cp1]);
    i++;
    j++;
  }
  return p;
};


function gen_map(p){
  var vertical = MAP(p[0])(dom2D);
  for(i=1;i<15;i++){
    v = MAP(p[i])(dom2D);
    vertical = STRUCT([vertical,v]);
  }
  return vertical;
};



function map_spessore(p){
  var s = MAP(p[0])(dom2D);
  for(i=1;i<30;i++){
    v = MAP(p[i])(dom2D);
    s = STRUCT([s,v]);
  }
  return s;
};

////////////////////////////////////////////////////

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



////////////////////////////////////////////////////////////////////////////////////////////////////////////

black_0 = [0.35,0.35,0.35,1]
black_1 = [0.2,0.2,0.2,1]
brown_0 = [0.8,0.521,0.247,1]
brown_1 = [0.7,0.421,0.147,1]

a0 = [[4.8,-3.3,0],[6,-3.3,1.95],[6.35,-3.3,3.35],[6.05,-3.3,4.25],[4.85,-3.3,4.85],[3.6,-3.3,4.85],[2.55,-3.3,4.8],[1.7,-3.3,5.2],[1.1,-3.3,6.3],[0.7,-3.3,7],[0.5,-3.3,8]]

a1 = [[4.4,-3.3,0],[5.6,-3.3,1.95],[5.95,-3.3,3.35],[5.7,-3.3,4.1],[4.85,-3.3,4.45],[3.6,-3.3,4.45],[2.55,-3.3,4.4],[1.25,-3.3,5.05],[0.7,-3.3,6.3],[0.3,-3.3,7],[0.1,-3.3,8]]

p0_0 = vertical(a0)
p0_1 = width(p0_0)
p0 = COLOR(brown_1)(gen_map(p0_1))

p1_0 = vertical(a1)
p1_1 = width(p1_0)
p1 = COLOR(brown_1)(gen_map(p1_1))

p2_0 = depth(a0,a1)
p2 = COLOR(brown_0)(map_spessore(p2_0))

p3_0 = T([0,1,2])([0.1,-3.3,8])(CUBOID([0.4,0.227586207,0.05]))
p3 = COLOR(brown_1)(repeat_horizontal(p3_0))

vertical_0 = STRUCT([p0,p1,p2,p3])

a2 = [[4.9,-3.3,4.45],[4.9,-3.2,4.05],[4.9,3.2,4.05],[4.9,3.3,4.45]]

a3 = [[4.9,-3.3,4.45],[4.9,-3.2,3.85],[4.9,3.2,3.85],[4.9,3.3,4.45]]

p4 = COLOR(brown_1)(horizontal(a2,a3))

p5_0 = T([0,1,2])([0.1,-3.3,7.772])(R([0,2])([-PI/20])(CUBOID([0.1,6.6,0.227])))
p5_1 = T([0,2])([0.69,-0.454])(R([0,2])([-PI/40])(p5_0))
p5_2 = T([0,2])([0.49,-0.452])(R([0,2])([-PI/60])(p5_1))
p5 = COLOR(brown_1)(T([2])([0.04])(STRUCT([p5_0,p5_1,p5_2])))

p6_0 = EXTRUDE([0.4])(arc(PI,3.073,3.3))
p6 = COLOR(brown_0)(T([0])([4.65])(R([0,1])([PI/2])(p6_0)))

horizontal_0 = STRUCT([p4,p5,p6])

pantosh_chair = STRUCT([vertical_0,horizontal_0])

//DRAW(pantosh_chair)