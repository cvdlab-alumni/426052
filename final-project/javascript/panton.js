var dom1D = INTERVALS(1)(32)
var dom2D = PROD1x1([dom1D,dom1D])

function nodes (points) { //funzione che mi calcola i nodi a partire dai punti
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



function trasla_fix(points,k1,k2,t,d){ 
  var p = points;
  for(i=1; i<12; i++){
    for(k=k1;k<k2;k++){
      p[i][k][d] = p[i][k][d] + t;
    }
  }
  return p;
};

function trasla_y(points){ //trasla di una quantità fissa
  var p = [];
  for(i=0; i<13; i++){
    p[i] = trasla(points,[0,i*(0.333333333),0]);
  }
  return p;
};


function trasla_x0(points,v){
  var p = points;
  for(i=1; i<12; i++){
    p[i][0][0] = p[i][0][0] - v[i];
  }
  return p;
};

function trasla_z12(points,v){
  var p = points;
  for(i=1; i<12; i++){
    p[i][12][2] = p[i][12][2] + v[i];
    p[i][13][2] = p[i][13][2] + v[i];
  }
  return p;
};



function panton(p){
  var v0 = [0,2.3,3.1,3.55,3.8,4,4.1,4,3.8,3.55,3.1,2.3,0];           //vettore per la base
  var v1 = [0,0.5,0.6,0.7,0.8,0.85,0.85,0.85,0.8,0.7,0.6,0.5,0];      //vettore per la parte alta
  var p0 = trasla_y(p);                  //genero le curve tutte alla stessa y a partire da quella data
  var p1 = trasla_x0(p0,v0);             //traslo le x dei primi punti di ogni curva per ottenere la base circolare
  var p2 = trasla_fix(p1,1,2,-2,0);      //modifico le x e le z di alcuni punti  
  var p3 = trasla_fix(p2,2,3,-0.5,0);     
  var p4 = trasla_fix(p3,2,4,0.5,2);      
  var p5 = trasla_fix(p4,7,10,-0.25,2);     
  var p6 = trasla_z12(p5,v1);           // traslo le z degli ultimi punti per ottenere la parte alta tondeggiante
  var p7 = trasla_fix(p6,12,13,-0.15,2);
  var p8 = trasla_fix(p7,9,11,-0.3,0);          
  return p8;
};



/////////////////////////////////////////////////

function gen_nubs(points){
  var p = [];
  for(i=0; i<13; i++){
    knots0 = nodes(points[i]); 
    p[i] = NUBS(S0)(2)(knots0)(points[i]);
  }
  return p;
};

function gen_map(p){
  var vertical = MAP(p[0])(dom1D);
  for(i=1;i<13;i++){
    v = MAP(p[i])(dom1D);
    vertical = STRUCT([vertical,v]);
  }
  return vertical;
};


/////////////////////////////////////////////////////////////////////////////////////////

c = [0,0.5,1,1]
g = [1,0.84,0,1]
r = [1,0,0,1]

p_in = [[5.05,-2.3,0],[5.55,-2.3,1.05],[5.9,-2.3,2.2],[6,-2.3,3.15],[5.65,-2.3,4.05],[4.95,-2.3,4.525],[3.65,-2.3,4.5],[2.4,-2.3,4.3],[1.85,-2.3,4.55],[1.45,-2.3,5.2],[1.15,-2.3,6.4],[0.6,-2.1,6.9]]
knotspin = nodes(p_in) 
cpin = NUBS(S0)(2)(knotspin)(p_in)

p_fin = [[5.05,2.3,0],[5.55,2.3,1.05],[5.9,2.3,2.2],[6,2.3,3.15],[5.65,2.3,4.05],[4.95,2.3,4.525],[3.65,2.3,4.5],[2.4,2.3,4.3],[1.85,2.3,4.55],[1.45,2.3,5.2],[1.15,2.3,6.4],[0.6,2.1,7]]
knotspfin = nodes(p_fin) 
cpfin = NUBS(S0)(2)(knotspfin)(p_fin)

p0 = [[4.75,-2,0],[5.25,-2,1.05],[5.6,-2,2.2],[5.7,-2,3.15],[5.35,-2,3.9],[4.95,-2,4.275],[3.65,-2,4.4],[2.4,-2,4.3],[1.85,-2,4.55],[1.45,-2,5.2],[1.15,-2,6.4],[0.9,-2,7.3],[0.55,-2,7.5],[0.45,-2,7.1]]

p_0 = panton(p0)
p_1 = gen_nubs(p_0)
p_2 = BEZIER(S1)([cpin,p_1[0],p_1[1],p_1[2],p_1[3],p_1[4],p_1[5],p_1[6],p_1[7],p_1[8],p_1[9],p_1[10],p_1[11],p_1[12],cpfin])

panton_0 = MAP(p_2)(dom2D)

p1 = [[4.9,-2.5,0],[5.35,-2.5,1.05],[5.65,-2.5,2.2],[5.55,-2.5,3.15],[5.25,-2.5,3.7],[4.95,-2.5,4.075],[3.65,-2.5,4.15],[2.4,-2.5,4.1],[1.85,-2.5,4.4],[1.15,-2.5,5.2],[0.8,-2.5,6.4],[0.6,-2.1,6.9]]
knotsp1 = nodes(p1) 
cp1 = NUBS(S0)(2)(knotsp1)(p1)

p2 = [[4.9,2.5,0],[5.35,2.5,1.05],[5.65,2.5,2.2],[5.55,2.5,3.15],[5.25,2.5,3.7],[4.95,2.5,4.075],[3.65,2.5,4.15],[2.4,2.5,4.1],[1.85,2.5,4.4],[1.15,2.5,5.2],[0.8,2.5,6.4],[0.6,2.1,6.9]]
knotsp2 = nodes(p2) 
cp2 = NUBS(S0)(2)(knotsp2)(p2)

p_3 = BEZIER(S1)([cpin,cp1])

panton_1 = MAP(p_3)(dom2D)

p_4 = BEZIER(S1)([cpfin,cp2])

panton_2 =  MAP(p_4)(dom2D)

panton_chair = COLOR(r)(STRUCT([panton_0,panton_1,panton_2]))


DRAW(panton_chair)