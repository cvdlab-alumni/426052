
//////// EXERCISE 6 ////////////////////////////////////////////////////////////////////////////////////////////////


function lar_to_obj(model){
	s = "";
	v = model[0];
	fv = model[1];
	for(i=0;i<v.length;i++){
		s += "v\t" + v[i][0]+" "+v[i][1]+" 0\n";
	} 
	for(j=0;j<fv.length;j++){
		s += "f\t";
		for(k=0;k<fv[j].length;k++){
			s += fv[j][k]+" "; 
		}
	s += "\n";
	}
	return s;
}

// Per testare la funzione ho utilizzato gli stessi dati dell'esempio fatto in aula dal professore


V =[[0,6],[0,0],[3,0],[6,0],[0,3],[3,3],[6,3],[6,6],[3,6]];


FV = [[5,6,7,8],[0,5,8],[0,4,5],[1,2,4,5],[2,3,5,6],[0,8,7],[3,6,7],[1,2,3],[0,4,1]];

m = [V,FV];

a = lar_to_obj(m);