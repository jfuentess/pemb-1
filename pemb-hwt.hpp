
//Funcion para responder consulta de contencion, sobre esta funcion estoy midiendo tiempos
bool Inside(unsigned int x, unsigned int levelx, unsigned int y, unsigned int levely){
	unsigned int comprobar = goLevel(x,levelx,levely);
	if(comprobar == y)return true;
	else return false;
}

// Funcion para mover region a otro nivel
unsigned int goLevel(unsigned int x, unsigned int levels, unsigned int levelt){
	unsigned int regions = regionID(goUp(x,levels),0);
	unsigned int regiont = regionID(goDown(regions,levelt),levelt);
	return regiont;
}


unsigned int goDown(unsigned int x, unsigned int level){
    int posB = m_B_select1[0](x+1);
    if(plano[posB] >= level){
      	int retornar = rankwt2(posB,level);
		bool pab = (m_B[level][retornar] == 1) ? true: false; //Compruebo si parentesis es abierto o cerrado 
		if(pab)return retornar;
		else {
    		retornar = m_B_st[level].find_open(retornar); // obtengo el parentesis que abre
    		return retornar;
		   }
      	}
    else{
    	int sumar2 = rankwt2(posB,level);
    	int parentesiscercano =selectwt(level,sumar2+1); //Posicion siguiente parentesis en B
    	sumar2 = rankwt2(parentesiscercano,level);
    	int retornar = sumar2; // Posicion del parentesis siguiente en m_b
    	bool abierto = (m_B[level][retornar] == 1) ? true: false;  //Compruebo si es abierto o cerrado
    	if(abierto /*&& sumar3 != sumar4 */){
    		retornar = m_B_st[level].parent_t(retornar); //en caso de ser abierto obtengo el padre en m_B
    		return retornar; //retorno directamente
    	}
    	else {
    		retornar = m_B_st[level].find_open(retornar); // Caso de ser cerrado obtengo parentesis que abre
    		return retornar;		
    	}
    }
}



unsigned int goUp(unsigned int x, unsigned int level){
	int posicion = m_B_select1[level](x+1); // posicion en M_B donde esta el x-ecimo parentesis abierto
	posicion = selectwt(level,posicion+1);
	return posicion;
}

//Funcion select propio
unsigned int selectwt(unsigned int x,unsigned int posx){
	int pos = 0,vactual = 0,romper = 0;
	long long int mov;
	//Bloque de codigo para la busqueda en mi arbol resumen, comienzo de la raiz y me muevo hasta encontrar la hoja
	//correspondiente
	for(int i = 0; i < m_height; i++){
		mov = pos<<1;
		if(resumen[x][mov+1]+vactual >= posx){
			pos = mov+1;
		}
		else{
			vactual+= resumen[x][mov+1];
			pos = mov+2;
		}
	}

	int actual = vactual;
	int posactual =  (pos-posinicial)*m_bs;
	pos = posactual;
    unsigned int r = posactual+ m_bs,l = posactual,m,totalentro = 0;
    //Bloque en donde realizo la busqueda dentro de la hoja, utilizo busqueda binaria
    while(r > l){
        m = (l+r)/2;
        unsigned int mitad = rankwt2(m+1,x);
        if(mitad >= posx)r = m;
        else l = m+1;
    }
//Este bloque comentado es lo mismo que el anterior, realizo la busqueda en la hoja, pero esta vez con busqueda secuencial    
/*
int maximo = B[0].size();
for(int i = 0; i < m_bs && posactual < maximo; i++,posactual++){
	if(actual == posx)break;
	if(plano[posactual] >= x){
		actual++;
		pos = posactual;
	}
//	romper++;
}
*/
	return l;
}

//Funcion rank propia
unsigned int rankwt2(unsigned int posx, unsigned int x){
	double time;
	unsigned t0, t1,t2;
	int pos = 0,vactual = 0,vactual2 = 0,romper = 0,divisor = 2, total =pow(2,m_height)*m_bs; 
	//Bloque de codigo para la busqueda en mi arbol resumen, comienzo de la raiz y me muevo hasta encontrar la hoja
	//correspondiente
	for(int i = 0; i < m_height; i++){
        if(total/divisor +vactual2 >= posx){
            pos = (pos<<1)+1;
        }
        else{
        vactual+= resumen[x][2*pos+1];
        vactual2+= total/divisor;
        pos = (pos<<1)+2;
        }
		divisor= divisor << 1;
    }
	int actual = vactual;
	int posactual =  (pos-posinicial)*m_bs;
	pos = posactual;
	//Bloque en donde realizo la busqueda dentro de la hoja, utilizo busqueda secuencial
	for( /* int i = 0; i < m_bs && */ ; posactual < posx;/* i++,*/posactual++){
        if(plano[posactual] >= x){
            actual++;
            pos = posactual;
        }
	}
    return actual;
}