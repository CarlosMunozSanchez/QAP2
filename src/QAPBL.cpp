/* 
 * File:   QAPBL.cpp
 * Author: carlos
 * 
 * Created on 5 de abril de 2023, 13:18
 * Rework on 19 de mayo de 2023, 13:00
 */

#include "QAPBL.h"
#include "random.hpp"
#include "funciones.h"
#include <iostream>

using namespace std;
// get base random alias which is auto seeded and has static API and internal state
// it is not threads secure, you can also use ::random_thread_local
using Random = effolkronium::random_static;


QAPBL::QAPBL(const vector<vector<int>> & flujos, const vector<vector<int>> & distancias, int seed) {
    //inicializar la máscara dlb
    dlb.assign(flujos.size(), 0);
    
    solucion.resize(dlb.size());
    
    //inicializar el vector solucion
    for(int i = 0; i < solucion.size(); i++){
        solucion[i] = i;
    }    
    //inicialización aleatoria
    Random::seed(seed);   
    Random::shuffle(solucion);

    float aux = 1;
    coste = evaluarSolucion(solucion, flujos, distancias, aux);
        
    busquedaLocal(flujos, distancias);
}


QAPBL::QAPBL(const std::vector<int> & inicial, const vector<vector<int>> & flujos, 
        const vector<vector<int>> & distancias, int seed) {
    //inicializar la máscara dlb
    dlb.assign(flujos.size(), 0);
    
    //inicializar el vector solucion
    solucion = inicial;

    float aux = 1;
    coste = evaluarSolucion(solucion, flujos, distancias, aux);
        
    busquedaLocal(flujos, distancias);
}

int QAPBL::comprobarMovimiento(int i, int j, const vector<vector<int>> & flujos, const vector<vector<int>> & distancias){
    int diferencia = 0;
    
    //factorización del cambio en el coste debido a cambio
    for(int k = 0; k < solucion.size(); k++){
        if(k != i and k != j){
            diferencia +=
            flujos[i][k] * (distancias[solucion[j]][solucion[k]] - distancias[solucion[i]][solucion[k]]) +
            flujos[j][k] * (distancias[solucion[i]][solucion[k]] - distancias[solucion[j]][solucion[k]]) +
            flujos[k][i] * (distancias[solucion[k]][solucion[j]] - distancias[solucion[k]][solucion[i]]) +
            flujos[k][j] * (distancias[solucion[k]][solucion[i]] - distancias[solucion[k]][solucion[j]]);
        }
    }
    
    return diferencia;
}
        
void QAPBL::aplicarMovimiento(int i, int j){
    int aux = solucion[i];
    //permutar i y j
    solucion[i] = solucion[j];
    solucion[j] = aux;    
}

void QAPBL::busquedaLocal(const vector<vector<int>> & flujos, const vector<vector<int>> & distancias){
    const int MAX_EVAL = 400;
    int iter = 0;
    bool hay_mejora;
   
    //Para llevar la cuenta de cuántos bits de dlb están a 0
    int dlbStatus = dlb.size();
    
    //Se explora el espacio de soluciones mientras se haya al menos un bit a 0
    //en dlb y no se supere el umbral MAX_EVAL
    while(iter < MAX_EVAL and dlbStatus > 0){
        
        //para cada unidad
        for(int i = 0; i < solucion.size(); i++){
            //se comprueba la máscara
            if(!dlb[i]){
                hay_mejora = false;
                //se busca la primera solución que mejore
                for(int j = 0; j < solucion.size(); j++){
                    int aux = comprobarMovimiento(i, j, flujos, distancias);
                    iter++;
                    //cuando encontremos mejora, la realizamos
                    if(aux < 0){
                        hay_mejora = true;
                        aplicarMovimiento(i, j);
                        
                        //actualizamos los bits
                        if(dlb[i] != false){
                            dlb[i] = false;
                            dlbStatus++;
                        }
                        if(dlb[j] != false){
                            dlb[j] = false;
                            dlbStatus++;
                        }
                    }
                }
                //si no hay mejora
                if(!hay_mejora){
                    dlb[i] = true;
                    dlbStatus--;
                }
            }            
        }
    }

    float aux = 1;
    coste = evaluarSolucion(solucion, flujos, distancias, aux);
}
