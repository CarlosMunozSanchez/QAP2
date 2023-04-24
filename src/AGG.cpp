/* 
 * File:   AGG.cpp
 * Author: carlos
 * 
 * Created on 17 de abril de 2023, 19:23
 */

#include "AGG.h"
#include "random.hpp"
// get base random alias which is auto seeded and has static API and internal state
// it is not threads secure, you can also use ::random_thread_local
using Random = effolkronium::random_static;
using namespace std;

AGG::AGG(int n_genes, int tipo, const vector<vector<int>> & flujos, 
        const vector<vector<int>> & distancias, int seed) {
    
    this->tipo = tipo;
    
    //inicialización aleatoria
    Random::seed(seed); 
    
    //generar población inicial
    //me genero un vector auxiliar
    vector<int> aux;
    for(int i = 0; i < n_genes; i++){
        aux.push_back(i);
    }
    
    //la poblacion inicial es una permutacion de este vector
    for(int i = 0; i < POBLACION_SIZE; i++){
        poblacion[i] = Random::shuffle(aux);        
        
        //calculo el coste de cada solución para saber cuál es la mejor
        /*float f = 1;
        int c = evaluarSolucion(poblacion[i], flujos, distancias, f);
        if(c < coste_mejor){
            coste_mejor = c;
            mejor_solucion = i;
        }*/
    }
    
}
    
vector<int> AGG::seleccion(vector<int> & c1, vector<int> & c2);
    
pair<vector<int>, vector<int>> AGG::cruce(vector<int> & c1, 
    vector<int> & c2, int tipo);

void AGG::mutacion(vector<int> & cromosoma);
