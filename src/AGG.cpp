/* 
 * File:   AGG.cpp
 * Author: carlos
 * 
 * Created on 17 de abril de 2023, 19:23
 */

#include "AGG.h"
#include "random.hpp"
#include "funciones.h"
#include <unordered_map>
#include <iostream>
// get base random alias which is auto seeded and has static API and internal state
// it is not threads secure, you can also use ::random_thread_local
using Random = effolkronium::random_static;
using namespace std;

AGG::AGG(int tipo, int n_genes, const vector<vector<int>> & flujos, 
        const vector<vector<int>> & distancias, int seed) {
    
    this->tipo = tipo;
    
    poblacion.resize(POBLACION_SIZE);
    fitnessPoblacion.resize(POBLACION_SIZE);
        
    //inicialización aleatoria
    Random::seed(seed); 
    
    //generar población inicial
    //me genero un vector auxiliar
    vector<int> aux;
    for(int i = 0; i < n_genes; i++){
        aux.push_back(i);
    }
    
    //la poblacion inicial está compuesta por 
    //permutaciones de este vector
    for(int i = 0; i < POBLACION_SIZE; i++){
        poblacion[i] = aux;
        Random::shuffle(poblacion[i]);        
    }
    
    simularEvolucion(flujos, distancias, tipo);
    
}
    
pair<vector<int>, vector<int>> AGG::cruce(const vector<int> & c1, 
    const vector<int> & c2, int tipo){
    
    pair<vector<int>, vector<int>> hijos;

    switch (tipo){    
        case 0:{ //cruce posición
            vector<int> aux;
            aux.resize(c1.size());
            vector<int> noIguales;
            
            //comparo los dos cromosomas
            for(int i = 0; i < c1.size(); i++){
                //si la asignación coincide, la guardo
                if(c1[i] == c2[i]){
                    aux[i] = c1[i];
                }
                //Marco con -1 si la posición si no coincide
                else{
                    aux[i] = -1;
                    noIguales.push_back(c1[i]);
                }
            }
            //los dos hijos tienen las posiciones cuya asignación es igual
            hijos.first = hijos.second = aux;
            
            //ahora, mezclo aleatoriamente los indices que no eran iguales
            aux = noIguales;
            Random::shuffle(aux);
            Random::shuffle(noIguales);
            
            //Último paso: asignar estos índices a las posiciones marcadas con -1
            int j = 0;
            for(int i = 0; i < hijos.first.size() and j < aux.size(); i++){
                if(hijos.first[i] == -1){
                    hijos.first[i] = aux[j];
                    hijos.second[i] = noIguales[j];
                    j++;
                }
            }    
        }
        break;
        
        case 1:{ //cruce PMX
            hijos.first.resize(c1.size());
            hijos.second.resize(c1.size());
            
            //tamaño de la subcadena central
            int n = c1.size() / 3;
            
            //Estructuras para almacenar las correspondencias
            unordered_map<int, int> corr1_2;
            unordered_map<int, int> corr2_1;            
            
            //copio la cadena central
            for(int i = n; i < 2*n; i++){
                hijos.first[i] = c2[i];
                hijos.second[i] = c1[i];
                
                //almaceno las correspondencias
                corr1_2.insert(pair<int, int>(c1[i], c2[i]));
                corr2_1.insert(pair<int, int>(c2[i], c1[i]));
            }
            
            //Asigno el resto de posiciones siguiendo las correspondencias
            for(int i = 0; i < c1.size(); i++){                
                //Asigno el elemento que indica el otro padre
                hijos.first[i] = c1[i];                
                //Si hay conflictos, los resuelvo con las correspondencias
                while(corr2_1.count(hijos.first[i]) > 0){
                    hijos.first[i] = corr2_1.find(hijos.first[i])->second;
                }
                
                hijos.second[i] = c2[i];                
                //Si hay conflictos, los resuelvo con las correspondencias
                while(corr1_2.count(hijos.second[i]) > 0){
                    hijos.second[i] = corr1_2.find(hijos.second[i])->second;
                }
                
                //me salto la subcadena central, que ya la tengo asignada
                if(i == n-1)
                    i += n;
            }
            
        }        
        break;
    }
    
    return hijos;
}

void AGG::mutacion(vector<int> & cromosoma){
    //Obtengo dos números aleatorios distintos
    int g1 = Random::get(0, (int)cromosoma.size()-1);
    int g2;
    
    do{
        g2 = Random::get(0, (int)cromosoma.size()-1);
    }while(g1 == g2);
    
    //Permuto estas asignaciones
    int aux = cromosoma[g1];
    cromosoma[g1] = cromosoma[g2];
    cromosoma[g2] = aux;

}

void AGG::simularEvolucion(const vector<vector<int> >& flujos, 
        const vector<std::vector<int> >& distancias, int tipo){
    
    //épocas
    int T = 0;
    // evaluaciones
    int evaluaciones = 0;
    while(evaluaciones < MAX_EVAL){
        float aux = 1;
        
        //1. Calcular el fitness de toda la población actual
        
        //Si es la primera época del algoritmo, lo calculo todo
        //en caso contrario, sólo para los nuevos individuos
        if(T == 0){
            //Recorro cada individuo
            for(int i = 0; i < POBLACION_SIZE; i++){
                //Calculo su fitness
                int c = evaluarSolucion(poblacion[i], flujos, distancias, aux);
                fitnessPoblacion[i] = c;
                evaluaciones++;
                
                //tengo que separar a los mejores padres del conjunto
                //Los primeros, entran directamente al grupo
                if(i < POBLACION_SIZE - ESPERANZA_CRUCES*2){
                    mejoresPadres.insert(i);
                    fitnessMejoresPadres.insert(pair<int, int>(c, i) );
                }
                //para los siguientes, tengo que comprobar si son mejores que los que ya hay
                else{
                    //otengo el padre con el peor fitness, que está al final del multiset
                    multimap<int, int>::iterator it = fitnessMejoresPadres.end();
                    it--;
                    
                    //si el individuo es mejor, lo inserto y borro el peor
                    if(c < (*it).first){
                        mejoresPadres.erase((*it).second);
                        mejoresPadres.insert(i);
                        fitnessMejoresPadres.erase(it);
                        fitnessMejoresPadres.insert(pair<int, int>(c, i));
                        
                    }                
                }                
            }
        }
        //si no es la primera época del algoritmo, tengo que actualizar todos los fitness
        //menos los de los mejores padres de la generación anterior, pero eso lo hago al final
            
        // 2. SELECCIÓN
        int listaPadres[ESPERANZA_CRUCES*2];
        
        for(int i = 0; i < ESPERANZA_CRUCES*2; i++){
            //selecciono dos padres y los someto a torneo binario
            int p1 = Random::get(0, POBLACION_SIZE-1);
            int p2 = Random::get(0, POBLACION_SIZE-1);
            
            //me quedo con el mejor
            listaPadres[i] = fitnessPoblacion[p1] > fitnessPoblacion[p2] ? p1 : p2; 
            
        }
        
        // 3. CRUCE
        //Aquí guardaré los hijos
        vector<vector<int>> hijos;
        hijos.resize(ESPERANZA_CRUCES*2);
        
        pair<vector<int>, vector<int>> combinacion;
        for(int i = 0; i < ESPERANZA_CRUCES; i++){
            //obtengo la combinación de cada par de padres
            //pair<vector<int>, vector<int>> combinacion = 
            //            cruce(poblacion[listaPadres[i*2]], poblacion[listaPadres[i*2+1]], tipo);
            combinacion = cruce(poblacion[listaPadres[i*2]], poblacion[listaPadres[i*2+1]], tipo);
            
            if(i < ESPERANZA_MUTACION){
                //para los primeros ESPERANZA_MUTACION padres, muto a sus hijos
                mutacion(combinacion.first);
                mutacion(combinacion.second);
            }
            hijos[i*2] = combinacion.first;
            hijos[i*2+1] = combinacion.second;
        }
        
        // 4. REEMPLAZO
        set<int>::iterator m = mejoresPadres.begin();
        int h = 0;     
        //recorro la poblacion
        for(int i = 0; i < POBLACION_SIZE; i++){
            if(m != mejoresPadres.end()){
                //si i no corresponde a uno de los mejores padres, inserto en su lugar
                //uno de los hijos generados
                if(i != *m){
                    poblacion[i] = hijos[h];
                    h++;
                    //me guardo su fitness
                    for(int j = 0; j < poblacion[i].size(); j++){
                    }
                    int c = evaluarSolucion(poblacion[i], flujos, distancias, aux);
                    fitnessPoblacion[i] = c;
                    evaluaciones++;
                }
                //si este es uno de los mejores padres, lo conservo
                else{
                    m++;
                }
            }
            //si ya he pasado todos los mejores padres, el resto de individuos los reemplazo
            else{
                poblacion[i] = hijos[h];
                h++;
                //me guardo su fitness
                int c = evaluarSolucion(poblacion[i], flujos, distancias, aux);
                fitnessPoblacion[i] = c;
                evaluaciones++;
            }
        }
        //ahora que he actualizado la población, actualizo también la lista
        //de los mejores padres, para la próxima iteración
        mejoresPadres.clear();
        fitnessMejoresPadres.clear();
        //Recorro cada individuo
        for(int i = 0; i < POBLACION_SIZE; i++){
            //tengo que separar a los mejores padres del conjunto
            //Los primeros, entran directamente al grupo
            if(i < POBLACION_SIZE - ESPERANZA_CRUCES*2){
                mejoresPadres.insert(i);
                fitnessMejoresPadres.insert(pair<int, int>(fitnessPoblacion[i], i) );
            }
            //para los siguientes, tengo que comprobar si son mejores que los que ya hay
            else{
                //otengo el padre con el peor fitness, que está al final del multiset
                multimap<int, int>::iterator it = fitnessMejoresPadres.end();
                it--;

                //si el individuo es mejor, lo inserto y borro el peor
                if(fitnessPoblacion[i] < (*it).first){
                    mejoresPadres.erase((*it).second);
                    mejoresPadres.insert(i);
                    fitnessMejoresPadres.erase(it);
                    fitnessMejoresPadres.insert(pair<int, int>(fitnessPoblacion[i], i));
                }                
            }                
        }
    
        T++;
    }
    
    //Fin de la evolución, escojo la mejor solución
    mejor_solucion = poblacion[fitnessMejoresPadres.begin()->second];
    
}
