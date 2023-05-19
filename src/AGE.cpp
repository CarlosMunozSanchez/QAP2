/* 
 * File:   AGE.cpp
 * Author: carlos
 * 
 * Created on 17 de abril de 2023, 19:23
 */

#include "AGE.h"
#include "random.hpp"
#include "funciones.h"
#include <unordered_map>
#include <iostream>
// get base random alias which is auto seeded and has static API and internal state
// it is not threads secure, you can also use ::random_thread_local
using Random = effolkronium::random_static;
using namespace std;

AGE::AGE(int tipo, int n_genes, const vector<vector<int>> & flujos, 
        const vector<vector<int>> & distancias, int seed) {
    
    this->tipo = tipo;
    coste_mejor = 999999999;
    
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
    
pair<vector<int>, vector<int>> AGE::cruce(const vector<int> & c1, 
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

void AGE::mutacion(vector<int> & cromosoma){
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

void AGE::simularEvolucion(const vector<vector<int> >& flujos, 
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
                
                fitnessOrdenados.insert(pair<int, int>(c, i));
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
            
            //muto 1 de cada diez cruces
            if(T % (int)(PMUTACION*100)){
                //para los primeros ESPERANZA_MUTACION padres, muto a sus hijos
                mutacion(combinacion.first);
                mutacion(combinacion.second);
            }
            hijos[i*2] = combinacion.first;
            hijos[i*2+1] = combinacion.second;
        }
        
        // 4. REEMPLAZO
        //Los hijos compiten con la población para entrar
        
        //Evaluo a los hijos;
        int h0 = evaluarSolucion(hijos[0], flujos, distancias, aux);
        int h1 = evaluarSolucion(hijos[1], flujos, distancias, aux);
        evaluaciones += 2;
        
        //también los ordeno
        if(h1 < h0){
            swap(h0, h1);
            swap(hijos[0], hijos[1]);
        }
        
        //obtengo los peores padres
        multimap<int, int>::iterator peor, segundoPeor = fitnessOrdenados.end();
        segundoPeor--;
        peor = segundoPeor--;
        
        //Compiten hijos vs peores padres y entran los mejores
        
        //mejor hijo v mejor padre
        if(h0 < segundoPeor->first){
            //sale el peor padre, entra el mejor hijo
            poblacion[peor->second] = hijos[0];
            fitnessPoblacion[peor->second] = h0;
            fitnessOrdenados.insert(pair<int, int>(h0, peor->second));
            //2ºhijo v peor padre
            if(h1 < peor->first){
                //entra también el segundo hijo, sale el otro padre
                poblacion[segundoPeor->second] = hijos[1];
                fitnessPoblacion[segundoPeor->second] = h1;
                fitnessOrdenados.insert(pair<int, int>(h1, segundoPeor->second));
                fitnessOrdenados.erase(segundoPeor);
            }
            //else -> entra sólo el mejor hijo y sale el peor padre
            fitnessOrdenados.erase(peor);
        }
        //El mejor hijo pierde contra mejor padre
        //mejor hijo v peor padre
        else if(h0 < peor->first){
            //sale el peor padre, entra el mejor hijo
            poblacion[peor->second] = hijos[0];
            fitnessPoblacion[peor->second] = h0;
            fitnessOrdenados.insert(pair<int, int>(h0, peor->second));
            fitnessOrdenados.erase(peor);
        }
        //else -> no entra ningún hijo
        
        T++;
    }
    
    //Fin de la evolución, escojo la mejor solución
    mejor_solucion = poblacion[fitnessOrdenados.begin()->second];
    
}
