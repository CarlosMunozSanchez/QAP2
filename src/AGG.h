/* 
 * File:   AGG.h
 * Author: carlos
 *
 * Created on 17 de abril de 2023, 19:23
 */

#ifndef AGG_H
#define AGG_H

#include <vector>
#include <set>
#include <map>

class AGG {
private:
    const int POBLACION_SIZE = 50;
    const int MAX_EVAL = 50000;
    const float PCRUCE = 0.7;
    const float PMUTACION = 0.1;
    //Esta es la cantidad de padres que tengo que obtener y cruzar
    const int ESPERANZA_CRUCES = PCRUCE * POBLACION_SIZE / 2;
    //Esta es la cantidad de padres cuyos hijos sufrirán una mutación
    const int ESPERANZA_MUTACION = ESPERANZA_CRUCES * POBLACION_SIZE;
    
    int tipo; //0 -> posicion 1 -> PMX
    std::vector<int> mejor_solucion; //indice de la mejor solución
    int coste_mejor;
    
    std::vector<std::vector<int>> poblacion;
    std::vector<int> fitnessPoblacion;
    
    //en cada generación, habrá ESPERANZA_CRUCES cruces tal que la cantidad
    //de hijos obtenidos será menor al tamaño de la población
    //La nueva población incluye estos hijos. El resto de individuos hasta
    //llegar a POBLACION_SIZE serán los mejores individuos de la generación anterior

    //creo sets para tener a los mejores padres ordenados por su fitness
    //y con su indice asociado en el array poblacion
    
    //En este set, tengo los índices ordenados por posición
    std::set<int> mejoresPadres;
    //en este multimap, tengo los fitness ordenados de mejor a peor, con su indice asociado
    std::multimap<int, int> fitnessMejoresPadres;
        
    std::pair<std::vector<int>, std::vector<int>> cruce(const std::vector<int> & c1, 
        const std::vector<int> & c2, int tipo);
    
    void mutacion(std::vector<int> & cromosoma);
    
    void simularEvolucion(const std::vector<std::vector<int>> & flujos, 
                    const std::vector<std::vector<int>> & distancias, int tipo);
    
public:
    AGG(int tipo, int n_genes, const std::vector<std::vector<int>> & flujos, 
            const std::vector<std::vector<int>> & distancias, int seed = 42);
    
    inline std::vector<int> getSolucion(){
        return mejor_solucion;
    }

    

};

#endif /* AGG_H */

