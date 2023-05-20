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
    //Esta es la cantidad de cruces que voy a realizar
    const int ESPERANZA_CRUCES = PCRUCE * POBLACION_SIZE / 2;
    //Esta es la cantidad de padres cuyos hijos sufrirán una mutación
    const int ESPERANZA_MUTACION = ESPERANZA_CRUCES * PMUTACION;
    
    int tipo; //0 -> posicion 1 -> PMX
    std::vector<int> mejor_solucion; //mejor solución
    
    std::vector<std::vector<int>> poblacion;
    std::vector<int> fitnessPoblacion;
    
    //en cada generación, habrá ESPERANZA_CRUCES cruces tal que la cantidad
    //de hijos obtenidos será menor al tamaño de la población
    //La nueva población incluye estos hijos. El resto de individuos hasta
    //llegar a POBLACION_SIZE serán los mejores individuos de la generación anterior

    //creo estructuras para tener a los mejores padres ordenados por su fitness
    //y con su indice asociado en el array poblacion
    
    //En este set, tengo los índices ordenados por posición
    std::set<int> mejoresPadres;
    //en este multimap, tengo los fitness ordenados de mejor a peor, con su indice asociado
    std::multimap<int, int> fitnessMejoresPadres;
    
    /**
     * @brief Operador de cruce. Dados dos cromosomas, calcula la descendencia 
     * que estos producen.
     * @param c1 Primer cromosoma.
     * @param c2 Segundo cromosoma.
     * @param tipo Indica el tipo de cruce. 0 -> posicion, 1 -> PMX.
     * @return Los dos hijos generados por c1 y c2.
     */
    std::pair<std::vector<int>, std::vector<int>> cruce(const std::vector<int> & c1, 
        const std::vector<int> & c2, int tipo);
    
    /**
     * @brief Operador de mutación. Dado un cromosoma, realiza una permutación 
     * aleatoria de dos genes. 
     * @param cromosoma Cromosoma que sufre la mutación.
     */
    void mutacion(std::vector<int> & cromosoma);
    
    /**
     * @brief Método que lleva el control de flujo del algoritmo. Mantiene y
     * actualiza la población utilizando los operadores anteriores. Implementa
     * el esquema generacional.
     * @param flujos Matriz de flujos asociada al problema.
     * @param distancias Matriz de distancias asociada al problema.
     * @param tipo Indicador del tipo de cruce. 0 -> posición, 1 -> PMX.
     */
    void simularEvolucion(const std::vector<std::vector<int>> & flujos, 
                    const std::vector<std::vector<int>> & distancias, int tipo);
    
public:
    /**
     * @brief Constructor de la clase.
     * @param tipo Tipo de operador de cruce. 0 -> posición, 1 -> cruce.
     * @param n_genes Nº de genes en cada cromosoma
     * @param flujos Matriz de flujos asociada.
     * @param distancias Matriz de distancias asociada.
     * @param seed Semilla para el generador de números aleatorios
     */
    AGG(int tipo, int n_genes, const std::vector<std::vector<int>> & flujos, 
            const std::vector<std::vector<int>> & distancias, int seed = 42);
    
    /**
     * @brief Obtiene la mejor solución
     * @return vector<int> con la mejor solución encontrada tras la simulación.
     */
    inline std::vector<int> getSolucion(){
        return mejor_solucion;
    }

    

};

#endif /* AGG_H */

