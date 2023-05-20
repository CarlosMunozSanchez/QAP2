/* 
 * File:   AGE.h
 * Author: carlos
 *
 * Created on 17 de abril de 2023, 19:23
 */

#ifndef AGE_H
#define AGE_H

#include <vector>
#include <set>
#include <map>

class AGE {
private:
    const int POBLACION_SIZE = 50;
    const int MAX_EVAL = 50000;
    //const float PCRUCE = 0.7;
    const float PMUTACION = 0.1;
    //Esta es la cantidad de cruces que voy a realizar
    const int ESPERANZA_CRUCES = 1;
    //Esta es la cantidad de padres cuyos hijos sufrirán una mutación
    //voy a mutar 1 pareja de hijos de cada 10 que genero
    const int ESPERANZA_MUTACION = ESPERANZA_CRUCES * PMUTACION;
    
    int tipo; //0 -> posicion 1 -> PMX
    std::vector<int> mejor_solucion; //indice de la mejor solución
    
    std::vector<std::vector<int>> poblacion;
    std::vector<int> fitnessPoblacion;
    
    //Me guardo los fitness pero ordenados para poder acceder siempre a los peores
    std::multimap<int, int> fitnessOrdenados;
    
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
     * el esquema estacionario.
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
    AGE(int tipo, int n_genes, const std::vector<std::vector<int>> & flujos, 
            const std::vector<std::vector<int>> & distancias, int seed = 42);
    
    /**
     * @brief Obtiene la mejor solución
     * @return vector<int> con la mejor solución encontrada tras la simulación.
     */
    inline std::vector<int> getSolucion(){
        return mejor_solucion;
    }

    

};

#endif /* AGE_H */

