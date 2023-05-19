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
    int coste_mejor;
    
    std::vector<std::vector<int>> poblacion;
    std::vector<int> fitnessPoblacion;
    
    //Me guardo los fitness pero ordenados para poder acceder siempre a los peores
    std::multimap<int, int> fitnessOrdenados;
    
    std::pair<std::vector<int>, std::vector<int>> cruce(const std::vector<int> & c1, 
        const std::vector<int> & c2, int tipo);
    
    void mutacion(std::vector<int> & cromosoma);
    
    void simularEvolucion(const std::vector<std::vector<int>> & flujos, 
                    const std::vector<std::vector<int>> & distancias, int tipo);
    
public:
    AGE(int tipo, int n_genes, const std::vector<std::vector<int>> & flujos, 
            const std::vector<std::vector<int>> & distancias, int seed = 42);
    
    inline std::vector<int> getSolucion(){
        return mejor_solucion;
    }

    

};

#endif /* AGE_H */

