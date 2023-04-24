/* 
 * File:   AGG.h
 * Author: carlos
 *
 * Created on 17 de abril de 2023, 19:23
 */

#ifndef AGG_H
#define AGG_H

#include <vector>
#include "funciones.h"

class AGG {
private:
    const int POBLACION_SIZE = 50;
    const int MAX_EVAL = 50000;
    const float PCRUCE = 0.7;
    const float PMUTACION = 0.1;
    
    int tipo; //0 -> posicion 1 -> PMX
    int mejor_solucion = 0;
    int coste_mejor = 99999999;
    
    std::vector<std::vector<int>> poblacion;
    
    std::vector<int> seleccion(std::vector<int> & c1, std::vector<int> & c2);
    
    std::pair<std::vector<int>, std::vector<int>> cruce(std::vector<int> & c1, 
        std::vector<int> & c2, int tipo);
    
    void mutacion(std::vector<int> & cromosoma);
    
public:
    AGG(int tipo, int seed);
    
    

    

};

#endif /* AGG_H */

