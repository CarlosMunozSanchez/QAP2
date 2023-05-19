/* 
 * File:   QAPBL.h
 * Author: carlos
 *
 * Created on 5 de abril de 2023, 13:18
 */

#include <vector>

#ifndef QAPBL_H
#define QAPBL_H

class QAPBL {

    private:
        std::vector<int> solucion;
        std::vector<bool> dlb;
        
        int coste;
        
        /**
         * @brief Calcula el cambio en el coste al intercambiar los elementos i
         * y j del vector solución, pero sin realizar ese cambio.
         * @param i Elemento i.
         * @param j Elemento j.
         * @param flujos Matriz de flujos asociada.
         * @param distancias Matriz de distancias asociada.
         * @return El cambio en el coste que se produce si se permutan i y j.
         */
        int comprobarMovimiento(int i, int j, const std::vector<std::vector<int>> & flujos, const std::vector<std::vector<int>> & distancias);
        
        /**
         * @brief Intercambia los elementos i y j del vector solución.
         * @param i Elemento i.
         * @param j Elemento j.
         */
        void aplicarMovimiento(int i, int j);
        
        /**
         * @brief Algoritmo de búsqueda local para el problema QAP.
         * @param flujos Matriz de flujos asociada.
         * @param distancias Matriz de distancias asociada.
         */
        void busquedaLocal(const std::vector<std::vector<int>> & flujos, const std::vector<std::vector<int>> & distancias);
        
    public:
        
        /**
         * @brief Constructor sin solución inicial (genera una aleatoriamente)
         * @param flujos
         * @param distancias
         * @param seed
         */
        QAPBL(const std::vector<std::vector<int>> & flujos, const std::vector<std::vector<int>> & distancias, int seed);
        /**
         * @brief Constructor con solución inicial
         * @param flujos
         * @param distancias
         * @param seed
         */
        QAPBL(const std::vector<int> & inicial, const std::vector<std::vector<int>> & flujos, 
            const std::vector<std::vector<int>> & distancias, int seed);
      
        std::vector<int> getSolucion()const {
            return solucion;
        }
        
        int getCoste() const{
            return coste;
        }
    


};

#endif /* QAPBL_H */

