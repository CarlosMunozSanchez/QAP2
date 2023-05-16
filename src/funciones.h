/* 
 * File:   funciones.h
 * Author: carlos
 *
 * Created on 5 de abril de 2023, 12:21
 */
#ifndef FUNCIONES_H
#define FUNCIONES_H

/**
 * @brief Función que lee los datos del fichero de entrada y actualiza.
 * los vectores flujo y distancia con dichos valores.
 * @param entrada Ruta del fichero de entrada.
 * @param flujos Matriz de flujos.
 * @param distancias Matriz de distancias.
 */
void leerDatos(const std::string & entrada , std::vector<std::vector<int>> & flujos, std::vector<std::vector<int>> & distancias);

/**
 * @brief Función para leer un archivo de solución óptima.
 * @param entrada Ruta del fichero de entrada.
 * @param solucion Vector en el que se guardará la solución.
 * @return El coste de la solución óptima
 */
int leerSolucion(const std::string & entrada, std::vector<int> & solucion);

/**
 * @brief Función para mostrar por pantalla una matriz.
 * @param matriz  Matriz que se desea mostrar.
 */
void mostrarMatriz(const std::vector<std::vector<int>> & matriz);

/**
 * @brief Función para mostrar por pantalla un vector.
 * @param vector Vector que se desea mostrar.
 */
void mostrarVector(const std::vector<int> & vector);

/**
 * @brief Función que calcula el coste de una solución. 
 * @param permutacion Solución propuesta al problema.
 * @param flujos Matriz de flujos asociada.
 * @param distancias Matriz de distancias asociada.
 * @param fitness Indica el costo de la solución óptima. Al acabar la función, 
 * se sobrescribe con el fitness
 * @return El coste asociado al vector solución proporcionado
 */
int evaluarSolucion(const std::vector<int> & permutacion, 
        const std::vector<std::vector<int>> & flujos, const std::vector<std::vector<int>> & distancias, 
        float & fitness);

#endif /* FUNCIONES_H */

