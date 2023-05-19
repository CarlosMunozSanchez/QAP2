/* 
 * File:   main.cpp
 * Author: carlos
 *
 * Created on 17 de abril de 2023, 18:46
 */
#include <iostream>
#include <chrono>
#include "AGG.h"
#include "AGE.h"
#include "AM.h"
#include "funciones.h"

using namespace std;
using namespace std::chrono;

/*
 * 
 */
int main(int argc, char** argv){
    
    if(argc < 3 or argc > 4){
        cerr << "Uso: " << argv[0] << " [archivo de datos] [semilla] [archivo de solucion]" << endl;
        exit(-1);
    }
    
    //lectura de datos de entrada
    string entrada(argv[1]);
    int seed = atoi(argv[2]);
    string sol;
    vector<int> optimo;
    int costeOptimo = 1;
    
    if(argc == 4){
        sol = argv[3];
        costeOptimo = leerSolucion(sol, optimo);
    }
    
    vector<vector<int>> flujos, distancias;
    
    leerDatos(entrada, flujos, distancias);
    
    vector<int> solucion;
    
    cout << "---------------------------------------------------------------" << endl;
    cout << "Resultados del archivo " << entrada << endl << endl;
    
    //AGG con cruce por posicion
    auto momentoInicio = high_resolution_clock::now();
    AM generacionalPosicion(1, 0, optimo.size(), flujos, distancias, seed);
    
    solucion = generacionalPosicion.getSolucion();
    auto momentoFin = high_resolution_clock::now();
    
    float fitness = costeOptimo;
    cout << "Solución AGG con cruce por posicion con coste " << evaluarSolucion(solucion, flujos, distancias, fitness) <<
            " y fitness = " << fitness << endl;    
    mostrarVector(solucion);
    
    microseconds tiempo = duration_cast<std::chrono::microseconds>(momentoFin - momentoInicio);
    cout << "Tiempo Pasado (ms): " << tiempo.count() / 1000.0 << endl;
    
    
    cout << "---------------------------------------------------------------" << endl;
    
    //AGG con cruce PMX
     momentoInicio = high_resolution_clock::now();
    AM generacionalPMX(1, 1, optimo.size(), flujos, distancias, seed);
    
    solucion = generacionalPMX.getSolucion();
    momentoFin = high_resolution_clock::now();
    
    fitness = costeOptimo;
    cout << "Solución AGG con cruce por posicion con coste " << evaluarSolucion(solucion, flujos, distancias, fitness) <<
            " y fitness = " << fitness << endl;    
    mostrarVector(solucion);
    
    tiempo = duration_cast<std::chrono::microseconds>(momentoFin - momentoInicio);
    cout << "Tiempo Pasado (ms): " << tiempo.count() / 1000.0 << endl;
    
    
    cout << "---------------------------------------------------------------" << endl;
    
}