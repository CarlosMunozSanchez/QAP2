#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void leerDatos(const string & entrada , vector<vector<int>> & flujos, vector<vector<int>> & distancias){
    
    flujos.clear();
    distancias.clear();
    
    fstream fs;
    char aux[6];
    
    fs.open(entrada);
    
    if(!fs.is_open()){
        cerr << "Error abriendo " << entrada << " ¿Es correcta la ruta?" << endl;
        exit(-1);
    }
    
    int tam;
    
    //leer el tamaño de la matriz
    fs.getline(aux, 64);
    tam = atoi(aux);
        
    //actualizar tamaños
    flujos.resize(tam);
    distancias.resize(tam);
    
    for(int i = 0; i < tam; i++){
        flujos[i].resize(tam);
        distancias[i].resize(tam);
    }
    
    //leer los datos para inicializar las matrices
    //matriz de flujos
    for(int i = 0; i < tam; i++){
        for(int j = 0; j < tam; j++){
            //me salto los espacios en blanco
            while(!isdigit(fs.peek())){
                fs.ignore();
            }
            
            string cadena;
            //escribo en la matriz
            while(isdigit(fs.peek())){
                cadena.push_back(fs.get());
            }
            flujos[i][j] = stoi(cadena);
        }
    }
    
    
    //matriz de distancias
    for(int i = 0; i < tam; i++){
        for(int j = 0; j < tam; j++){
            //me salto los espacios en blanco
            while(!isdigit(fs.peek())){
                fs.ignore();
            }
            
            string cadena;
            //escribo en la matriz
            while(isdigit(fs.peek())){
                cadena.push_back(fs.get());
            }
            distancias[i][j] = stoi(cadena);
        }
    }
            
    fs.close();
}

int leerSolucion(const string & entrada, vector<int> & solucion){
    int coste;
    int size;
    
    fstream fs;
    
    fs.open(entrada);
    
    if(!fs.is_open()){
        cerr << "Error abriendo " << entrada << " ¿Es correcta la ruta?" << endl;
        exit(-1);
    }
    
    //----------------------comienzo a leer-----------------------
    //me salto los espacios en blanco
    while(fs.peek() == ' ' or fs.peek() == '\n'){
        fs.ignore();
    }

    string cadena;
    //leo el tamaño
    while(fs.peek() != ' ' and fs.peek() != '\n'){
        cadena.push_back(fs.get());
    }
    
    size = stoi(cadena);    
    solucion.resize(size);
    
    for(int i = 0; i < size + 1; i++){
        //me salto los espacios en blanco
        while(!isdigit(fs.peek())){
            fs.ignore();
        }

        cadena.clear();
        //leo el tamaño
        while(isdigit(fs.peek())){
            cadena.push_back(fs.get());
        }
        
        if(i == 0){
            coste = stoi(cadena);
        }
        else{
            solucion[i-1] = stoi(cadena) - 1;
        }
    }
    
    return coste;    
}

void mostrarMatriz(const vector<vector<int>> & matriz){
    cout << endl;
    cout << "Tamaño de la matriz: " << matriz.size() << " x " << matriz[0].size() << endl;
    
    for(int i = 0; i < matriz.size(); i++){
        for(int j = 0; j < matriz[i].size(); j++){
            cout << matriz[i][j] << " ";
        }
        cout << endl;
    }
}

void mostrarVector(const vector<int> & v){
    cout << endl;
    cout << "Tamaño del vector: " << v.size() << endl;
    
    for(int j = 0; j < v.size(); j++){
        cout << v[j] << " ";
    }
    cout << endl;
    
}

int evaluarSolucion(const vector<int> & permutacion,    
        const vector<vector<int>> & flujos, const vector<vector<int>> & distancias, 
        float & fitness){
    int sol = 0;
    
    for(int i = 0; i < permutacion.size(); i++){
        for(int j = 0; j < permutacion.size(); j++){
        
            sol += flujos[i][j] * distancias[permutacion[i]][permutacion[j]];        
        }        
    }
    
    //calcular el fitness
    fitness = (sol - fitness) / fitness;
    
    return sol;
}