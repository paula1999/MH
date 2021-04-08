#include <iostream>
#include "PAR.h"
#include "random.h"
#include <string>
#include <chrono>
#include <fstream>

using namespace std;

void calcular_PAR (const string datos, const string restricciones, const int num_clusters, const string algoritmo, int semilla){
    vector<int> solucion;
    PAR par(datos, restricciones, num_clusters);
    string rutaSol;
    fstream fSol;

    Set_random(semilla);
    rutaSol = restricciones + "_" + algoritmo + "_" + to_string(Get_random()) + ".out";
    fSol.open(rutaSol, fstream::out);

    chrono::system_clock::time_point start = chrono::system_clock::now();
    
    if (algoritmo.compare("Greedy") == 0)
        solucion = par.greedy();
    else if (algoritmo.compare("BL") == 0)
        solucion = par.busquedaLocal();
    
    chrono::system_clock::time_point stop = chrono::system_clock::now();
    chrono::milliseconds duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    
    if (fSol.is_open()){
        fSol << "Tiempo: " << duration.count() << " milisegundos\n";
        fSol << par;

        cout << "Resultados copiados en " << rutaSol << "\n";

        fSol.close();
    }
    else   
        cerr << "Error al abrir el fichero " << rutaSol << "\n";   
}

int main (int argc, char ** argv){
    string datos, restricciones, ruta_data = "data/";
    int num_clusters;
    unsigned long semilla;

    if (argc < 5){
        cerr << "ERROR en los argumentos.\n";
        cerr << "Formato: " << argv[0] << " <fichero_datos.dat> <fichero_restricciones.const> <numero_clusters> <semilla>\n";

        exit(-1);
    }

    datos = ruta_data + argv[1];
    restricciones = ruta_data + argv[2];
    num_clusters = atoi(argv[3]);
    semilla = atoi(argv[4]);

    cout << "\nEjecutando algoritmo Greedy...\n";

    calcular_PAR (datos, restricciones, num_clusters, "Greedy", semilla);

    cout << "\nEjecutando algoritmo Busqueda Local...\n";

    calcular_PAR (datos, restricciones, num_clusters, "BL", semilla);
}