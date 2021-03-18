#include <iostream>
#include "PAR.h"
#include "random.h"
#include <string>
#include <chrono>

using namespace std;

void calcular_PAR (const string datos, const string restricciones, const int num_clusters, const string algoritmo){
    vector<int> solucion;
    PAR par(datos, restricciones, num_clusters);
    double tiempo = 0.0;

    // GREEDY
    solucion = par.greedy();
    //cout << "\n\nTIEMPO: " << tiempo;
    //par.imprimirRestricciones();
    
    
    cout << "\n\nLista de clusters:\n";

    for (int i = 0; i < solucion.size(); i++)
        cout << " " << solucion[i];
    
    cout << "\nRestricciones incumplidas: " << par.infeasibilityGreedy(solucion) << endl;

    cout << "\nCENTROIDES:\n";
    par.imprimirCentroides();
    cout << "\n\nCLUSTERS:\n";
    par.imprimirClusters ();

    cout << "\n\nDISTANCIA INTRA-CLUSTER:\n";

    for (int i = 0; i < 7; i++)
        cout << "\nDistancia intra-cluster " << i << ": " << par.distanciaIntracluster(i);
    par.calcularDistancias();
    par.fitness(solucion);
    //par.imprimirDistancias();

    //par.imprimirCentroides();
    //cout << par.desviacionGeneral(solucion);
    
   // BUSQUEDA LOCAL
    //solucion = par.busquedaLocal();
    //cout << "\nBL:\nLista de clusters:\n";

    //for (int i = 0; i < solucion.size(); i++)
        //cout << " " << solucion[i];

    //cout << "\nfitnessBL: " << par.fitness(solucion) << endl;
    cout << endl;
    
}

int main (int argc, char ** argv){
    string datos, restricciones;
    int num_clusters;

    if (argc < 4){
        cerr << "ERROR en los argumentos.\n";
        cerr << "Formato: " << argv[0] << " <fichero_datos.dat> <fichero_restricciones.const> <numero_clusters> \n";

        exit(-1);
    }

    datos = argv[1];
    restricciones = argv[2];
    num_clusters = atoi(argv[3]);

    Set_random(unsigned(3242355));

    calcular_PAR (datos, restricciones, num_clusters, "GREEDY");

}