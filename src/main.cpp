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
    //rutaSol = restricciones + "_" + algoritmo + "_" + to_string(Get_random()) + ".out";
    rutaSol = restricciones + "_" + algoritmo + ".out";
    //fSol.open(rutaSol, fstream::out);
    fSol.open(rutaSol, fstream::app);

    chrono::system_clock::time_point start = chrono::system_clock::now();
    
    if (algoritmo.compare("Greedy") == 0)
        solucion = par.greedy();
    else if (algoritmo.compare("BL") == 0)
        solucion = par.busquedaLocal(par.solAleatoria(), 100000);
    else if (algoritmo.compare("AGG-UN") == 0)
        solucion = par.algoritmoGenetico(50, "G", "UN", 0.7);
    else if (algoritmo.compare("AGG-SF") == 0)
        solucion = par.algoritmoGenetico(50, "G", "SF", 0.7);
    else if (algoritmo.compare("AGE-UN") == 0)
        solucion = par.algoritmoGenetico(50, "E", "UN", 1.0);
    else if (algoritmo.compare("AGE-SF") == 0)
        solucion = par.algoritmoGenetico(50, "E", "SF", 1.0);
    else if (algoritmo.compare("AM-(10,1.0)") == 0)
        solucion = par.algoritmoMemetico(50, "1.0", 10, 0.7);
    else if (algoritmo.compare("AM-(10,0.1)") == 0)
        solucion = par.algoritmoMemetico(50, "0.1", 10, 0.7);
    else if (algoritmo.compare("AM-(10,0.1mej)") == 0)
        solucion = par.algoritmoMemetico(50, "0.1mej", 10, 0.7);
    else if (algoritmo.compare("ES") == 0)
        solucion = par.ES(par.solAleatoria(), 100000);
    else if (algoritmo.compare("BMB") == 0)
        solucion = par.BMB(10, 10000);
    else if (algoritmo.compare("ILS") == 0)
        solucion = par.ILS(10000, 10, "ILS", 10);
    else if (algoritmo.compare("ILS-ES") == 0)
        solucion = par.ILS(10000, 10, "ILS-ES", 10);
    
    chrono::system_clock::time_point stop = chrono::system_clock::now();
    chrono::milliseconds duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    
    if (fSol.is_open()){
        //fSol << "Tiempo: " << duration.count() << " milisegundos\n";
        //fSol << par << endl;
        fSol << semilla << " : " << par << duration.count() << endl;
        //fSol << duration.count() << endl;

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

    ///////////////////////////////// PR??CTICA 1 /////////////////////////////////
    /*
    cout << "\nEjecutando algoritmo Greedy...\n";

    calcular_PAR (datos, restricciones, num_clusters, "Greedy", semilla);

    cout << "\nEjecutando algoritmo Busqueda Local...\n";

    calcular_PAR (datos, restricciones, num_clusters, "BL", semilla);
    */

    ///////////////////////////////// PR??CTICA 2 /////////////////////////////////
    /*
    cout << "\nEjecutando algoritmo AGG-UN ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "AGG-UN", semilla);

    cout << "\nEjecutando algoritmo AGG-SF ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "AGG-SF", semilla);

    cout << "\nEjecutando algoritmo AGE-UN ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "AGE-UN", semilla);

    cout << "\nEjecutando algoritmo AGE-SF ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "AGE-SF", semilla);

    cout << "\nEjecutando algoritmo AM-(10,1.0) ...\n";

    calcular_PAR (datos, restricciones, num_clusters, "AM-(10,1.0)", semilla);

    cout << "\nEjecutando algoritmo AM-(10,0.1) ...\n";

    calcular_PAR (datos, restricciones, num_clusters, "AM-(10,0.1)", semilla);

    cout << "\nEjecutando algoritmo AM-(10,0.1mej) ...\n";

    calcular_PAR (datos, restricciones, num_clusters, "AM-(10,0.1mej)", semilla);
    */

    ///////////////////////////////// PR??CTICA 3 /////////////////////////////////
    
    cout << "\nEjecutando algoritmo ES ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "ES", semilla);

    cout << "\nEjecutando algoritmo BMB ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "BMB", semilla);

    cout << "\nEjecutando algoritmo ILS ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "ILS", semilla);

    cout << "\nEjecutando algoritmo ILS-ES ...\n";
    
    calcular_PAR (datos, restricciones, num_clusters, "ILS-ES", semilla);
}