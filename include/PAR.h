#ifndef PAR_H
#define PAR_H

#include <iostream>
#include <vector>
#include <list>
#include <tuple>

using namespace std;

class PAR{
    public:
        PAR (const string f_datos, const string f_restricciones, const int k);
        void leerFicheros (const string fDatos, const string fRestricciones);
        void leerFicheros2 (const string fDatos, const string fRestricciones);
        double fitness (vector<int> solucion);
        double fitnessBL (vector<int> solucion);
        double distanciaIntracluster (int cluster, vector<int> solucion);
        double desviacionGeneral (vector<int> solucion);
        double distancia (int posicionPunto, int cluster);
        int infeasibilityGreedy (vector<int> solucion);
        int infeasibilityGreedy (int posicionPunto, int cluster);
        int infeasibilityBL (vector<int> solucion);
        void actualizarCentroides ();
        bool cambioCluster (vector<int> & solucion, int indice, int cluster);
        bool clusterVacio (vector<int> solucion);
        vector<int> greedy ();
        vector<int> busquedaLocal ();
        void imprimirCentroides ();
        void imprimirClusters ();
        void imprimirDistancias ();
        void calcularDistancias ();
        void imprimirRestricciones ();
        bool parValido (pair<int, int> par, vector<int> solucion);

    private:
        vector< vector<int> > restricciones; // Matriz de restricciones
        vector< vector<int> > listaRestricciones; // Lista de restricciones
        
        vector< vector<double> > datos; // Matriz de datos
        vector<int> clusters; // Contiene los indices a los clusters de cada punto
        vector< vector<double> > centroides; // Matriz con los centros de cada cluster

        
        vector< vector<int> > restriccionesML;
        vector< vector<int> > restriccionesCL;

        vector< vector<double> > distancias;

        
        void imprimirDatos ();
        

        int num_clusters;
};

#endif