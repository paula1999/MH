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
        double fitness (vector<int> solucion);
        double distanciaIntracluster (int cluster);
        double desviacionGeneral (vector<int> solucion);
        double distancia (int posicionPunto, int cluster);
        int infeasibilityGreedy (vector<int> solucion);
        int infeasibilityGreedy (int posicionPunto, int cluster);
        int infeasibilityBL (vector<int> solucion);
        void actualizarCentroides ();
        void cambioCluster (vector<int> solucion, int indice, int cluser);
        bool clusterVacio (vector<int> solucion);
        vector<int> greedy ();
        vector<int> busquedaLocal ();
        void imprimirCentroides ();
        void imprimirClusters ();
        void imprimirDistancias ();
        void calcularDistancias ();

    private:
        vector< vector<int> > restricciones; // Matriz de restricciones
        
        vector< vector<double> > datos; // Matriz de datos
        vector<int> clusters; // Contiene los indices a los clusters de cada punto
        vector< vector<double> > centroides; // Matriz con los centros de cada cluster

        list <tuple <int, int, int> > listaRestricciones;
        vector<int> restriccionesML;
        vector<int> restriccionesCL;

        vector< vector<double> > distancias;

        void imprimirRestricciones ();
        void imprimirDatos ();
        

        int num_clusters;
};

#endif