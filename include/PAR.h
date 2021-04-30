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
        void leerFicheroDatos (const string fDatos);
        void leerFicheroRestricciones (const string fRestricciones);
        double fitnessGreedy ();
        double fitnessBL ();
        double distanciaIntracluster (int cluster, const vector<int> solucion);
        double desviacionGeneral ();
        double distancia (int posicionPunto, int cluster) const;
        int infeasibilityGreedy () const;
        int infeasibilityCluster (int posicionPunto, int cluster) const;
        int infeasibilityBL () const;
        int infeasibilityBL (int posicionPunto, int cluster) const;
        void actualizarCentroides ();
        bool cambioCluster (vector<int> & solucion, int indice, int cluster);
        bool clusterVacio (vector<int> solucion) const;
        vector<int> greedy ();
        vector<int> busquedaLocal ();
        void imprimirCentroides () const;
        void imprimirClusters () const;
        void imprimirDistancias () const;
        void calcularDistancias ();
        void imprimirRestricciones () const;
        void imprimirDatos () const;
        bool parValido (pair<int, int> par, vector<int> solucion) const;
        void calcularDistMax();
        vector<int> busquedaLocalSuave (vector<int> solucion, const int nFallosMAX);
        vector<int> algoritmoGenetico (int M, const string evolucion, const string operadorCruce, const double probCruce, const double probMutacion);
        vector<int> algoritmoMemetico (const int M, const double probCruce, const double probMutacion);
        int calcularMejorCromosoma (const vector<double> pFitness);
        int calcularPeorCromosoma (const vector<double> pFitness);
        int operadorSeleccion (const int c1, const int c2, const vector<double> pFitness);
        vector<int> operadorCruceUN (const vector<int> padre1, const vector<int> padre2);
        vector<int> operadorCruceSF (const vector<int> padre1, const vector<int> padre2);
        vector<int> operadorMutacionUN (vector<int> cromosoma);
        void repararCromosoma (vector<int> & cromosoma);

        friend ostream & operator << (ostream & flujo, const PAR & par);

    private:
        vector< vector<int> > restricciones; // Matriz de restricciones
        vector< vector<int> > restriccionesML; // Lista de restricciones Must Link
        vector< vector<int> > restriccionesCL; // Lista de restricciones Cannot Link
        vector< vector<double> > datos; // Matriz de datos
        vector<int> clusters; // Contiene los indices a los clusters de cada punto
        vector< vector<double> > centroides; // Matriz con los centros de cada cluster        
        vector< vector<double> > distancias; // Matriz de distancias
            
        int num_clusters, infeasibility, restrMax;
        double distMax, lambda, desvGeneral, fitness;
        double distOptima;
};

#endif