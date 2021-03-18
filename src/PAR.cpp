#include <iostream>
#include "PAR.h"
#include <fstream>
#include <sstream>
#include "random.h"
#include <algorithm>
#include <math.h>
#include <tuple>

using namespace std;

void PAR::imprimirRestricciones (){
    for (int i = 0; i < restricciones.size(); i++){
        
        for (int j = 0; j < restricciones[i].size(); j++)
            cout << restricciones[i][j] << " ";
        
        cout << "\n";
    }
}

void PAR::imprimirDatos (){
    for (int i = 0; i < datos.size(); i++){
        for (int j = 0; j < datos[i].size(); j++)
            cout << datos[i][j] << " ";
        
        cout << "\n";
    }
}

void PAR::imprimirCentroides (){
    for (int i = 0; i < centroides.size(); i++){
        cout << "\nCentroide " << i << ": ";
        for (int j = 0; j < centroides[i].size(); j++)
            cout << centroides[i][j] << " ";
    }
}

void PAR::imprimirClusters (){
    for (int i = 0; i < num_clusters; i++){
        cout << "\nCluster " << i << ": ";
        for (int j = 0; j < clusters.size(); j++){
            if (clusters[j] == i)
                cout << j << " ";
        }
    }
}

void PAR::imprimirDistancias (){
    for (int i = 0; i < distancias.size(); i++){
        cout << "\n";
        for (int j = 0; j < distancias[i].size(); j++){
            cout << distancias[i][j] << " ";
        }
    }
}

PAR::PAR (const string fDatos, const string fRestricciones, const int k){
    num_clusters = k;

    leerFicheros (fDatos, fRestricciones);
    clusters.resize(datos.size(), -1);

    // Inicializo lista de restricciones
    for (int i = 0; i < restricciones.size(); i++)
        for (int j = i+1; j < restricciones[i].size(); j++)
            listaRestricciones.push_back(make_tuple(i, j, restricciones[i][j]));

    // Inicializo matriz de distancias
    vector<double> aux(datos.size(), 0.0);
    for (int i = 0; i < datos.size(); i++)
            distancias.push_back(aux);


}

void PAR::calcularDistancias (){
    double dist = 0.0;

    for (int i = 0; i < datos.size(); i++){
        for (int j = i+1; j < datos.size(); j++){
            dist = 0.0;
            for (int k = 0; k < datos[i].size(); k++)
                dist += (datos[i][k] - datos[j][k]) * (datos[i][k] - datos[j][k]);

            distancias[i][j] = sqrt(dist);
        }
    }
}

void PAR::leerFicheros (const string fDatos, const string fRestricciones){
    string linea, valor;
    ifstream ifDatos(fDatos), ifRestricciones(fRestricciones);
    int contador = 0;

    // Cargo el fichero de datos
    if (ifDatos.is_open()){
        datos.resize(contador);

        while (getline(ifDatos, linea)){
            stringstream cadena(linea);
            datos.resize(contador + 1);
        
            while (getline(cadena, valor, ','))
                datos[contador].push_back(stod(valor));
        
            contador++;
        }

        ifDatos.close();
    }
    else   
        cerr << "Error al abrir el fichero " << fDatos << "\n";

    contador = 0;

    // Cargo el fichero de restricciones
    if (ifRestricciones.is_open()){
        restricciones.resize(contador);

        while (getline(ifRestricciones, linea)){
            stringstream cadena(linea);
            restricciones.resize(contador + 1);
        
            while (getline(cadena, valor, ','))
                restricciones[contador].push_back(stod(valor));
        
            contador++;
        }

        ifRestricciones.close();
    }
    else   
        cerr << "Error al abrir el fichero " << fRestricciones << "\n";

    
}

double PAR::fitness (vector<int> solucion){
    double lambda, dist = 0.0, distmax = 0.0;
    int restrmax = 0;

    // DUDA hacer funcion para calcular lambda y si eso meterlo en una variable
    // Calculo mayor distancia
    for (int i = 0; i < datos.size(); i++){
        for (int j = i+1; j < datos.size(); j++){
            dist = 0.0;
            for (int k = 0; k < datos[i].size(); k++)
                dist += (datos[i][k] - datos[j][k]) * (datos[i][k] - datos[j][k]);

            dist = sqrt(dist);
            
            if (dist > distmax)
                distmax = dist;
        }
    }

    cout << "\nDistancia maxima: " << distmax;

    // Calculo numero de restricciones maximas
    for (int i = 0; i < restricciones.size(); i++)
        for (int j = i+1; j <= restricciones[i].size(); j++)
            if (restricciones[i][j] == 1 || restricciones[i][j] == -1)
                restrmax++;

    cout << "\nRestricciones maxima: " << restrmax;

    lambda = distmax/restrmax;

    cout << "\nDesviacionGeneral: " << desviacionGeneral(solucion);

    return (desviacionGeneral(solucion) + lambda * infeasibilityBL(solucion));
}

double PAR::distanciaIntracluster (int cluster){
    double suma = 0.0, aux = 0.0;
    int tam = 0;

    for (int i = 0; i < clusters.size(); i++){
        if (clusters[i] == cluster){
            for (int j = 0; j < datos[i].size(); j++)
                aux += (datos[i][j]-centroides[cluster][j]) * (datos[i][j]-centroides[cluster][j]);

            suma += sqrt(aux);
            tam++;
        }
    }

    return suma/tam;
}

double PAR::desviacionGeneral (vector<int> solucion){
    double suma = 0.0;

    for (int i = 0; i < solucion.size(); i++)
        suma += distanciaIntracluster(solucion[i]);

    return suma/solucion.size();
}

double PAR::distancia (int posicionPunto, int cluster){
    double suma = 0.0;

    for (int i = 0; i < datos[posicionPunto].size(); i++)
        suma += (datos[posicionPunto][i] - centroides[cluster][i]) * (datos[posicionPunto][i] - centroides[cluster][i]);

    return sqrt(suma);
}

int PAR::infeasibilityGreedy (vector<int> solucion){
    int infeasibility = 0;

    for (int i = 0; i < restricciones.size(); i++)
        for (int j = i+1; j < restricciones[i].size(); j++){
            if (restricciones[i][j] == 1 && solucion[i] != solucion[j])
                infeasibility++;
            else if (restricciones[i][j] == -1 && solucion[i] == solucion[j])
                infeasibility++;
        }

    return infeasibility;
}


int PAR::infeasibilityGreedy (int posicionPunto, int cluster){
    int infeasibility = 0;

    /*
    // Primera forma: 
    vector<int> aux(clusters); // aux = clusters
    
    aux[posicionPunto] = cluster;

    for (int i = 0; i < restricciones.size(); i++)
        for (int j = i+1; j < restricciones[i].size(); j++){
            if (restricciones[i][j] == 1 && aux[i] != aux[j])
                infeasibility++;
            else if (restricciones[i][j] == -1 && aux[i] == aux[j])
                infeasibility++;   
        }
    */
    // Segunda forma
    
    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] == cluster)
            if (restricciones[posicionPunto][i] == -1)
                infeasibility++;

    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] != cluster)  
            if (restricciones[posicionPunto][i] == 1)
                infeasibility++;
    
    return infeasibility;
}

int PAR::infeasibilityBL (vector<int> solucion){
    int infeasibility = 0;
    list< tuple <int, int, int> >::iterator it = listaRestricciones.begin();

    /*
    for (int i = 0; i < solucion.size(); i++){
        while (get<0>(*it) != i && it != listaRestricciones.end() )
            it++;


    }
    */
    return infeasibility;
}

void PAR::actualizarCentroides (){
    vector<int> contador_clusters(num_clusters, 0);
    vector< vector <double> > aux (centroides); // DUDA BORRAR
    vector<double> aux2;

    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] == -1)
            return;
    
    for (int i = 0; i < clusters.size(); i++)
        contador_clusters[clusters[i]]++;

    for (int i = 0; i < centroides.size(); i++)
        for (int j = 0; j < centroides[i].size(); j++)
            centroides[i][j] = 0.0;
 
    for (int i = 0; i < clusters.size(); i++)
        for (int j = 0; j < centroides[clusters[i]].size(); j++){
            centroides[clusters[i]][j] += datos[i][j];
            contador_clusters[clusters[i]]++;
        }

    for (int k = 0; k < contador_clusters.size(); k++){
        if (contador_clusters[k] == 0){
            centroides.resize(0);

            for (int i = 0; i < num_clusters; i++){
                aux2.clear();
                
                //cerr << "\nCentroide " << i << ": ";
                
                for (int j = 0; j < datos[0].size(); j++){
                    aux2.push_back(Randfloat(0, 1));
                    //cout << " " << aux[j];
                }

                centroides.push_back(aux2);
            }

            return;
        }
    }

    for (int i = 0; i < num_clusters; i++)
        for (int j = 0; j < centroides[i].size(); j++)
            centroides[i][j] = centroides[i][j]/contador_clusters[i];
}

void PAR::cambioCluster (vector<int> solucion, int indice, int cluser){
    // si el indice hace que un cluster se quede vacio,  NO CAMBIAR CLUSTER
}

bool PAR::clusterVacio (vector<int> solucion){
    vector<int> contador_clusters(num_clusters, 0);
    bool estaVacio = false;

    for (int i = 0; i < solucion.size(); i++)
            contador_clusters[solucion[i]]++;

    for (int i = 0; i < contador_clusters.size(); i++)
        if (contador_clusters[i] == 0)
            estaVacio = true;
    
    return estaVacio;
}

vector<int> PAR::greedy (){
    vector<int> indices;
    vector<double> aux;
    vector<int>::iterator it;
    bool hay_cambio_C = true;
    const int dim = datos[0].size();
    int cluster_min, incremento, incremento_min, num_iteraciones = 0; 
    double dist_min, min = 1.0, max = 0.0; // Los datos estan normalizados
    
    // Calculamos el dominio
    for (int i = 0; i < datos.size(); i++)
        for (int j = 0; j < datos[i].size(); j++){
            if (datos[i][j] < min)
                min = datos[i][j];
            if (datos[i][j] > max)
                max = datos[i][j];
        }

    centroides.resize(0);

    for (int i = 0; i < num_clusters; i++){
        aux.clear();
        
        for (int j = 0; j < dim; j++)
            aux.push_back(Randfloat(min, max));

        centroides.push_back(aux);
    }
    
    // Inicializo indices
    for (int i = 0; i < datos.size(); i++)
        indices.push_back(i);

    // Barajo los indices
    random_shuffle(indices.begin(), indices.end()); 

    do{
        hay_cambio_C = false;
        num_iteraciones++;
        
        for (it = indices.begin(); it != indices.end(); ++it){
            cluster_min = clusters[*it];
            incremento_min = 999999999;
            dist_min = 999.0;

            for (int i = 0; i < num_clusters; i++){
                incremento = infeasibilityGreedy(*it, i);

                if (incremento < incremento_min){
                    incremento_min = incremento;
                    cluster_min = i;
                    dist_min = distancia (*it, i);
                }
                else if (incremento == incremento_min){
                    if (distancia (*it, i) < dist_min){
                        cluster_min = i;
                        dist_min = distancia (*it, i);
                    }
                }       
            }

            if (cluster_min != clusters[*it]){
                clusters[*it] = cluster_min;
                hay_cambio_C = true;
            } 
        }
        
        // Actualizar el centroide i promediando las instancias de su cluster asociado i
        actualizarCentroides ();
    } while(hay_cambio_C);

    cout << "\nNumero de iteraciones: " << num_iteraciones << endl;

    return clusters;
}


vector<int> PAR::busquedaLocal (){
    vector<int> indicesDatos, indicesClusters;
    vector<int>::iterator itDatos, itClusters;
    vector<int> solucion(clusters.size(), -1);
    bool recalcular, hay_mejora;
    const int nEvaluacionesMAX = 100000;
    int nEvaluaciones = 0;

    // Genero la solucion inicial
    do{
        recalcular = false;

        for (int i = 0; i < solucion.size(); i++)
            solucion[i] = Randfloat(0, num_clusters);

        // Comprobamos que ningun cluster se queda vacio
        recalcular = clusterVacio (solucion);
    } while(recalcular);

    double fit = fitness(solucion);

    cout << "\n\nFIT: " << fit ;
    /*

    // GENERACION DE VECINOS

    // Inicializo indices a los datos
    for (int i = 0; i < datos.size(); i++)
        indicesDatos.push_back(i);

    // Barajo los indices a los datos
    random_shuffle(indicesDatos.begin(), indicesDatos.end()); 

    // Inicializo indices a los clusters
    for (int i = 0; i < num_clusters; i++)
        indicesClusters.push_back(i);

    // Barajo los indices a los clusters
    random_shuffle(indicesClusters.begin(), indicesClusters.end()); 

    do{
        hay_mejora = false;

        for (itDatos = indicesDatos.begin(); itDatos != indicesDatos.end() && !hay_mejora; ++itDatos){
            for (itClusters = indicesClusters.begin(); itClusters != indicesClusters.end() && !hay_mejora; ++itClusters){
                if (*itClusters != clusters[*itDatos]){

                }

            }
        
        }
    } while (hay_mejora && nEvaluaciones < nEvaluacionesMAX);
    */
    return solucion;
}