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
            cout << i << " " << j << " :" << restricciones[i][j] << "\n";
        
        cout << "\n";
    }

    cout <<"\n\n\nLISTA ML:\n";


    for (int i = 0; i < restriccionesML.size(); i++){   
        for (int j = 0; j < restriccionesML[i].size(); j++)
                cout << i << " " << j << " :" << restriccionesML[i][j] << "\n";
        
        cout << "\n";
    }

    cout <<"\n\n\nLISTA CL:\n";


    for (int i = 0; i < restriccionesCL.size(); i++){   
        for (int j = 0; j < restriccionesCL[i].size(); j++)
                cout << i << " " << j << " :" << restriccionesCL[i][j] << "\n";
        
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

        for (int j = 0; j < distancias[i].size(); j++)
            cout << distancias[i][j] << " ";
    }
}

PAR::PAR (const string fDatos, const string fRestricciones, const int k){
    vector<int> aux2;

    num_clusters = k;

    leerFicheros (fDatos, fRestricciones);
    leerFicheros2 (fDatos, fRestricciones);
    
    clusters.resize(datos.size(), -1);

    // Inicializo lista de restricciones
    for (int i = 0; i < restricciones.size(); i++){
        for (int j = i+1; j < restricciones[i].size(); j++){
            aux2.clear();
            aux2.push_back(i);
            aux2.push_back(j);
            
            if (restricciones[i][j] == 1)
                restriccionesML.push_back(aux2);
            else if (restricciones[i][j] == -1)
                restriccionesCL.push_back(aux2);
        }
    }

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

        
}

void PAR::leerFicheros2 (const string fDatos, const string fRestricciones){
        string linea, valor;
    ifstream ifDatos(fDatos), ifRestricciones(fRestricciones);
    int contador = 0;

    // Cargo el fichero de restricciones
    if (ifRestricciones.is_open()){
        restricciones.resize(contador);

        while (getline(ifRestricciones, linea)){
            stringstream cadena(linea);
            restricciones.resize(contador + 1);

            while (getline(cadena, valor, ',')){
           
                restricciones[contador].push_back(stod(valor));
                
            }
          
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
    double agregado;

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
    
    // Calculo numero de restricciones maximas
    for (int i = 0; i < restricciones.size(); i++)
        for (int j = i+1; j < restricciones[i].size(); j++)
            if (restricciones[i][j] == 1 || restricciones[i][j] == -1)
                restrmax++;
    
    lambda = distmax/restrmax;
    agregado = (desviacionGeneral(solucion) + lambda * infeasibilityGreedy(solucion));
    /*
    cout << "\nLambda: " << lambda;
    cout << "\nDistancia maxima: " << distmax;
    cout << "\nRestricciones maxima: " << restrmax;
    cout << "\nError distancias: " << desviacionGeneral(solucion);
    cout << "\nAgregado: " << agregado;
    */
    return agregado;
}

double PAR::distanciaIntracluster (int cluster, vector<int> solucion){
    double suma = 0.0, aux = 0.0;
    int tam = 0;

    for (int i = 0; i < solucion.size(); i++){
        aux = 0.0;

        if (solucion[i] == cluster){
            for (int j = 0; j < datos[i].size(); j++){
                aux += (datos[i][j]-centroides[cluster][j]) * (datos[i][j]-centroides[cluster][j]);
            }

            suma += sqrt(aux);
            tam++;
        }
    }

    return suma/tam;
}

double PAR::desviacionGeneral (vector<int> solucion){
    double suma = 0.0;

    for (int i = 0; i < num_clusters; i++)
        suma += distanciaIntracluster(i, solucion);

    return suma/num_clusters;
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
            if (restricciones[posicionPunto][i] == -1){
                //cout << "\nrestriccion -1:" << posicionPunto << " " << i << " " << clusters[i] << " " << cluster << " " << restricciones[posicionPunto][i];
                infeasibility++;
            }

    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] != cluster && clusters[i] != -1)  
            if (restricciones[posicionPunto][i] == 1 && posicionPunto != i){
                //cout << "\nrestriccion 1:" << posicionPunto << " " << i;
                infeasibility++;
            }
            
    return infeasibility;
}

int PAR::infeasibilityBL (vector<int> solucion){
    int infeasibility = 0;

    for (int i = 0; i < restriccionesML.size(); i++)
        if (solucion[restriccionesML[i][0]] != solucion[restriccionesML[i][1]])
            infeasibility++;

    for (int i = 0; i < restriccionesCL.size(); i++)
        if (solucion[restriccionesCL[i][0]] == solucion[restriccionesCL[i][1]])
            infeasibility++;
    
    return infeasibility;
}

void PAR::actualizarCentroides (){
    vector<int> contador_clusters(num_clusters, 0);
    vector< vector <double> > aux (centroides); // DUDA BORRAR
    vector<double> aux2;
    vector<double> aux3;

    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] == -1){
            cerr << "\n\nESTE MENSAJE NO DEBERIA SALIR\n\n";
            return;
        }
    
    for (int i = 0; i < clusters.size(); i++)
        contador_clusters[clusters[i]]++;
    
    centroides.resize(num_clusters);

    for (int i = 0; i < centroides.size(); i++)
        centroides[i].resize(datos[0].size(), 0.0);
 
    for (int i = 0; i < clusters.size(); i++)
        for (int j = 0; j < centroides[clusters[i]].size(); j++){
            centroides[clusters[i]][j] += datos[i][j];
            //cerr << "\n\tit: " << i << " problema.datos[(*it)][j]: " << datos[i][j];
            //contador_clusters[clusters[i]]++;
        }
//cout << "\nPFFFFFFFF2" << endl;
    for (int i = 0; i < num_clusters; i++)
        for (int j = 0; j < centroides[i].size(); j++){
            centroides[i][j] = centroides[i][j]/contador_clusters[i]*1.0;
            //cout << " " << centroides[i][j];
        }
    
//cout << endl;
}

bool PAR::cambioCluster (vector<int> & solucion, int indice, int cluster){
    // si el indice hace que un cluster se quede vacio,  NO CAMBIAR CLUSTER
    if (parValido(make_pair(indice, cluster), solucion)){
        solucion[indice] = cluster;
        
        return true;
    }
    else
        return false;
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

    //cout << "\nINICIO:";

    //imprimirCentroides();
    
    // Inicializo indices
    for (int i = 0; i < datos.size(); i++)
        indices.push_back(i);

    // Barajo los indices
    random_shuffle(indices.begin(), indices.end()); 

    do{
        hay_cambio_C = false;
        num_iteraciones++;
        
        for (it = indices.begin(); it != indices.end(); ++it){
            //cout << "\nPOSICION " << *it << " ESTA EN " << clusters[*it];
            cluster_min = clusters[*it];
            incremento_min = 999999999;
            dist_min = 999.0;

            for (int i = 0; i < num_clusters; i++){
                incremento =  infeasibilityGreedy(*it, i);
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
	            //cout << "\nCAMBIO: "  << "POSICION " << *it << " RESTRICCIONES " << incremento_min << " DISTANCIA " << dist_min << " CLUSTER " << cluster_min;
                clusters[*it] = cluster_min;
                hay_cambio_C = true;
            } 
        }
        
        // Actualizar el centroide i promediando las instancias de su cluster asociado i
        //cout << "\n\nCAMBIO: ACTUALIZO CENTROIDES";
        actualizarCentroides ();
        //imprimirCentroides();
    } while(hay_cambio_C);

    cout << "\nNumero de iteraciones: " << num_iteraciones << endl;

    return clusters;
}

bool PAR::parValido (pair<int, int> par, vector<int> solucion){
    vector<int> aux(solucion);

    aux[par.first] = par.second;
    
    return !clusterVacio(aux);
}

// DUDA REPITO CODIGO
double PAR::fitnessBL (vector<int> solucion){
    double agregado = 0.0, dist, distmax = 0.0, lambda;
    int restrmax;
    
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
    
    // Calculo numero de restricciones maximas
    restrmax = restriccionesML.size() + restriccionesCL.size();
    
    lambda = distmax/restrmax;
    agregado = (desviacionGeneral(solucion) + lambda * infeasibilityBL(solucion));
    //cout << "\ndesviacionGeneral(solucion): "<< desviacionGeneral(solucion);
    //cout << "\nlambda: " << lambda;
    //cout << "\ninfeasibilityBL(solucion): " << infeasibilityBL(solucion);

    return agregado;
}


vector<int> PAR::busquedaLocal (){
    vector<int> indicesDatos, indicesVecindarios;
    vector<int>::iterator itDatos, itVecindarios;
    vector<int> solucion(clusters.size());
    bool recalcular, hay_mejora;
    const int nEvaluacionesMAX = 100000;
    int nEvaluaciones = 0;
    vector< pair<int, int> > vecindarioVirtual;
    pair<int, int> par;
    double fit_min, fit;

    // Genero la solucion inicial
    do{
        recalcular = false;

        for (int i = 0; i < clusters.size(); i++)
            clusters[i] = Randint(0, num_clusters-1);

        // Comprobamos que ningun cluster se queda vacio
        recalcular = clusterVacio (clusters);
    } while(recalcular);
 
    actualizarCentroides();

    vector<int> nuevaSolucion(clusters);
    
    // GENERACION DE VECINOS

    // Inicializo indices a los datos
    for (int i = 0; i < datos.size(); i++)
        indicesDatos.push_back(i);

    // Barajo los indices a los datos
    random_shuffle(indicesDatos.begin(), indicesDatos.end()); 

    // Vecindario virtual
    for (int i = 0; i < clusters.size(); i++)
        for (int j = 0; j < num_clusters; j++)
            if (clusters[i] != j){
                par = make_pair(i, j);

                if (parValido(par, clusters)){
                    vecindarioVirtual.push_back(par);
                    //cout << "\nCLUSTER: " << clusters[i] << " par (" << par.first << "," << par.second << ")";
                }
            }

    // Inicializo indices al vecindario
    for (int i = 0; i < vecindarioVirtual.size(); i++)
        indicesVecindarios.push_back(i);

    // Barajo los indices a los clusters
    random_shuffle(indicesVecindarios.begin(), indicesVecindarios.end()); 

    double agregado = 0.0, dist, distmax = 0.0, lambda;
    int restrmax;
    
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
    restrmax = restriccionesML.size() + restriccionesCL.size();
    
    lambda = distmax/restrmax;

    //imprimirCentroides();
    fit_min = fitnessBL(clusters);
    cout << "\nINICIO fitness min: " << fit_min;
    cout << "\ndesviacionGeneral(solucion): " << desviacionGeneral(clusters);
    cout << "\nlambda: " << lambda;
    cout << "\ninfeasibilityBL(solucion): " << infeasibilityBL(clusters);

    cout << "\n\n\nAAAAAAAAAAAAANTES LISTA CLUSTERS:";

    for (int i = 0; i < clusters.size(); i++)
        cout << " " << clusters[i];
     cout << "\n" << infeasibilityBL(clusters);
     
    do{
        random_shuffle(indicesDatos.begin(), indicesDatos.end()); 
        random_shuffle(indicesVecindarios.begin(), indicesVecindarios.end()); 

        hay_mejora = false;
        //fit_min = fitnessBL(clusters);
        nEvaluaciones++;

        for (itDatos = indicesDatos.begin(); itDatos != indicesDatos.end() && !hay_mejora; ++itDatos){
            for (itVecindarios = indicesVecindarios.begin(); itVecindarios != indicesVecindarios.end() && !hay_mejora; ++itVecindarios){
                if (vecindarioVirtual[*itVecindarios].first == *itDatos){
                    //cout << "\nCOMPARO vecindarioVirtual[*itVecindarios].first " << vecindarioVirtual[*itVecindarios].first << " con *itDatos " << *itDatos << " vecindarioVirtual[*itVecindarios].second " << vecindarioVirtual[*itVecindarios].second;
                    
                    if (cambioCluster(nuevaSolucion, *itDatos, vecindarioVirtual[*itVecindarios].second)){
                        cambioCluster(clusters, *itDatos, vecindarioVirtual[*itVecindarios].second);
                        actualizarCentroides();
                        fit = fitnessBL(nuevaSolucion);
                        nEvaluaciones++;

                        //cout << "\nCAMBIO CLUSTER: indice " << *itDatos << " cluster " << vecindarioVirtual[*itVecindarios].second << " fitness " << fit << " fitness_min " << fit_min << " cluster_min " << clusters[vecindarioVirtual[*itVecindarios].first] << " ?????? " << vecindarioVirtual[*itVecindarios].first;


                        if (fit < fit_min){
                            //cout << "\nHAY MEJORA: " << fit << " < " << fit_min;
                            fit_min = fit;
                            hay_mejora = true;
                            //if (nEvaluaciones % 100 == 0)
                                //cerr << " " << nEvaluaciones;
                            //cambioCluster(clusters, *itDatos, vecindarioVirtual[*itVecindarios].second);
                            //actualizarCentroides();
                            //imprimirCentroides();
                        }
                        else{
                            cambioCluster(clusters, *itDatos, clusters[vecindarioVirtual[*itVecindarios].first]);
                            actualizarCentroides();
                        }
                        
                    }
                    
                }
                

            }
        
        }
    } while (hay_mejora && nEvaluaciones < nEvaluacionesMAX);
    
    cout << endl << nEvaluaciones;

    cout << "\n\n\nDESPUEEEEEEEES LISTA CLUSTERS:";

    for (int i = 0; i < clusters.size(); i++)
        cout << " " << clusters[i];

    cout << "\n" << infeasibilityBL(clusters);
    
    return clusters;
}