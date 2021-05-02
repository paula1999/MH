#include <iostream>
#include "PAR.h"
#include <fstream>
#include <sstream>
#include "random.h"
#include <algorithm>
#include <math.h>
#include <tuple>

using namespace std;

ostream & operator << (ostream & os, const PAR & par){
    os << "\nTABLA: \t (infeas. \t | \t desvGeneral \t | \t fitness \t ) \n";
    os << "\t\t\t" << par.infeasibility << " \t\t & \t\t " << par.desvGeneral << " \t & \t " << par.fitness << " \t &";

    os << "\nLista de asignaciones a los clusters:";

    for (int i = 0; i < par.clusters.size(); i++)
        os << " " << par.clusters[i];
            
    os << "\n\nInfeasable: " << par.infeasibility;
    os << "\n\nAgregado: " << par.fitness;
    os << "\n\nRestricciones maximas: " << par.restrMax;
    os << "\n\nDistancia maxima: " << par.distMax;

    os << "\n\nCLUSTERS:\n";
    
    for (int i = 0; i < par.num_clusters; i++){
        os << "\nCluster " << i << ": ";
        
        for (int j = 0; j < par.clusters.size(); j++)
            if (par.clusters[j] == i)
                os << j << " ";
    }

    os << "\n\nCENTROIDES:\n";
    
    for (int i = 0; i < par.centroides.size(); i++){
        os << "\nCentroide " << i << ": ";
        
        for (int j = 0; j < par.centroides[i].size(); j++)
            os << par.centroides[i][j] << " ";
    }

    os << "\n\nDESVIACION GENERAL: " << par.desvGeneral;
    //os << "\n\nError_Dist: " << abs(par.desvGeneral - par.distOptima);

    return os;
}

void PAR::imprimirRestricciones () const{
    cout <<"\n\n\nMATRIZ:\n";
    
    for (int i = 0; i < restricciones.size(); i++){
        for (int j = 0; j < restricciones[i].size(); j++)
            cout << restricciones[i][j] << " ";
        
        cout << "\n";
    }

    cout <<"\n\n\nLISTA ML:\n";

    for (int i = 0; i < restriccionesML.size(); i++)
        cout << "(" << restriccionesML[i][0] << "," << restriccionesML[i][1] << "), ";

    cout <<"\n\n\nLISTA CL:\n";

    for (int i = 0; i < restriccionesCL.size(); i++)
        cout << "(" << restriccionesCL[i][0] << "," << restriccionesCL[i][1] << "), ";
}

void PAR::imprimirDatos () const{
    for (int i = 0; i < datos.size(); i++){
        for (int j = 0; j < datos[i].size(); j++)
            cout << datos[i][j] << " ";
        
        cout << "\n";
    }
}

void PAR::imprimirCentroides () const{
    for (int i = 0; i < centroides.size(); i++){
        cout << "\nCentroide " << i << ": ";
        
        for (int j = 0; j < centroides[i].size(); j++)
            cout << centroides[i][j] << " ";
    }
}

void PAR::imprimirClusters () const{
    for (int i = 0; i < num_clusters; i++){
        cout << "\nCluster " << i << ": ";
        
        for (int j = 0; j < clusters.size(); j++)
            if (clusters[j] == i)
                cout << j << " ";
    }
}

void PAR::imprimirDistancias () const{
    for (int i = 0; i < distancias.size(); i++){
        cout << "\n";

        for (int j = 0; j < distancias[i].size(); j++)
            cout << distancias[i][j] << " ";
    }
}

PAR::PAR (const string fDatos, const string fRestricciones, const int k){
    vector<int> aux2;

    // Inicializo la distancia óptima
    if (fDatos.find("zoo") != string::npos)
        distOptima = 0.904799856193481;
    else if (fDatos.find("glass") != string::npos)
        distOptima = 0.364290281975566;
    else if (fDatos.find("bupa") != string::npos)
        distOptima = 0.220423749236421;
    else
        distOptima = 0;

    num_clusters = k;

    leerFicheroDatos (fDatos);
    leerFicheroRestricciones (fRestricciones);
    
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

    calcularDistMax();
    
    restrMax = restriccionesML.size() + restriccionesCL.size(); // Calculo numero de restricciones maximas
    lambda = distMax/restrMax;
    infeasibility = -1;
    desvGeneral = -1;
    fitness = -1;
}

void PAR::calcularDistancias (){
    double dist = 0.0;

    for (int i = 0; i < datos.size(); i++)
        for (int j = i+1; j < datos.size(); j++){
            dist = 0.0;
            
            for (int k = 0; k < datos[i].size(); k++)
                dist += (datos[i][k] - datos[j][k]) * (datos[i][k] - datos[j][k]);

            distancias[i][j] = sqrt(dist);
        }
}

void PAR::leerFicheroDatos (const string fDatos){
    string linea, valor;
    ifstream ifDatos(fDatos);
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
}

void PAR::leerFicheroRestricciones (const string fRestricciones){
    string linea, valor;
    ifstream ifRestricciones(fRestricciones);
    int contador = 0;

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

void PAR::calcularDistMax(){
    double dist = 0.0, dMax = 0.0;

    // Calculo mayor distancia
    for (int i = 0; i < datos.size(); i++){
        for (int j = i+1; j < datos.size(); j++){
            dist = 0.0;
            
            for (int k = 0; k < datos[i].size(); k++)
                dist += (datos[i][k] - datos[j][k]) * (datos[i][k] - datos[j][k]);

            dist = sqrt(dist);
            
            if (dist > dMax)
                dMax = dist;
        }
    }

    distMax = dMax;
}


double PAR::fitnessGreedy (){
    fitness = desviacionGeneral() + lambda * infeasibilityGreedy();
    
    return fitness;
}

double PAR::distanciaIntracluster (int cluster, const vector<int> solucion){
    double suma = 0.0, aux = 0.0;
    int tam = 0;

    for (int i = 0; i < solucion.size(); i++){
        aux = 0.0;

        if (solucion[i] == cluster){
            for (int j = 0; j < datos[i].size(); j++)
                aux += (datos[i][j]-centroides[cluster][j]) * (datos[i][j]-centroides[cluster][j]);

            suma += sqrt(aux);
            tam++;
        }
    }

    return suma/tam;
}

double PAR::desviacionGeneral (){
    double suma = 0.0;

    actualizarCentroides();
    
    for (int i = 0; i < num_clusters; i++)
        suma += distanciaIntracluster(i, clusters);

    desvGeneral = suma/num_clusters;

    return desvGeneral;
}

double PAR::distancia (int posicionPunto, int cluster) const{
    double suma = 0.0;

    for (int i = 0; i < datos[posicionPunto].size(); i++)
        suma += (datos[posicionPunto][i] - centroides[cluster][i]) * (datos[posicionPunto][i] - centroides[cluster][i]);

    return sqrt(suma);
}

int PAR::infeasibilityGreedy () const{
    int infeas = 0;

    for (int i = 0; i < restricciones.size(); i++)
        for (int j = i+1; j < restricciones[i].size(); j++){
            if (restricciones[i][j] == 1 && clusters[i] != clusters[j])
                infeas++;
            else if (restricciones[i][j] == -1 && clusters[i] == clusters[j])
                infeas++;
        }
    
    return infeas;
}

int PAR::infeasibilityCluster (int posicionPunto, int cluster) const{
    int infeas = 0;
    
    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] == cluster)
            if (restricciones[posicionPunto][i] == -1)
                infeas++;

    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] != cluster && clusters[i] != -1)  
            if (restricciones[posicionPunto][i] == 1 && posicionPunto != i)
                infeas++;
            
    return infeas;
}

int PAR::infeasibilityBL () const{
    int infeas = 0;
    
    for (int i = 0; i < restriccionesML.size(); i++)
        if (clusters[restriccionesML[i][0]] != clusters[restriccionesML[i][1]])
            infeas++;

    for (int i = 0; i < restriccionesCL.size(); i++)
        if (clusters[restriccionesCL[i][0]] == clusters[restriccionesCL[i][1]])
            infeas++;
    
    return infeas;
}

int PAR::infeasibilityBL (int posicionPunto, int cluster) const{
    int infeas = 0;

    for (int i = 0; i < restriccionesML.size(); i++){
        if (restriccionesML[i][0] == posicionPunto && cluster != clusters[restriccionesML[i][1]])
            infeas++;
        else if (restriccionesML[i][1] == posicionPunto && cluster != clusters[restriccionesML[i][0]])
            infeas++;
    }

    for (int i = 0; i < restriccionesCL.size(); i++){
        if (restriccionesCL[i][0] == posicionPunto && cluster == clusters[restriccionesCL[i][1]])
            infeas++;
        else if (restriccionesCL[i][1] == posicionPunto && cluster == clusters[restriccionesCL[i][0]])
            infeas++;
    }
    
    return infeas;
}

void PAR::actualizarCentroides (){
    vector<int> contador_clusters(num_clusters, 0);

    for (int i = 0; i < clusters.size(); i++)
        if (clusters[i] == -1)
            return;
    
    for (int i = 0; i < clusters.size(); i++)
        contador_clusters[clusters[i]]++;

    centroides.resize(0);
    centroides.resize(num_clusters);

    for (int i = 0; i < centroides.size(); i++)
        centroides[i].resize(datos[0].size(), 0.0);
 
    for (int i = 0; i < clusters.size(); i++)
        for (int j = 0; j < centroides[clusters[i]].size(); j++)
            centroides[clusters[i]][j] += datos[i][j];

    for (int i = 0; i < num_clusters; i++)
        for (int j = 0; j < centroides[i].size(); j++)
            centroides[i][j] = centroides[i][j]/contador_clusters[i]*1.0;
}

bool PAR::cambioCluster (vector<int> & solucion, int indice, int cluster){
    if (parValido(make_pair(indice, cluster), solucion)){
        solucion[indice] = cluster;
        
        return true;
    }
    else
        return false;
}

bool PAR::clusterVacio (vector<int> solucion) const{
    vector<int> contador_clusters(num_clusters, 0);
    bool estaVacio = false;
    
    for (int i = 0; i < solucion.size(); i++){
        //cout << solucion[i] << " ";
        contador_clusters[solucion[i]]++;
    }
    
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
    int cluster_min, incremento, incremento_min; 
    double dist_min, min = 0.0, max = 1.0; // Los datos estan normalizados

    centroides.clear();

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
        
        for (it = indices.begin(); it != indices.end(); ++it){
            cluster_min = clusters[*it];
            incremento_min = 999999999;
            dist_min = 999.0;

            for (int i = 0; i < num_clusters; i++){
                incremento =  infeasibilityCluster(*it, i);

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

    

    
    if (!clusterVacio(clusters)){
        desvGeneral = desviacionGeneral();
        fitness = fitnessGreedy();
        infeasibility = infeasibilityGreedy();
        
        return clusters;
    }
    else{
        return greedy();
    }
}

bool PAR::parValido (pair<int, int> par, vector<int> solucion) const{
    vector<int> aux(solucion);

    aux[par.first] = par.second;
    
    return !clusterVacio(aux);
}

double PAR::fitnessBL (){
    fitness = desviacionGeneral() + lambda * infeasibilityBL();

    return fitness;
}


vector<int> PAR::busquedaLocal (){
    vector<int> indicesDatos, indicesVecindarios;
    vector<int>::iterator itDatos, itVecindarios;
    bool recalcular, hay_mejora, nuevainf;
    const int nEvaluacionesMAX = 100000;
    int nEvaluaciones = 0, antiguainf, infeasibility_min, infeasibility_nueva, cluster_min;
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
    
    // Inicializo indices a los datos
    for (int i = 0; i < datos.size(); i++)
        indicesDatos.push_back(i);    

    fit_min = fitnessBL();
    infeasibility_min = infeasibilityBL();
    infeasibility_nueva = infeasibility_min;

    do{
        random_shuffle(indicesDatos.begin(), indicesDatos.end());

        hay_mejora = false;

        for (itDatos = indicesDatos.begin(); itDatos != indicesDatos.end() && !hay_mejora && nEvaluaciones < nEvaluacionesMAX; ++itDatos){
            // Vecindario virtual
            vecindarioVirtual.clear();

            for (int j = 0; j < num_clusters; j++)
                if (clusters[*itDatos] != j){
                    par = make_pair(*itDatos, j);

                    if (parValido(par, clusters))
                        vecindarioVirtual.push_back(par);
                }

            // Inicializo indices al vecindario
            indicesVecindarios.clear();
            
            for (int i = 0; i < vecindarioVirtual.size(); i++)
                indicesVecindarios.push_back(i);

            random_shuffle(indicesVecindarios.begin(), indicesVecindarios.end());
            
            for (itVecindarios = indicesVecindarios.begin(); itVecindarios != indicesVecindarios.end() && !hay_mejora && nEvaluaciones < nEvaluacionesMAX; ++itVecindarios){                    
                cluster_min = clusters[*itDatos];
                clusters[*itDatos] = vecindarioVirtual[*itVecindarios].second;
                infeasibility_nueva -= infeasibilityCluster(*itDatos, cluster_min);
                infeasibility_nueva += infeasibilityCluster(*itDatos, vecindarioVirtual[*itVecindarios].second);
                fit = desviacionGeneral() + lambda * infeasibility_nueva;
                nEvaluaciones++;

                if (fit < fit_min){
                    fit_min = fit;
                    hay_mejora = true;   
                    infeasibility_min = infeasibility_nueva;     
                }
                else{
                    clusters[*itDatos] = cluster_min;
                    infeasibility_nueva = infeasibility_min;
                }  
            }
        }
    } while (hay_mejora && nEvaluaciones < nEvaluacionesMAX);
    
    fitness = fitnessBL();
    infeasibility = infeasibilityBL();

    return clusters;
}



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// PRÁCTICA 2 //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


int PAR::busquedaLocalSuave (vector<int> & solucion, double & fit_min, const int nFallosMAX){
    vector<int> indicesDatos;
    const int n = datos.size();
    int nFallos, i, cluster_min, infeasibility_min, infeasibility_nueva, nEvaluaciones = 0;
    bool mejora;
    pair<int, int> par;
    double fit;

    // Inicializo indices a los datos
    for (i = 0; i < n; i++)
        indicesDatos.push_back(i); 

    random_shuffle(indicesDatos.begin(), indicesDatos.end());

    nFallos = 0;
    i = 0;
    mejora = true;
    infeasibility_min = infeasibilityBL();
    infeasibility_nueva = infeasibility_min;

    while ((mejora or nFallos < nFallosMAX) and i < n) {
        mejora = false;

        // Asignar el mejor valor posible a S_RSI[i]
        for (int j = 0; j < num_clusters; j++){
            par = make_pair(solucion[indicesDatos[i]], j);

            if (solucion[indicesDatos[i]] != j and parValido(par, solucion)){
                cluster_min = solucion[indicesDatos[i]];
                solucion[indicesDatos[i]] = j;
                infeasibility_nueva -= infeasibilityCluster(indicesDatos[i], cluster_min);
                infeasibility_nueva += infeasibilityCluster(indicesDatos[i], j);
                fit = desviacionGeneral() + lambda * infeasibility_nueva;
                nEvaluaciones++;

                if (fit < fit_min){
                    fit_min = fit;
                    mejora = true;   
                    infeasibility_min = infeasibility_nueva;     
                }
                else{
                    solucion[indicesDatos[i]] = cluster_min;
                    infeasibility_nueva = infeasibility_min;
                }
            }
        }

        if (!mejora)
            nFallos++;

        i++;
    }

    return nEvaluaciones;
}


vector<int> PAR::algoritmoGenetico (int M, const string evolucion, const string operadorCruce, const double probCruce){
    bool recalcular;
    double fit, fit_min, probMutacion, prob;
    const int nFitnessMAX = 100000;
    int m, n = clusters.size(), nCruces, nMutaciones, pos, c1, c2, nIterMAX, posCMejor, nFitness;
    vector< vector<int> > pActual, pSiguiente;
    vector<double> pActualFitness;
    vector<int> cromosoma, CMejor;

    probMutacion = 0.1/n;

    // TO DO: inicializar parametros para no tener tantos ifs distinguiendo evoluciones y cruces
    // TO DO: inicializar parametros para no tener tantos ifs distinguiendo evoluciones y cruces
    if (evolucion.compare("G") == 0)
        m = M;
    else if (evolucion.compare("E") == 0)
        m = 2;
    else{
        cerr << "ERROR. Parametro 'evolucion' incorrecto\n";
        exit(-1);
    }
    
    // Inicializar y evaluar la poblacion actual
    for (int i = 0; i < M; i++){
        // Calcular un cromosoma valido
        do{
            cromosoma.clear();

            for (int j = 0; j < n; j++)
                cromosoma.push_back(Randint(0, num_clusters-1));

            // Comprobar que ningun cluster se queda vacio
            recalcular = clusterVacio (cromosoma);
        } while(recalcular);

        pActual.push_back(cromosoma);
    }

    // Evaluar la poblacion actual
    // TO DO crear mejor una funcion para hacer esto
    for (int i = 0; i < M; i++){
        clusters = pActual[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
    }

    nFitness = 0;
    
    for (int t = 0; nFitness < nFitnessMAX; t++){
        
        // Guardo el mejor cromosoma
        posCMejor = calcularMejorCromosoma (pActualFitness);
        CMejor = pActual[posCMejor];
        //cerr << "Mejor cromosoma: " << pActualFitness[posCMejor] << endl;
        //cerr << nFitness << " ";
        ////////////////////////////////// SELECCIÓN //////////////////////////////////

        pSiguiente.clear();

        // Seleccionar P' desde P(t-1)
        for (int i = 0; i < m; i++){
            c1 = Randint(0, M-1);

            do{
                c2 = Randint(0, M-1);
            } while (c1 == c2);

            //cerr << "\n" << c1 << " " << c2;

            pos = operadorSeleccion(c1, c2, pActualFitness);
            //cerr << " " << pActualFitness[pos];

            pSiguiente.push_back(pActual[pos]);
        }


        /////////////////////////////////////
        /*
        cerr << "Seleccion acabada" << endl;  
        cerr << "---------------------" << endl;
        cerr << "\nPADRES\n++++++++++++++++++" << endl;
        cerr << "Tamaño poblacion: " << m << endl;
        cerr << "Mejor solucion: " << posCMejor << ", mejor objetivo: " << pActualFitness[posCMejor] <<  endl;
        

        pActualFitness.clear();

        for (int i = 0; i < m; i++){
        clusters = pSiguiente[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
         }
        cerr << "\nHIJOS\n++++++++++++++++++" << endl;
        cerr << "Tamaño poblacion: " << m << endl;
        cerr << "Mejor solucion: " << calcularMejorCromosoma (pActualFitness) << ", mejor objetivo: " << pActualFitness[calcularMejorCromosoma (pActualFitness)] <<  endl;
        
    pActualFitness.clear();
for (int i = 0; i < M; i++){
        clusters = pActual[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
    }
        ////////////////////////////////////////
*/

        ////////////////////////////////// CRUZAR //////////////////////////////////

        // Recombinar P'
        nCruces = probCruce * m/2.0;
        //cerr << nCruces << " ";

        // FIXME esto esta feo pero no se como arreglarlo todavia
        for (int i = 0; i < 2*nCruces; i += 2){
            if (operadorCruce.compare("UN") == 0){
                pSiguiente.push_back(operadorCruceUN(pSiguiente[i], pSiguiente[i+1]));           
                pSiguiente.push_back(operadorCruceUN(pSiguiente[i+1], pSiguiente[i]));
                
            }
            else if (operadorCruce.compare("SF") == 0){
                pSiguiente.push_back(operadorCruceSF(pSiguiente[i], pSiguiente[i+1]));
                pSiguiente.push_back(operadorCruceSF(pSiguiente[i+1], pSiguiente[i]));
            }
            else{
                cerr << "ERROR. Parametro 'operador cruce' incorrecto\n";
                exit(-1);
            }
        }
        
        for (int i = 0; i < 2*nCruces; i++)
            pSiguiente.erase(pSiguiente.begin());

        //cerr << "Cruzar acabada" << endl;  
        //cerr << "---------------------" << endl;
/*
        pActualFitness.clear();

        for (int i = 0; i < m; i++){
        clusters = pSiguiente[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
         }
        cerr << "\nHIJOS\n++++++++++++++++++" << endl;
        cerr << "Tamaño poblacion: " << m << endl;
        cerr << "Mejor solucion: " << calcularMejorCromosoma (pActualFitness) << ", mejor objetivo: " << pActualFitness[calcularMejorCromosoma (pActualFitness)] <<  endl;
        
    pActualFitness.clear();
for (int i = 0; i < M; i++){
        clusters = pActual[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
    }
*/
        ////////////////////////////////// MUTACIÓN //////////////////////////////////

        // Mutar P'
        if (evolucion.compare("G") == 0){
            nMutaciones = probMutacion * m * n;

            for (int i = 0; i < nMutaciones; i++){
                pos = Randint(0, m-1);
                pSiguiente[pos] = operadorMutacionUN (pSiguiente[pos]);
            }
        }
        else if (evolucion.compare("E") == 0){
            for (int i = 0; i < pSiguiente.size(); i++){
                prob = Randfloat(0.0, 1.0);

                if (prob < probMutacion*n){
                    pos = Randint(0, m-1);
                    pSiguiente[pos] = operadorMutacionUN (pSiguiente[pos]);
                }
            }
        }
/*
        cerr << "Mutar acabada" << endl;  
        cerr << "---------------------" << endl;

        pActualFitness.clear();

        for (int i = 0; i < m; i++){
        clusters = pSiguiente[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
         }
        cerr << "\nHIJOS\n++++++++++++++++++" << endl;
        cerr << "Tamaño poblacion: " << m << endl;
        cerr << "Mejor solucion: " << calcularMejorCromosoma (pActualFitness) << ", mejor objetivo: " << pActualFitness[calcularMejorCromosoma (pActualFitness)] <<  endl;
        
    pActualFitness.clear();
for (int i = 0; i < M; i++){
        clusters = pActual[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
    }
        
*/
        ////////////////////////////////// REEMPLAZAR ////////////////////////////////// 
        if (evolucion.compare("G") == 0){
            // Reemplazar P(t) a partir de P(t-1) y P'
            pActual = pSiguiente;

            // Si el mejor cromosoma no esta, lo añado al final
            if (find(pActual.begin(), pActual.end(), CMejor) == pActual.end())
                pActual[M-1] = CMejor;
            
        }
        // DUDA esto ponerlo mas bonito
        // FIXME nombres de variables contradictorios
        else if (evolucion.compare("E") == 0){
            for (int i = 0; i < m; i++){
                posCMejor = calcularPeorCromosoma (pActualFitness);
                clusters = pSiguiente[i];
                fit = fitnessBL();
                nFitness++;


                if (pActualFitness[posCMejor] > fit){
                    pActual[posCMejor] = pSiguiente[i];
                    pActualFitness[posCMejor] = fit;
                }
            }
        }

    /*
        cerr << "\njjPADRES\n++++++++++++++++++" << endl;
        cerr << "Tamaño poblacion: " << m << endl;
        cerr << "Mejor solucion: " << posCMejor << ", mejor objetivo: " << pActualFitness[posCMejor] <<  endl;
        */

        ////////////////////////////////// EVALUAR //////////////////////////////////
        
        // Evaluar P(t)
        // TO DO crear mejor una funcion para hacer esto
        if (evolucion.compare("G") == 0){
            pActualFitness.clear();
            
            for (int i = 0; i < M; i++){
                clusters = pActual[i];
                pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
                nFitness++;
            }
            
        }

    }

    // Actualizo atributos
    // TO DO: mejorar esto
    posCMejor = calcularMejorCromosoma (pActualFitness);
    clusters = pActual[posCMejor];
    desvGeneral = desviacionGeneral();
    fitness = fitnessGreedy();
    infeasibility = infeasibilityGreedy();

    return clusters;
}

// usar AGG-UN
vector<int> PAR::algoritmoMemetico (const int M, const string probSeleccion, const int nGeneracion, const double probCruce){
    bool recalcular;
    double fit, fit_min;
    int n = clusters.size(), nCruces, nMutaciones, pos, c1, c2, nIterMAX, posCMejor, nFitness;
    vector< vector<int> > pActual, pSiguiente;
    vector<double> pActualFitness;
    vector<int> cromosoma, CMejor;

    const int   nFitnessMAX = 100000;
    const double nFallosMAX = 0.1*M, // FIXME: n es M?
                 probMutacion = 1/n;
    
    // Inicializar y evaluar la poblacion actual
    for (int i = 0; i < M; i++){
        // Calcular un cromosoma valido
        do{
            cromosoma.clear();

            for (int j = 0; j < n; j++)
                cromosoma.push_back(Randint(0, num_clusters-1));

            // Comprobar que ningun cluster se queda vacio
            recalcular = clusterVacio (cromosoma);
        } while(recalcular);

        pActual.push_back(cromosoma);
    }

    // Evaluar la poblacion actual
    // TO DO crear mejor una funcion para hacer esto
    for (int i = 0; i < M; i++){
        clusters = pActual[i];
        pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
    }

    nFitness = 0;
    
    for (int t = 0; nFitness < nFitnessMAX; t++){
        if (t % 10 == 0){
            for (int i = 0; i < pActual.size(); i++){
                nFitness += busquedaLocalSuave(pActual[i], pActualFitness[i], nFallosMAX);
            }
        }

        // Guardo el mejor cromosoma
        posCMejor = calcularMejorCromosoma (pActualFitness);
        CMejor = pActual[posCMejor];

        ////////////////////////////////// SELECCIÓN //////////////////////////////////

        pSiguiente.clear();

        // Seleccionar P' desde P(t-1)
        for (int i = 0; i < M; i++){
            c1 = Randint(0, M-1);

            do{
                c2 = Randint(0, M-1);
            } while (c1 == c2);
            

            pos = operadorSeleccion(c1, c2, pActualFitness);

            pSiguiente.push_back(pActual[pos]);
        }

        ////////////////////////////////// CRUZAR //////////////////////////////////

        // Recombinar P'
        nCruces = probCruce * M/2;

        // FIXME esto esta feo pero no se como arreglarlo todavia
        for (int i = 0; i < 2*nCruces; i += 2){
            pSiguiente.push_back(operadorCruceUN(pSiguiente[i], pSiguiente[i+1]));           
            pSiguiente.push_back(operadorCruceUN(pSiguiente[i+1], pSiguiente[i]));
        }
        
        for (int i = 0; i < nCruces; i++){
            pSiguiente.erase(pSiguiente.begin());
        }

        ////////////////////////////////// MUTACIÓN //////////////////////////////////

        // Mutar P'
        nMutaciones = M * n * probMutacion;

        for (int i = 0; i < nMutaciones; i++){
            pos = Randint(0, M-1);
            pSiguiente[pos] = operadorMutacionUN (pSiguiente[pos]);
        }

        ////////////////////////////////// REEMPLAZAR //////////////////////////////////
        
        // Reemplazar P(t) a partir de P(t-1) y P'
        pActual = pSiguiente;

        // Si el mejor cromosoma no esta, lo añado al final
        if (find(pActual.begin(), pActual.end(), CMejor) == pActual.end())
            pActual[M-1] = CMejor;


        ////////////////////////////////// EVALUAR //////////////////////////////////
        
        // Evaluar P(t)
        // TO DO crear mejor una funcion para hacer esto
        pActualFitness.clear();
        
        for (int i = 0; i < M; i++){
            clusters = pActual[i];
            pActualFitness.push_back(fitnessBL()); // DUDA cambiar lo de BL
            nFitness++;
        }

    }

    // Actualizo atributos
    // TO DO: mejorar esto
    posCMejor = calcularMejorCromosoma (pActualFitness);
    clusters = pActual[posCMejor];
    desvGeneral = desviacionGeneral();
    fitness = fitnessGreedy();
    infeasibility = infeasibilityGreedy();

    return clusters;

}


int PAR::calcularMejorCromosoma (const vector<double> pFitness){
    int fitnessMin = pFitness[0];

    for (int i = 0; i < pFitness.size(); i++)
        if (pFitness[i] < pFitness[fitnessMin])
            fitnessMin = i;

    return fitnessMin;
}

int PAR::calcularPeorCromosoma (const vector<double> pFitness){
    int fitnessMin = pFitness[0];

    for (int i = 0; i < pFitness.size(); i++)
        if (pFitness[i] > pFitness[fitnessMin])
            fitnessMin = i;

    return fitnessMin;
}


int PAR::operadorSeleccion (const int c1, const int c2, const vector<double> pFitness){
    if (pFitness[c1] > pFitness[c2])
        return c2;
    else 
        return c1;
}


vector<int> PAR::operadorCruceUN (const vector<int> padre1, const vector<int> padre2){
    const int n = padre1.size();
    vector<int> posP1, cromosoma(n, -1);
    int pos;

    for (int i = 0; i < n/2; i++){
        do{
            pos = Randint (0, n-1);
        } while (find(posP1.begin(), posP1.end(), pos) != posP1.end());
        
        posP1.push_back(pos);

        cromosoma[pos] = padre1[pos];
    }
    
    for (int i = 0; i < n; i++)
        if (cromosoma[i] == -1)
            cromosoma[i] = padre2[i];

    repararCromosoma(cromosoma);

    return cromosoma;
}

vector<int> PAR::operadorCruceSF (const vector<int> padre1, const vector<int> padre2){
    int r, v, pos;
    const int n = padre1.size();
    vector<int> cromosoma(n);

    r = Randint(0, n-1);
    v = Randint(0, n);

    // Hacer v iteraciones empezando en r%n
    for (int i = 0; i < v; i++){
        pos = (r + i) % n;
        cromosoma[pos] = padre1[pos];
    }

    // Hacer (n-v) iteraciones empezando en (r+v)%n
    for (int i = 0; i < (n-v);  i++){
        pos = (r + v + i) % n; // DUDA
        cromosoma[pos] = padre2[pos];
    }

    repararCromosoma(cromosoma);

    return cromosoma;
}

vector<int> PAR::operadorMutacionUN (vector<int> cromosoma){
    const int n = cromosoma.size();
    int posGen, gen, antiguoC;
    
    posGen = Randint(0, n-1);
    gen = Randint(0, num_clusters-1);
    antiguoC = cromosoma[posGen];
    cromosoma[posGen] = gen;

    if (clusterVacio(cromosoma))
        cromosoma[posGen] = antiguoC;

    return cromosoma;
}

void PAR::repararCromosoma (vector<int> & cromosoma){
    vector<int> contador(num_clusters, 0);
    int pos;
    
    for (int i = 0; i < cromosoma.size(); i++)
        contador[cromosoma[i]]++;
    
    for (int i = 0; i < contador.size(); i++)
        if (contador[i] == 0){
            do{
                pos = Randint(0, cromosoma.size() - 1);
            } while (contador[cromosoma[pos]] - 1 == 0);

            cromosoma[pos] = i;
            contador[cromosoma[pos]]--;
            contador[i]++;
        }
}