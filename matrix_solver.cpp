#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <chrono>

using namespace std;



typedef vector<vector<double>> matriz;


vector<double> resolverSuperior(matriz A,vector<double> b){
    int n = A.size();
    vector<double> x (n);
    for(int i = n-1;i >= 0;i--){
        double sumaDerecha = 0;
        for(int j = i+1;j <= n-1;j++){
            sumaDerecha += x[j] * A[i][j];
        }
        x[i] = (b[i] - sumaDerecha)/(A[i][i]);
    }
    return x;
}
vector<double> resolverInferior(matriz A,vector<double> b){
    int n = A.size();
    vector<double> x (n);
    for(int i = 0;i < n;i++){
        double sumaDerecha = 0;
        for(int j = 0;j < i;j++){
            sumaDerecha += x[j] * A[i][j];
        }
        x[i] = (b[i] - sumaDerecha)/(A[i][i]);
    }
    return x;
}


void printMatriz(matriz A){
    for(auto f:A){
        for(auto x:f){
            cout << x << " ";
        }
        cout << endl;
    }
}

void printVector(vector<double> v){
    for(auto x:v) cout << x << " ";
    cout << endl;
}

void writeVector(vector<double> v, bool first){
    ofstream file;
    if(first) file.open("output.txt");
    else file.open("output.txt",fstream::app);
    for(auto x:v) file << x << endl;
    file.close();
}

void writeIsos(vector<int> isos, bool first){
    ofstream  file;
    if(first) file.open("output.txt");
    else file.open("output.txt",fstream::app);
    for(auto x:isos) file << x << endl;
    file.close();
}

void writeInitialValues(double r_i,double r_e,int m,int n,double iso,int ninst, bool first){
    ofstream  file;
    if(first) file.open("output.txt");
    else file.open("output.txt",fstream::app);
    file << r_i << " " << r_e << " " << m << " " << n << " " << iso << " " << ninst << endl;
    file.close();
}

matriz transpuesta(matriz A){
    for(int i=0; i < A.size();i++){
        for(int j=i+1;j < A.size();j++){
            swap(A[i][j],A[j][i]);
        }
    }

    return A;
}


void restarFilas(vector<double> F_i,vector<double> &F_j,double alpha){
    for(int k = 0; k < F_j.size();k++){
        F_j[k] -= alpha * F_i[k];
    }
}

void triangularSuperior(matriz &A,vector<double> &b){
    int n = A.size();
    for(int i = 0; i < n;i++){
        for(int j=i+1;j < n;j++){
            auto m_ij = A[j][i]/A[i][i];
            restarFilas(A[i],A[j],A[j][i]/A[i][i]);
            b[j] -= m_ij * b[i];
        }
    }
}


void triangularSuperiorBanda(matriz &A,vector<double> &b,int band_width){
    int n = A.size();
    for(int i = 0; i < n;i++){
        for(int j = i+1;j < min(i + band_width,n);j++){
            auto m_ij = A[j][i]/A[i][i];
            A[j][i] = 0;
            for(int k = i+1;k < min(i+band_width,n);k++){
                A[j][k] = A[j][k] - m_ij * A[i][k];
            }
            b[j] -= m_ij * b[i];
        }
    }
}

pair<matriz,matriz> calcularLUBanda(matriz A,int band_width){
    matriz L (A.size(),vector<double> (A.size()));
    int n = A.size();
    for(int i = 0; i < n;i++){
        for(int j = i+1;j < min(i + band_width,n);j++){
            auto m_ij = A[j][i]/A[i][i];
            A[j][i] = 0;
            L[j][i] = m_ij;
            for(int k = i+1;k < min(i+band_width,n);k++){
                A[j][k] = A[j][k] - m_ij * A[i][k];
            }

        }
    } 
    for(int i = 0; i < A.size();i++){
        L[i][i] = 1;
    }

    return {L,A};
}



pair<matriz,matriz> calcularLU(matriz A){
    matriz L (A.size(),vector<double> (A.size()));
    int n = A.size();
    for(int i = 0; i < n;i++){
        for(int j=i+1;j < n;j++){
            double m_ij = A[j][i]/A[i][i];
            L[j][i] = m_ij;
            restarFilas(A[i],A[j],m_ij);
        }
    } 
    for(int i = 0; i < A.size();i++){
        L[i][i] = 1;
    }

    return {L,A};
}

vector<double> resolverGauss(matriz A,vector<double> b){
    triangularSuperior(A,b);
   
    return resolverSuperior(A,b);
}

vector<double> resolverLU(matriz L,matriz U,vector<double> b){
    auto y = resolverInferior(L,b);
    auto x = resolverSuperior(U,y);
    return x;
}

vector<double> matrizVectorMulti(matriz A,vector<double> x){
    vector<double> y (x.size());
    for(int i = 0; i < A.size();i++){
        double suma = 0;
        for(int j = 0; j < A.size();j++){
            suma += A[i][j] * x[j];
        }
        y[i] = suma;
    }
    return y;
}

matriz matrizMatrizMulti(matriz A,matriz B){
    int n = A.size();
    matriz AB (A.size(),vector<double> (A.size()));
    for(int i = 0; i<n; i++){
        for(int j = 0; j < n;j++){
            double suma = 0;
            for(int k = 0; k < n;k++){
                suma += A[i][k] * B[k][j];
            }
            AB[i][j] = suma;
        }
    }
    return AB;
}

matriz crearMatrizSistema(vector<double> r, int n, double delta_r,double delta_theta){
    int m = r.size() - 1;
    matriz A ((m-1)*n, vector<double> ((m-1) * n));
    for(int i = 1; i < m;i++){
        for(int j = 0; j < n;j++){
            
            int p_ij = (i-1)*n + j;
            int p_ij1 = (i-1)*n + (j + 1)%n;
            int p_ijm1 = (i-1)*n + (j-1+n)%n;
            int p_im1j = (i-2)*n + j;
            int p_i1j = (i)*n + j;

            
            
            double c_ij = 1/(r[i] * delta_r) - 2/(r[i]*r[i]*delta_theta*delta_theta) - 2/(delta_r*delta_r);
            double c_ij1 = 1/(r[i]*r[i]*delta_theta*delta_theta);
            double c_ijm1 = c_ij1;
            double c_im1j = 1/(delta_r * delta_r) - 1/(r[i]*delta_r);
            double c_i1j = 1/(delta_r*delta_r);
            int fila = (i-1)*n + j;
           
            A[fila][p_ij] = c_ij;
            A[fila][p_ij1] = c_ij1;
            A[fila][p_ijm1] = c_ijm1;
            if(i != 1){
                A[fila][p_im1j] = c_im1j;
            }
            if(i != m-1){
                A[fila][p_i1j] = c_i1j;
            }
            
        }
    }
    
    return A;
}

vector<double> crearBSistema(double r_i,int m,vector<double> temperatura_interna, vector<double> temperatura_externa,double delta_r){
    int n = temperatura_interna.size();
    vector<double> b ((m-1)*n);
    for(int j = 0; j < n; j++){
        double c_im1j = 1/(delta_r * delta_r) - 1/(r_i*delta_r);
        double c_i1j = 1/(delta_r*delta_r);
        b[j] = -c_im1j * temperatura_interna[j];
        b[(m-2)*n + j] = -c_i1j * temperatura_externa[j];
    }

    return b;

}

pair<vector<double>,double> crearRs(double r_i,double r_e, int particiones){
    vector<double> r (particiones + 1);
    r[0] = r_i;
    double width = (r_e - r_i) / particiones;
    for(int k = 1; k < particiones;k++){
        r[k] = r[k-1] + width;
    }
    r[particiones] = r_e;
    return {r,width};
}
/*
vector<double> resolverSistema(double r_i,double r_e,int particion_r,int particion_theta,vector<double> temperaturas_internas,vector<double> temperaturas_externas){
    auto R = crearRs(r_i,r_e,particion_r);
    auto r = R.first;
    double delta_r = R.second;
    double delta_theta = 2 * M_PI / particion_theta;
    matriz A = crearMatrizSistema(r,temperaturas_internas.size(),delta_r,delta_theta);
    vector<double> b = crearBSistema(r_i,r.size() - 1,temperaturas_internas,temperaturas_externas,delta_r);
    triangularSuperiorBanda(A,b,particion_theta+1);
    vector<double> x = resolverSuperior(A,b);
    cout << r_i << " " << r_e << " " << particion_r << " " << particion_theta << endl;

    return x;
}*/

vector<int> calcularIso(vector<double> x,double iso,int n,int m){
    vector<int> isos (n);
    int min;
    for(int j = 0; j < n;j++){
        min = 1;
        for(int i = 1; i < m;i++){
            if(abs(x[(min-1)*n + j] - iso) > abs(x[(i-1)*n + j]- iso)) min = i;
        }
        isos[j] = min;
    }
    return isos;
}

void resolverSistemaGauss(double r_i,double r_e,int particion_r,int particion_theta,vector<vector<double>> temperaturas_internas,vector<vector<double>> temperaturas_externas,double iso){
    int ninst = temperaturas_internas.size();
    auto R = crearRs(r_i,r_e,particion_r);
    auto r = R.first;
    double delta_r = R.second;
    double delta_theta = 2 * M_PI / particion_theta;
    matriz A = crearMatrizSistema(r,temperaturas_internas[0].size(),delta_r,delta_theta);
    matriz Apri = A;
    for(int i = 0; i < ninst;i++){

        vector<double> b = crearBSistema(r_i,r.size()-1,temperaturas_internas[i],temperaturas_externas[i],delta_r);
        triangularSuperiorBanda(Apri,b,particion_theta+1);
        vector<double> x = resolverSuperior(Apri,b);
        vector<int> isos = calcularIso(x,iso,particion_theta,particion_r);
        writeVector(x,false);
        writeIsos(isos,false);
        Apri = A;


    }
}

void resolverSistemaLU(double r_i,double r_e,int particion_r,int particion_theta,vector<vector<double>> temperaturas_internas,vector<vector<double>> temperaturas_externas,double iso){
    int ninst = temperaturas_internas.size();
    auto R = crearRs(r_i,r_e,particion_r);
    auto r = R.first;
    double delta_r = R.second;
    double delta_theta = 2 * M_PI / particion_theta;
    matriz A = crearMatrizSistema(r,temperaturas_internas[0].size(),delta_r,delta_theta);
    auto p = calcularLUBanda(A,particion_theta+1);
    matriz L = p.first;
    matriz U = p.second;
    for(int i = 0; i < ninst;i++){

        vector<double> b = crearBSistema(r_i,r.size()-1,temperaturas_internas[i],temperaturas_externas[i],delta_r);
        //vector<double> y = resolverInferior(L,b);
        vector<double> x = resolverSuperior(U,resolverInferior(L,b));
        vector<int> isos = calcularIso(x,iso,particion_theta,particion_r);
        writeVector(x,false);
        writeIsos(isos,false);
    }
}

void resolverSistemaGaussSinBanda(double r_i,double r_e,int particion_r,int particion_theta,vector<vector<double>> temperaturas_internas,vector<vector<double>> temperaturas_externas,double iso){
    int ninst = temperaturas_internas.size();
    auto R = crearRs(r_i,r_e,particion_r);
    auto r = R.first;
    double delta_r = R.second;
    double delta_theta = 2 * M_PI / particion_theta;
    matriz A = crearMatrizSistema(r,temperaturas_internas[0].size(),delta_r,delta_theta);
    matriz Apri = A;
    for(int i = 0; i < ninst;i++){

        vector<double> b = crearBSistema(r_i,   r.size()-1,temperaturas_internas[i],temperaturas_externas[i],delta_r);
        triangularSuperior(Apri,b);
        vector<double> x = resolverSuperior(Apri,b);
        vector<int> isos = calcularIso(x,iso,particion_theta,particion_r);
        writeVector(x,false);
        writeIsos(isos,false);
        Apri = A;


    }
}

void resolverSistemaLUSinBanda(double r_i,double r_e,int particion_r,int particion_theta,vector<vector<double>> temperaturas_internas,vector<vector<double>> temperaturas_externas,double iso){
    int ninst = temperaturas_internas.size();
    auto R = crearRs(r_i,r_e,particion_r);
    auto r = R.first;
    double delta_r = R.second;
    double delta_theta = 2 * M_PI / particion_theta;
    matriz A = crearMatrizSistema(r,temperaturas_internas[0].size(),delta_r,delta_theta);
    auto p = calcularLU(A);
    matriz L = p.first;
    matriz U = p.second;
    for(int i = 0; i < ninst;i++){

        vector<double> b = crearBSistema(r_i,r.size()-1,temperaturas_internas[i],temperaturas_externas[i],delta_r);
        //vector<double> y = resolverInferior(L,b);
        vector<double> x = resolverSuperior(U,resolverInferior(L,b));
        vector<int> isos = calcularIso(x,iso,particion_theta,particion_r);
        writeVector(x,false);
        writeIsos(isos,false);
    }
}




void ingresoDatos(string test_file,int tipo,bool banda){
    fstream test;
    test.open("tests\\" + test_file + ".txt",ios::in);
    double r_i,r_e,iso;
    int n,m,ninst;
    test >> r_i >> r_e >> m >> n >> iso >> ninst;
    writeInitialValues(r_i,r_e,m-1,n,iso,ninst,true);
    vector<vector<double>> temperaturas_internas (ninst,vector<double>(n));
    vector<vector<double>> temperaturas_externas (ninst,vector<double>(n));
    for(int j = 0; j < ninst;j++){
        for(int i = 0; i < n;i++){
            test >> temperaturas_internas[j][i];
        }
        for(int i = 0; i < n;i++){
            test >> temperaturas_externas[j][i];
        }
    }
    
    test.close();

    if(tipo == 0){
        if(banda){
            resolverSistemaGauss(r_i,r_e,m-1,n,temperaturas_internas,temperaturas_externas,iso);
        }else{
            resolverSistemaGaussSinBanda(r_i,r_e,m-1,n,temperaturas_internas,temperaturas_externas,iso);
        }
    }else{
        if(banda){
            resolverSistemaLU(r_i,r_e,m-1,n,temperaturas_internas,temperaturas_externas,iso);
        }else{
            resolverSistemaLUSinBanda(r_i,r_e,m-1,n,temperaturas_internas,temperaturas_externas,iso);
        }
    }

}


int main(){
    /*fstream time_file;
    time_file.open("time_tests.txt",fstream::app);*/
    string test_file;
    cin >> test_file;
    ingresoDatos(test_file,1,true);
    /*
    chrono::steady_clock::time_point LU_begin = chrono::steady_clock::now();
    for(int i = 0; i < 1;i++){
        ingresoDatos(test_file,1,true);
    }
    chrono::steady_clock::time_point LU_end = chrono::steady_clock::now();
    cout << "Tiempo transcurrido para LU: " << chrono::duration_cast<chrono::milliseconds>(LU_end - LU_begin).count() << "[ms]" << endl;
    chrono::steady_clock::time_point gauss_begin = chrono::steady_clock::now();
    for(int i = 0; i < 1;i++){
        ingresoDatos(test_file,0,true);
    }
    chrono::steady_clock::time_point gauss_end = chrono::steady_clock::now();
    cout << "Tiempo transcurrido para Gauss: " << chrono::duration_cast<chrono::milliseconds>(gauss_end - gauss_begin).count() << "[ms]" << endl;*/
    /*bool gauss_lenta_N = false,gauss_lenta_M = false,LU_lenta_N = false,LU_lenta_M = false;
    for(int i = 20; i >= 1;i--){
        gauss_lenta_N = false; LU_lenta_N = false;
        for(int j = 1;j <= 20;j++){
            int m = 5*i;
            int n = 5*j;
            string file = "size_test_" + to_string(m) + "_" + to_string(n);
            time_file << "M = " << m << " , N = " << n << " :" << endl;
            if(!gauss_lenta_N && !gauss_lenta_M){
                chrono::steady_clock::time_point gauss_begin = chrono::steady_clock::now();
                ingresoDatos(file,0,false);
                chrono::steady_clock::time_point gauss_end = chrono::steady_clock::now();
                int gauss_time = chrono::duration_cast<chrono::milliseconds>(gauss_end - gauss_begin).count();
                if(gauss_time > 120000){
                    if(j == 1){
                        gauss_lenta_M = true;
                    }else{
                        gauss_lenta_N = true;
                    }
                    time_file << "Gauss: >120000 [ms]" << endl;
                }else{ 
                    time_file << "Gauss: " << gauss_time << "[ms]" << endl;
                }
            }else{
                time_file << "Gauss: >120000 [ms]" << endl;       
            }
            if(!LU_lenta_M && !LU_lenta_N){
                chrono::steady_clock::time_point LU_begin = chrono::steady_clock::now();
                ingresoDatos(file,1,false);
                chrono::steady_clock::time_point LU_end = chrono::steady_clock::now();
                int lu_time = chrono::duration_cast<chrono::milliseconds>(LU_end - LU_begin).count();
                if(lu_time > 120000){
                    if(j == 1){
                        LU_lenta_M = true;
                    }else{
                        LU_lenta_N = true;
                    }
                    time_file << "LU: >120000 [ms]" << endl;
                }else{ 
                    time_file << "LU: " << lu_time << "[ms]" << endl;
                }
            }else{
                time_file << "LU: >120000 [ms]" << endl;
            }
            chrono::steady_clock::time_point gaussB_begin = chrono::steady_clock::now();
            ingresoDatos(file,0,true);
            chrono::steady_clock::time_point gaussB_end = chrono::steady_clock::now();
            time_file << "Gauss (Banda): " << chrono::duration_cast<chrono::milliseconds>(gaussB_end - gaussB_begin).count() << "[ms]" << endl;
            chrono::steady_clock::time_point LUB_begin = chrono::steady_clock::now();
            ingresoDatos(file,1,true);
            chrono::steady_clock::time_point LUB_end = chrono::steady_clock::now();
            time_file << "LU (Banda): " << chrono::duration_cast<chrono::milliseconds>(LUB_end - LUB_begin).count() << "[ms]" << endl;
            
        }
    }
    chrono::steady_clock::time_point LUB_begin = chrono::steady_clock::now();

    ingresoDatos("size_test_100_100",0,true);
    chrono::steady_clock::time_point LUB_end = chrono::steady_clock::now();
    cout << "LU (Banda): " << chrono::duration_cast<chrono::milliseconds>(LUB_end - LUB_begin).count() << "[ms]" << endl;*/

    
}