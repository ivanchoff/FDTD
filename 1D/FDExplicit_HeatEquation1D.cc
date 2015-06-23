#include<bits/stdc++.h>
#include<armadillo>
/*
 * compile : g++ FDHeatEquation1D.cc -o t -O2 -larmadillo
 *
 * excercise
 * read from stdio xa,xb,ta,tb,N,T
 * N: # of nodes in space.
 * T: # of nodes in time.
 * Nota: si el factor de difusividad (gamma) es una constante la matriz A
 *       solo se debe hallar una vez.
 */

using namespace std;
using namespace arma;
#define db double
#define PI 3.14159265

//initial condition, it depends of each problem
db f(db i, db deltax){
    db x = i*deltax;
    return 2*sin(2*PI*x); //initial condition
}

//print temperature for each node
void print(mat &u, db deltax){
    int size =  u.n_cols;
    vec x =  vec(size);
    mat p = mat(size,2);

    for(int i=0;i<size;i++)
        x(i)=i*deltax;

    for(int i=0; i<u.n_rows; ++i){
        p.submat(0,0,size-1,0) = x;
        p.submat(0,1,size-1,1) = u.row(i).t();
        p.print();
    }
}

int main(){
    db deltaX,deltaT,xa,xb,ta,tb,gamma,r;  //gamma, el valor que multiplica la segunda derivada
    int N,T;
    cin>>xa>>xb>>ta>>tb>>N>>T>>gamma;         // initial condition

    deltaX = (xb-xa) / float(N+1);
    deltaT = (tb-ta) / float(T+1);
    r = (gamma*deltaT)/(deltaX*deltaX);

    if(r<=0  || r>=0.5 || deltaT>(deltaX*deltaX)/2.0){
        cout<<"r:"<<r<<" k"<<deltaT<<endl;
        cout<<"error"<<endl;
        return 0;
    }
    mat U = mat(T+2,N+2);   //matriz solution

    //condiciones iniciales
    for(int i=0; i<N+2; i++){
        U(0,i) = f(i,deltaX);
    }
    for(int t=1; t<T+2; ++t){       //all nodes time
        U(t,0) = 0;
        U(t,N+1) = 0;              //condicion inicial;
        for(int pos=1; pos<=N; ++pos){
            U(t,pos) = r*U(t-1,pos-1) + U(t-1,pos)*(1-(2*r)) + r*U(t-1,pos+1);
        }
    }
    print(U,deltaX);
    return 0;
}
