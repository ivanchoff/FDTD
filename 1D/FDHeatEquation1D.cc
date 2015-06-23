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

//when gamma is constant matriz A can be calculate only one time.
void load_matriz_A(mat &A, int N, db deltaX2, db deltaT, db gamma){
    for(int x=0; x<N; x++){
        A(x,x) = deltaX2+ 2*deltaT*gamma; //Corregido
        if(x==0){ //time 0
            A(x,x+1) = -deltaT*gamma;
        }
        else if(x==N-1){ //time N
            A(x,x-1) = -deltaT*gamma;
        }
        else{
            A(x,x-1) = -deltaT*gamma;  //Corregido
            A(x,x+1) = -deltaT*gamma;  //Corregido
        }
    }

}

int main(){
    db deltaX,deltaT,xa,xb,ta,tb,gamma,deltaX2;  //gamma, el valor que multiplica la segunda derivada
    int N,T;
    cin>>xa>>xb>>ta>>tb>>N>>T>>gamma;         // initial condition

    deltaX = (xb-xa) / float(N+1);
    deltaT = (tb-ta) / float(T+1);
    deltaX2=deltaX*deltaX;

    mat U = mat(T+2,N+2);   //matriz solution
    mat A = mat(N,N);       //A
    vec B = vec(N,1);      //B
    vec X = vec(N);          // X from AX=B;

    //condiciones iniciales
    for(int i=0; i<N+2; i++){
        U(0,i) = f(i,deltaX);
    }
    load_matriz_A(A, N, deltaX2, deltaT, gamma);
    for(int t=1; t<T+2; ++t){       //all nodes time
        U(t,0) = 0;
        U(t,N+1) = 0;              //condicion inicial;
        for(int x=0; x<N; x++){
            B(x) = deltaX2 * U(t-1,x+1); //x+1 offset in U
        }
        X = solve(A,B); // solver for AX=B
        U.submat(t,1,t,N)= X.t(); //inser X in matrix solution
    }
    print(U,deltaX);
    return 0;
}
