#include<bits/stdc++.h>
//#include<armadillo>
/*
 * excercise
 * read from stdio xa,xb,ta,tb,N,T
 * N: # of nodes in space.
 * T: # of nodes in time.
 */
using namespace std;
//using namespace arma;

#define db double
#define Vd vector<double>
#define VVd vector<vector<double> >
#define PI 3.14159265

void printM(VVd &m){
    for(int i=0; i<m.size(); ++i){
        for(int j=0; j<m[0].size(); j++)
            cout<<m[i][j]<<"| ";
        cout<<endl;
    }
    cout<<endl;
}

void printV(Vd &x){
    for(int i=0; i<x.size(); ++i)
        cout<<x[i]<<"|";
    cout<<endl;
}

int main(){
    db deltaX,deltaT,xa,xb,ta,tb;
    int N,T;
    cin>>xa>>xb>>ta>>tb>>N>>T;
    cout<<xa<<xb<<ta<<tb<<N<<T<<endl;
    deltaX = (xb-xa) / float(N+1);
    deltaT = (tb-ta) / float(T+1);

    VVd U(T+2,Vd (N+2,-1));
    VVd A(N,Vd (N,0));
    Vd B (N,-1);
    //mat U1(T+2,N+2,fill::zeros);
    //mat A1(N,N,fill::zeros);

    //condiciones iniciales
    for(int i=0; i<N+2; i++)
        U[0][i] = 2*sin(2*PI*i*deltaX);

    //for(int t=1; t<T+2; ++t){       //it over time
    for(int t=1; t<2; ++t){         //this is for test original line is previous line
        U[t][0] = 0;                //condicion inicial;
        U[t][N+1] = 0;
        for(int x=0; x<N; x++){
            A[x][x] = (deltaX*deltaX) + deltaT/8.0;
            if(x==0){
                A[x][x+1] = -deltaT/16.0;
                B[x] =  (deltaX*deltaX)*U[t-1][x+1];
            }
            else if(x==N-1){
                A[x][x-1] = -deltaT/16.0;
                B[x] = deltaX*U[t-1][x+1];
            }
            else{
                A[x][x-1] = -deltaT/16.0;
                A[x][x+1] = -deltaT/16.0;
                B[x] = (deltaX*deltaX)*U[t-1][x+1];
                }
        }
    }
    printM(U);
    printM(A);
    printV(B);
    return 0;
}
