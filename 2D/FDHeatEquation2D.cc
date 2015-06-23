#include<bits/stdc++.h>
#include<armadillo>

/*compile : g++ FDHeatEquation2D.cc -o t -O2 -larmadillo
 * Note: for A matrix, Nx and Ny must be mayor than five
 */

using namespace std;
using namespace arma;

#define db double
#define PI 3.14159265

//if gamma is constant A matrix is loaded just one time.
void load_Matrix_A(mat &A, int Nx, int Ny,db sx, db sy){
    db tmp = 1 + (2*sy) + (2*sx);
    int node,upper,lower;
    for(int x=1; x<Ny-1; ++x){        // iterate over rows
        for(int y=1; y<Nx-1; ++y){    //iterate over cols
            node = ( (x-1)*(Nx-2) ) + (y-1);
            upper = node + (Nx-2);
            lower = node - (Nx-2);
            A(node, node) =  tmp;
            if(x==1){
                A(node, upper) = -sy;
                if(y==1){                  //bottom-left
                    A(node, node+1) = -sx;
                }
                else if(y==Nx-2){         //bottom-right
                    A(node, node-1) =  -sx;
                }else{                    //bottom-middle
                    A(node, node-1) = -sx;
                    A(node, node+1) = -sx;
                }
            }
            else  if(x==Ny-2){
                A(node, lower) =  -sy;
                if(y==1){                 //top-left
                    A(node, node+1) = -sx;
                }
                else if(y==Nx-2){        //top-right
                    A(node, node-1) = -sx;
                }else{                   //top-middle
                    A(node, node-1) = -sx;
                    A(node, node+1) = -sx;
                }
            }
            else{
                A(node, lower) = -sy;
                A(node, upper) = -sy;
                if(y==1){                //left
                    A(node, node+1) = -sx;
                }
                else if(y==Nx-2){       //right
                    A(node, node-1) = -sx;
                }
                else{    //complete equation central points
                    A(node, node-1) = -sx;
                    A(node, node+1) = -sx;
                }
            }
        }
    }
}

// Load Matriz U with boundary and initial conditions
void load_Matrix_U(mat &U, int Nx, int Ny, int T){
    int nodes = Nx*Ny;
    U(0,(Ny/2)*Nx+(Nx/2))=100; //initial condition
    for(int t=0; t<T; ++t ){
        for(int i=0; i<=Nx; ++i){ //horizo
            U(t,i) = 0;
            U(t,(nodes-1)-i) = 0;
        }
        for(int k=(2*Nx)-1; k<nodes-Nx; k+=Nx){ //vertical
            U(t,k) = 0;
            U(t,k+1) = 0;
        }
    }

}
//insert x solution of ax=b in matrix solution U
void insert_X_in_U(mat &U, vec &X, int Nx, int Ny, int t){
    int posX=0,posU=0;
    int row = Nx-2;
    for(int i=1; i<Ny-1; ++i){
        posU = (i*Nx)+1;
        posX = row*(i-1);
        U.submat(t,posU, t,posU+row-1) =  X.subvec(posX,posX+row-1).t();
    }
}

//print heat at time 0 and T
void printU(mat &U, int Nx, int Ny, db deltaX, db deltaY, int T){
    for(int t=0; t<T; t+=T-1){
        for(db y=0; y<Ny; ++y)
            for(db x=0; x<Nx; ++x)
                cout<<x*deltaX<<" "<<y*deltaY<<" "<<U(t,y*Nx+x)<<endl;
        cout<<endl;
    }
}
int main(){
    db deltaX,deltaY,deltaT,xa,xb,ya,yb,ta,tb,k;
    int Nx,Ny,T,nodes,nodesA;
    cin>>xa>>xb>>ya>>yb>>ta>>tb>>Nx>>Ny>>T>>k;

    deltaX = (xb-xa) / float(Nx-1);
    deltaY = (yb-ya) / float(Ny-1);
    deltaT = (tb-ta) / float(T+1);
    nodes = Nx*Ny;
    nodesA = (Nx-2)*(Ny-2);

    mat U = mat(T,nodes);      //matriz solution
    mat A = mat(nodesA,nodesA);  //A
    vec B = vec(nodesA);        //B
    vec X = vec(nodesA);        // X from AX=B;

    db sx = (k*deltaT)/(deltaX*deltaX);
    db sy = (k*deltaT)/(deltaY*deltaY);

    //load matriz A (only if k is constant!)
    load_Matrix_A(A, Nx, Ny, sx, sy);
    load_Matrix_U(U, Nx, Ny, T);

    for(int t=0; t<T-1; ++t){ //each time
        int pos = 0;
        for(int i=Nx+2; i<(nodes-Nx); i++)
            if(i%Nx!=0 && (i-1)%Nx!=0)
                B(pos++)=U(t,i-1);
        X = solve(A,B);
        insert_X_in_U(U,X, Nx, Ny, t+1);
    }
    //U.print();
    printU(U, Nx, Ny, deltaX, deltaY, T);
    return 0;
}
