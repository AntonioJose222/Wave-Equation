#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <tgmath.h>
#define _USE_MATH_DEFINES
# define M_PI 3.14159265358979323846 
using namespace std;

class Matriz_zerada{

    public:
        Matriz_zerada(int l, int c){
            nl = l;
            nc = c;
            vet = new double[nl*nc]; //alocando o vetor de vetores

            for (int i = 0; i < nl*nc; i++){
                vet[i] = 0.0;
            }

        }
        ~Matriz_zerada() {
            delete [] vet; //desaloca a memoria alocada no construtor
        }

        double get(int i, int j){

            int k = detInd(i, j);

            if (k != -1){
                return vet[k];
            } else {
                //cout << "indice invalido!" << endl;
                return 0;
            }

        }

        void set(int i, int j, double val){

            int k = detInd(i, j);

            if (k != -1){
                vet[k] = val;
            } else {
                //cout << "indice invalido!" << endl;
            }

        }

    private:

        int nl;      //numero de linhas
        int nc;      //numero de colunas
        double *vet; //vetor de tamanho nl*nc

        int detInd(int i, int j){

            if (i >= 0 && i < nl && j >= 0 && j < nc){
                return i*nc + j;
            } else {
                return -1;
            }

        }
};

double vel(int x, float a, float c){

    if (0 <= x && x <= a){
        return -c*(3.1415*sin(2*3.1415*x/a))/a;       //-c*f'(x)
    } else {
        return 0.0;
    }

}

int main()
{
    cout.precision(17);
    float X = 10.0;                                     // Largura do percurso
    float T = 1;                                        // Tempo final do percurso
    float c = 2.0;                                      // Celeridade da onda
    double dx = 0.01;                                    // Passo no espaço
    double dt = 0.0005;                                 // Passo no tempo
    int Nx = int(X/dx);                                 // Numero de iteracoes no espaço
    int Nt = int(T/dt);                                 // Numero de iteracoes no tempo

    double C = c*dt/dx;                                 // Numero de Courant

    Matriz_zerada u(Nt, Nx);

    for (int i = 0; i < Nx - 1; i++){
        if (0 <= i*dx && i*dx <= 1){
            u.set(0, i, 2*(i*dx)*(1 - i*dx));                  //definindo os valores iniciais no espaco para t = 0
        }
    }
 
    for (int n = 1; n < Nt - 2; n++){
        for (int i = 1; i < Nx - 1; i++){
            u.set(n, i, 0.5 * ((u.get(n - 1, i + 1) + u.get(n - 1, i - 1)) - C*(u.get(n - 1, i + 1) - u.get(n - 1, i - 1))));
        }
    }
    
    ofstream myfile;
    myfile.open("data-lax.dat");

    //myfile << Nx << " " << Nt << endl << endl;
    
    for (int i = 0; i < Nt; i += 50){
        for (int j = 0; j < Nx; j++){

            myfile << std::fixed << std::setprecision(20) << j*dx << " " << u.get(i, j) << "\n";

        }
        myfile << "\n\n";
    }
    myfile.close();

    return 0;
}