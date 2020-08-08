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
            delete [][] vet; //desaloca a memoria alocada no construtor
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
    float X = 10.0;                                        // Largura do percurso
    float a = 3.32;                                        // Comprimento da Onda
    float T = 1.8;                                         // Tempo final do percurso
    float c = 1.0;                                         // Celeridade da onda
    int Nx = 1500+1;                                       // Numero de iteracoes no espaço
    int Nt = 2720+1;                                       // Numero de iteracoes no tempo

    double dx = 1.0/150;                                    // Passo no espaço
    double dt = 1.8/272;                                    // Passo no tempo

    double C = c*dt/dx;                                     // Numero de Courant
    double C2 = C*C;

    Matriz_zerada u(Nt, Nx);

    for (int i = 0; i < Nx - 1; i++){
        //cout << i*dx << endl;
        if (a <= i*dx && i*dx <= 2*a){
            u.set(0, i, 0.5*(1 - (cos(2*M_PI*i*dx/a))));                  //definindo os valores iniciais no espaco para t = 0
        }
    }
    
    for (int i = 1; i < Nx - 1; i++){
        u.set(1, i, (u.get(0, i) + C2*( u.get(0, i - 1) - 2 * u.get(0, i) + u.get(0, i + 1) ))); //t = 1
    }
 
    for (int n = 1; n < Nt - 2; n++){
        for (int i = 1; i < Nx - 1; i++){
            u.set(n + 1, i, (-u.get(n - 1, i) + 2 * u.get(n, i) + C2*(u.get(n, i - 1) - 2 * u.get(n, i) + u.get(n, i + 1))) );
        }
    }
    
    ofstream myfile;
    myfile.open("data.dat");

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