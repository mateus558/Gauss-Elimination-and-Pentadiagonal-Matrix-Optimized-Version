#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <limits>

using namespace std;

const double a = 10.0;

/* run this program using the console pauser or add your own getch, system("pause") or input loop */
template <typename T>
void print_matrix(T **A, int N);
template <typename T>
void print_array(T *A, int N);
template <typename T>
int read_matrix(char *fname, T ***matrix, int *N);
template <typename T>
void free_matrix(T **A, int N);
template <typename T>
void free_array(T* A);
double** assemble_system(double** A, double *b, int N, double dA, double dB, double Dp, double dC, double dD);
double** swap_lines(double **A, int N, int L1, int L2);
double LU_decomposition(double** A, double ***L, double ***U, int N);
int is_singular(double **A, int N);
double* back_substitution(double** A, int N, size_t *ops);
double* back_substitution_penta(double** A, int N, size_t *ops);
double* gauss_elimination(double **A, int N, size_t *ops, int verbose);
double* gauss_elimination_penta(double **A, int N, size_t *ops, int verbose);

void waitUserAction(void){
    cout << "\nPressione ENTER para continuar..." << endl;
    cin.ignore(numeric_limits<streamsize>::max(),'\n');
    cin.get();
}

int main(int argc, char** argv) {
	int i, j, N, giveB = 0, printSystem = 0, printSol = 0, optmized = 0;
	double **A;
	double *x = NULL, *b = NULL;
	double dA, dB, Dp, dC, dD;
	size_t ops = 0;
	
	if(argc == 1){
		cout << "Entre a dimensao do sistema: ";	
		cin >> N;
		cout << "Entre o coeficiente dA: ";
		cin >> dA;
		cout << "Entre o coeficiente dB: ";
		cin >> dB;
		cout << "Entre o coeficiente Dp: ";
		cin >> Dp;
		cout << "Entre o coeficiente dC: ";
		cin >> dC;
		cout << "Entre o coeficiente dD: ";
		cin >> dD;
		cout << "Fornecer vetor b? [1|0]: ";
		cin >> giveB;

		if(giveB){
			b = new double[N];

			for(i = 0; i < N; i++){
				cout << "b[" << i << "]: "; 
				cin >> b[i];
			}
		}

		A = assemble_system(A, b, N, dA, dB, Dp, dC, dD);
	}else{
		read_matrix<double>(argv[1], &A, &N);
	}
	
	cout << "Mostrar matriz de coeficientes? [0|1]: ";
	cin >> printSystem;

	if(printSystem)	print_matrix<double>(A, N);	
	
	cout << "Qual versao da eliminacao de gauss? [otimizado = 1 | nao otimi = 0]: ";
	cin >> optmized;
	
	if(optmized){
		x = gauss_elimination_penta(A, N, &ops, 0);
	}else{
		x = gauss_elimination(A, N, &ops, 0);
	}

	cout << "Mostrar solucao? [0|1]: ";
	cin >> printSol;
	if(printSol) print_array<double>(x, N);

	cout << ops << " operacoes." << endl;

	free_matrix<double>(A, N);
	free_array<double>(x);
	free_array<double>(b);
	
	waitUserAction();
	
	return 0;
}

template <typename T>
int read_matrix(char *fname, T*** matrix, int *N){
	int i, j, lines, cols, size;
	double val;
	FILE *file = NULL;
	
	file = fopen(fname, "r");
	if(!file){
		cout << "The file could not be opened!\n";
		return -1;
	}

	fscanf(file, "%d", N);
	lines = (*N);
	cols = lines+1;
	(*matrix) = new T*[lines];	
	for(i = 0; i < lines; ++i){ (*matrix)[i] = new T[cols]; }
	
	for(i = 0; i < lines; ++i){
		for(j = 0; j < cols; ++j){
			fscanf(file, "%lf", &(*matrix)[i][j]);
		}	
	}
	
	fclose(file);	
	return cols;
}

template <typename T>
void free_matrix(T **A, int N){
	int i;

	for(i = 0; i < N; i++){
		delete[] A[i];
	}
	delete[] A;
}

template <typename T>
void free_array(T* A){
	delete[] A;
}

double** assemble_system(double **A, double *b, int N, double dA, double dB, double Dp, double dC, double dD){
	int i, j;

	A = new double*[N];
	for(i = 0; i < N; i++){
		A[i] = new double[N+1];
		for(j = 0; j < N; j++){
			if ((i-j) == 2){
				A[i][j] = dA;
			}else if ((i-j) == 1){
				A[i][j] = dB;
			}else if ((i-j) == 0){
				A[i][j] = Dp;
			}else if ((i-j) == -2){
				A[i][j] = dC;
			}else if ((i-j) == -1){
				A[i][j] = dD;
			}else{
				A[i][j] = 0;
			}
		}
	}

	if(b == NULL){
		srand (time(NULL));

		for(i = 0; i < N; i++){
			A[i][N] = ((double)rand()/(double)(RAND_MAX)) * a;
		}
	}else{
		for(i = 0; i < N; i++){
			A[i][N] = b[i];
		}
	}

	return A;
}

template <typename T>
void print_array(T *A, int N){
	int i;

	for(i = 0; i < N-1; i++){
		cout <<"X[" << i << "]: " << A[i] << endl;
	}
	cout <<"X[" << i << "]: " << A[i] << endl;
}

template <typename T>
void print_matrix(T **A, int N){
	int i, j;
	
	for(i = 0; i < N; i++){
    	for(j = 0; j <= N; j++){
    		if(j == N) cout << " | ";
    		cout << A[i][j] << " ";
    	}
    	cout << endl;	
   }
}

double LU_decomposition(double** A, double ***L, double ***U, int N){
	int i, j, k;
	double **Lk = NULL, **Uk = NULL;
	double sumu = 0.0, suml = 0.0, det = 1.0;
	
	Lk = new double*[N];	
	Uk = new double*[N];	
	
	for(i = 0; i < N; ++i){ Uk[i] = new double[N]; }
	for(i = 0; i < N; ++i){ Lk[i] = new double[N]; }
	for(i = 0; i < N; ++i){ Lk[i][i] = 1; }
	
	for(k = 0; k < N; ++k){
		Uk[k][k] = A[k][k];
		
		for(i = k+1; i < N; ++i){
			Lk[i][k] = A[i][k]/A[k][k];
			Uk[k][i] = A[k][i];
		}
		
		for(i = k+1; i < N; ++i){
			for(j = k+1; j < N; ++j){
				A[i][j] = A[i][j] - Lk[i][k]*Uk[k][j];
			}
		}
	}
	
	*U = Uk;
	*L = Lk;
	
	for(i = 0; i < N; ++i){
		det *= Uk[i][i];
	}

	return det;
}

int is_singular(double **A, int N){
	int i, j;
	double **L = NULL, **U = NULL, **X = NULL;
	
	if(A[0][0] == 0) return -1;
	
	for(i = 1; i <= N; ++i){
		if(LU_decomposition(X, &L, &U, i) == 0){
			return 1;	
		}
		for(j = 0; j < i; ++j){
			delete[] L[j];
			delete[] U[j];
		}
		delete[] L;
		delete[] U;
	}
	
	return 0;
}

double** swap_lines(double **A, int N, int L1, int L2){
	double *temp = NULL;

	temp = A[L1];
	A[L1] = A[L2];
	A[L2] = temp;
	
	return A;
}

double** conditioned_matrix(double **A, int N){
	int i, j, l = -1;
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < i; ++j){
			if(fabs(A[j][i]) > fabs(A[i][i]) && A[j][i] != 0){ l = j; break; }
		}
		
		if(l != -1){
			for(j = 0; j < i; ++j){
				if(fabs(A[j][i]) > fabs(A[l][i])) l = j;
			}
		}
		
		for(j = i+1; j < N; ++j){
			if(fabs(A[j][i]) > fabs(A[i][i]) && A[j][i] != 0){ l = j; break; }
		}
		
		if(l != -1){
			for(j = i+1; j < N; ++j){
				if(fabs(A[j][i]) > fabs(A[l][i])) l = j;
			}
			A = swap_lines(A, N, i, l);
		}
		l = -1;
	}
	
	return A;
}

double* back_substitution_penta(double** A, int N, size_t *ops){
	int i, j;
	double sum = 0.0, *x = NULL;
	
	x = new double[N];
	
	x[N-1] = A[N-1][N]/A[N-1][N-1];
	(*ops) += 1;

	for(i = N-1; i >= 0; --i){
		for(j = i+1, sum = A[i][N]; j < N; ++j){
			if(abs(i-j) > 2) continue;
			sum -= A[i][j]*x[j];
			(*ops) += 1;
		}
		x[i] = sum/A[i][i];
		(*ops) += 1;
	}
	
	return x;
}


double* back_substitution(double** A, int N, size_t *ops){
	int i, j;
	double sum = 0.0, *x = NULL;
	
	x = new double[N];
	
	x[N-1] = A[N-1][N]/A[N-1][N-1];
	(*ops) += 1;
	
	for(i = N-1; i >= 0; --i){
		for(j = i+1, sum = A[i][N]; j < N; ++j){
			sum -= A[i][j]*x[j];
			(*ops) += 1;
		}
		x[i] = sum/A[i][i];
		(*ops) += 1;
	}
	
	return x;
}

double* gauss_elimination(double **A, int N, size_t *ops, int verbose){
	int i, j, k, l = -1, flag = 1;
	double ratio, sum = 0.0, *x = NULL;
	
	if(verbose){	
		cout << "\nPivoted Matrix: \n";
		print_matrix<double>(A, N);
	}

	/*if(is_singular(A, N)){
		cerr << "\nSingular matrix.\n";
		return NULL;
	}*/

	x = new double[N];

	for(k = 0; k < N - 1; ++k){
		for(i = k + 1; i < N; ++i){
			if(A[k][k] == 0){
				for(j = 0; j < k; ++j){
					if(l != -1 && A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				for(j = k+1; j < N; ++j){
					if(l != -1 && A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				if(l != -1){
					A = swap_lines(A, N, k, l);
					(*ops) += 3;
				}
				l = -1;
			}
			ratio = A[i][k]/A[k][k];
			A[i][k] -= A[i][k];
			for(j = k + 1; j < N; ++j){
				A[i][j] -= ratio * A[k][j];
				(*ops) += 1;
			}
			A[i][N] -= ratio * A[k][N];
			(*ops) += 3;
		}
	}
	
	if(verbose){
		cout << "\nScalonated matrix:\n";
		print_matrix<double>(A, N);
	}

	for(i = 0; i < N; ++i){
		if(A[N-1][i] != 0){ flag = 0; break; }
	}	
	
	if(flag){
		cout << "\nInfinity solutions.\n";
		return NULL;
	}
	
	x = back_substitution(A, N, ops);
		
	return x;
}

double* gauss_elimination_penta(double **A, int N, size_t *ops, int verbose){
	int i, j, k, l = -1, flag = 1;
	double ratio, sum = 0.0, *x = NULL;
	
	if(verbose){	
		cout << "\nPivoted Matrix: \n";
		print_matrix<double>(A, N);
	}

	/*if(is_singular(A, N)){
		cout << "\nSingular matrix.\n";
		return NULL;
	}*/
	
	x = new double[N];

	for(k = 0; k < N - 1; ++k){
		for(i = k + 1; i < N; ++i){
			if(abs(i-k) > 2) continue;
			if(A[k][k] == 0){
				for(j = 0; j < k; ++j){
					if(l != -1 && A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				for(j = k+1; j < N; ++j){
					if(l != -1 && A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				if(l != -1){
					A = swap_lines(A, N, k, l);
					ops += 3;
				}
				l = -1;
			}
	
			ratio = A[i][k]/A[k][k];
	
			for(j = k + 1; j < N; ++j){
				if(abs(i-j) > 2) continue;
				A[i][j] -= ratio * A[k][j];
				*ops += 1;
			}
			A[i][N] -= ratio * A[k][N];
			*ops += 3;
		}
	}
	
	if(verbose){
		cout << "\nScalonated matrix:\n";
		print_matrix<double>(A, N);
	}

	for(i = 0; i < N; ++i){
		if(A[N-1][i] != 0){ flag = 0; break; }
	}	
	
	if(flag){
		cout << "\nInfinity solutions.\n";
		return NULL;
	}
	
	x = back_substitution_penta(A, N, ops);

	return x;
}
