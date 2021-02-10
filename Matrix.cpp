#include "Matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


Matrix::Matrix(int a = 0, int b = 0)
{ 
	this->a = a;
	this->b = b;
	if (a == 0 || b == 0) {
		V = nullptr;
	}
	else {
		V = (double**)malloc(sizeof(double*)*a);
		for (int i = 0; i < a; i++) {
			V[i] = (double*)malloc(sizeof(double)*b);
			for (int j = 0; j < b; j++) {
				V[i][j] = 0.0;
			}
		}
	}
}

int Matrix::max(int a, int b) {
	int Max = a;
	if (Max <= b) {
		Max = b;
	}
	return Max;

}

int Matrix::min(int a, int b) {
	int Min = a;
	if (Min >= b) {
		Min = b;
	}
	return Min;
}


Matrix::Matrix(const Matrix &v) {
	 a = v.a;
	 b = v.b;
	V = (double**)malloc(sizeof(double*) * a);
	for (int i = 0; i < a; i++) {
		V[i] = (double*)malloc(sizeof(double)*b);
		for (int j = 0; j < b; j++) {
			V[i][j] = v.V[i][j];
		}
		
	}
}

Matrix::~Matrix()
{
	for (int i = 0; i < a; i++) {
		free(V[i]);
	}
	free(V);
}


int Matrix::geta() const {
	return a;
}


int Matrix::getb() const {
	return b;
}


void Matrix::setElement(int i, int j, double c) {
	V[i][j] = c;
}


double Matrix::getElement(int i, int j) const {
	return V[i][j];
}


void Matrix::print_matrix() const {
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < b; j++) {
			printf("%6.5lg\t", V[i][j]);
		}
		printf("\n");
	}
}


bool Matrix::Sum_Matrix(Matrix &A, Matrix &B) {
	if (A.geta() == B.geta() &&  B.getb() == A.getb()) {
		int i, j;
		for (i = 0; i < A.geta(); i++) {
			for (j = 0; j < A.getb(); j++) {
				V[i][j] = A.getElement(i, j) + B.getElement(i, j);
			}
		}
		return true;
	}
	else {
		printf("matrix can't be add\n");
		return false;
	}
}


void Matrix::Mult_Matrix(Matrix &A, Matrix &B) {
	if (a == A.geta() && b == B.getb()) {
		if (A.getb() == B.geta()) {
			int i, j, k;
			double c = 0.0;
			for (i = 0; i < A.geta(); i++) {
				for (j = 0; j < B.getb(); j++) {
					c = 0;
					for (k = 0; k < B.geta(); k++) {
						c += A.getElement(i, k) * B.getElement(k, j);
					}
					V[i][j] = c;
				}
			}
		}
		else {
			printf("Matrix can't be multiple\n");
		}
	}
	else {
		printf("new wrong size matrix\n");
	}
}

void Matrix::MultOn_Scolar(double c) {
	int i, j;
	for (i = 0; i < a; i++) {
		for (j = 0; j < b; j++) {
			V[i][j] = V[i][j] * c;
		}
	}
}

void Matrix::Transposed_Matrix(Matrix &A) {
	int i, j;
	for (i = 0; i < b; i++) {
		for (j = 0; j < a; j++) {
			//V[i][j] = A.getElement(j, i);
			V[i][j] = A.V[i][j];
		}
	}
}

void Matrix::getMinor(int i, int j, Matrix &A) {
	int y, x;
	int di = 0, dj = 0;
	for (x = 0; x < b; x++) {
		if (x == i) {
			di = 1;
		}
		dj = 0;
		for (y = 0; y < b; y++) {
			if (y == j) {
				dj = 1;
			}
			V[x][y] = A.getElement(x + di, y + dj);
		}
	}
}



double Matrix::Det() { 
	
	if (a == b) {
		

		int i, j, k, c = 0;
		double det = 1.0,  d = 0.0, x = 0.0, R;
		
		Matrix G(a, a);
		//скопируем матрицу, чтобы не диагонализировать исходную
		for (i = 0; i < a; i++) {
			for (j = 0; j < a; j++) {
				G.V[i][j] = V[i][j];
			}
		}


		for (k = 0; k <= a - 2; k++) {
			//если диагональный элемент равен 0
			if (fabs(G.V[k][k]) < 1e-14) {
				d = G.V[k][k];
				c = k;
				while ((fabs(d) < 1e-14) && (c < a - 1)) {
					c += 1;
					d = G.V[c][k];
				}
				if (fabs(d) < 1e-14) { //если есть нулевой столбец det = 0
					printf("degenerate  matrix\n");
					return det = 0.0;
				}
				else {
					det *= pow(-1, c); //замена строк местами с учётом чётности номера строки 
					for (j = 0; j < a; j++) {
						x = G.V[k][j];
						G.V[k][j] = G.V[c][j];
						G.V[c][j] = x;
					}
				}
			}
			//обнуление столбца
			for (i = k + 1; i < a; i++) {
				R = G.V[i][k];
				for (j = k; j < a; j++) {
					G.V[i][j] = G.V[i][j]- R * G.V[k][j] / G.V[k][k];
				}
			}
		}
		//перемножение диагональных членов
		for (k = 0; k < a; k++) {
			det *= G.V[k][k];
		}
		printf("\n");
		G.print_matrix();
		printf("\n\n");
	
		return det;
		
		
	}
	else {
		printf("matrix not scuare\n");
		return 0.0;
	}
}

Matrix Matrix::LU(MatrixType *out, Vector *X, Vector &f) {

	Matrix G(a, a);
	int k, i, j, d;
	double sum = 0.0;

	for (k = 0; k < a; k++) {
		for (j = k; j < min(k + 2, a); j++) {

			sum = 0.0;
			for (d = 0; d < k; d++) {
				sum += G.V[k][d] * G.V[d][j];
			}
			G.V[k][j] = V[k][j] - sum;
		}

		for (i = k + 1; i < min(k + 2, a); i++) {

			if (fabs(G.V[k][k]) < 1e-22) {
				printf("degenerate system\n");
				*out = M_SINGULAR;
				return Matrix();
			}
			sum = 0.0;
			for (d = 0; d < k; d++) {
				sum += G.V[i][d] * G.V[d][k];
			}
			G.V[i][k] = (V[i][k] - sum) / G.V[k][k];
		}
	}

	Vector Y(a);
	double c = 0.0;
	for (i = 0; i < a; i++) {
		sum = 0.0;

		for (j = 0; j <= i - 1; j++) {
			if (i == j) {
				c = 1.0;
			}
			else {
				c = G.V[i][j];
			}
			sum += c * Y.getElement(j);
		}

		Y.setElement(i, f.getElement(i) - sum);
		
	}
	for (i = a - 1; i >= 0; i--) {
		sum = 0.0;

		for (j = a - 1; j >= i + 1; j--) {
			sum += G.getElement(i, j)  * X->getElement(j);
		}

		c = (Y.getElement(i) - sum) / G.getElement(i, i);
		X->setElement(i, c);
	}
	*out = M_OK;
	
	return G;
}

double Matrix::Det1() {
	double det = 1.0;
    
	Matrix L(a, a);
	Matrix U(a, a);
	
	int k = 0, i = 0, j = 0, d = 0;
	double sum1 = 0.0, sum2 = 0.0;

	for (k = 0; k < a; k++) {
	
		for (j = k; j < a; j++) {
			
				sum1 = 0.0;
				for (d = 0; d < k ; d++) {
					sum1 += L.V[k][d] * U.V[d][j];
				}

				U.V[k][j] = V[k][j] - sum1;
			}
			
		for (i = k + 1; i < a; i++) {

			if (fabs(U.V[k][k]) < 1e-14) {
				printf("degenerate system\n");
				return 0.0;
			}

			sum2 = 0.0;
			for (d = 0; d < k; d++) {
				sum2 += L.V[i][d] * U.V[d][k];
			}

			L.V[i][k] = (V[i][k] - sum2) / U.V[k][k];
		}
	}
	
	for (i = 0; i < a; i++) {
		det *= U.V[i][i];
	}
	return det;
}

Matrix &Matrix ::operator=(const Matrix & v) {
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < b; j++) {
			V[i][j] = v.V[i][j];
		}
	}
	return *this;
}