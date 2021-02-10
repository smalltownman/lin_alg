#define _CRT_SECURE_NO_WARNINGS
#include "Matrix.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "Vector.h"
#include <omp.h>

const char gnuplot[] = "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist";

//теплопроводность
double lambda(double T);
//функция, расписанная через конечные разности
void F(double T0, double T1, Vector *T, Vector *val);
//якобиан
void jacobi(double T0, double T1, Vector* T, Vector *val, Vector *v1, Vector *v2, Vector *v3);
void Jacobi(double T0, double T1, Vector *T, Vector *val, Matrix *Val);

int max(int a, int b);
int min(int a, int b);

//метод прогонки
void solveMatrix(Vector *x, Vector b, Vector c, Vector a, Vector f);

int main(int argc, char * argv[]) {
	
	int N = 7;
	Vector T(N - 2);

	double T0 = 300.0, T1 = 1500.0, dx = 1.0 / (N - 1), E = 1e-4, dT = (T1 - T0) / (N - 1);
	int i, j;

	//массив координат
	Vector X(N);
	Vector T2(N);
	for (i = 0; i < N; i++) {
		X.setElement(i, dx * i);
		T2.setElement(i, T0 + i * dT);
	}

	//начальные данные для вычисления
	for (i = 0; i < N - 2; i++) {
		T.setElement(i, T0 + (i + 1) * dT);
	}
	
	Vector val(N - 2);
	Matrix L_U(N - 2, N - 2);
	Matrix::MatrixType out;
	Matrix Val(N - 2, N - 2);
	Vector H(N - 2); 
	
	clock_t begin = clock();

	//метод Ньютона
	F(T0, T1,&T, &val);

	int y = 0;
	Vector v1(N - 2), v2(N - 2), v3(N - 2);

	while (val.Norma() > E) {
		jacobi(T0, T1, &T, &val, &v1, &v2, &v3);
		//Jacobi(T0, T1, &T, &val, &Val);
		
		
		//L_U = Val.LU(&out, &H, val);
		//if (out) {
		//	return 0;
		//}
		
	
		solveMatrix(&H, v1, v2, v3, val);
		
		T = H - T; //T(k + 1);
		F(T0, T1, &T, &val);
		y += 1;
	}

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time = %g\n", time_spent);
	printf("y = %d\n", y);
	Vector newT(N);
	
	#pragma region  writing to file

	double a = -19.0 / 360000.0, b = 53.0 / 300.0, c1 = -185.0 / 4.0;
	double C2 = a * T0 * T0 * T0 / 3 + b / 2 * T0 * T0 + c1*T0;
	double C1 = a * T1 * T1 * T1 / 3 + b / 2 * T1 * T1 + c1 * T1 - C2;
	double x = 0.0;
	
	FILE * f = fopen("one.txt", "w");

	for (i = 0; i < N; i++) {
		if (i == 0) {
			newT.setElement(i, T0);
		}
		else if (i == N - 1) {
			newT.setElement(i, T1);
		}
		else {
			newT.setElement(i, T.getElement(i - 1));
		}
		x = (a / 3 * T2.getElement(i) * T2.getElement(i) * T2.getElement(i) + b / 2\
			* T2.getElement(i) *T2.getElement(i) + c1 * T2.getElement(i) - C2) / C1;
		fprintf(f, "%16.15lg\t%16.15lg\t%16.15lg\t%16.15lg\n", X.getElement(i), newT.getElement(i), x, T2.getElement(i));
	}

	fclose(f);
	#pragma endregion
	
	#pragma region drawing a graph
	
		
	FILE * pipe = _popen(gnuplot, "w");
	fprintf(pipe,
		"set terminal pngcairo size 15cm, 10cm font \'Times, 10\'\n"
		"set output \'curve1.png\'\n"
		"set yrange [300: 1500]\n"
		"set xrange [0: 1]\n"
		"set grid lw 1\n"
		"plot	\'one.txt\' using 1:2 title \'pr\' w l lc rgb\'#FA8072\',\\\n"
					" \'one.txt\' using 3:4 title \'th\' w l lc rgb\'#FA8000\'\n");
	_pclose(pipe);
#pragma endregion

	return 0;
}

double lambda(double T) {
	double a = - 19.0 / 360000.0, b = 53.0 / 300.0, c = -185.0 / 4.0;
	return a * T * T + b * T + c;
}

void F(double T0, double T1, Vector *T, Vector *val) {
	
	int N = T->get();
	double c = 0.0;
	for (int i = 0; i < N; i++) {
		if (i == 0) {
			  c = lambda((T->getElement(i + 1) + T->getElement(i)) / 2)*(T->getElement(i + 1) - T->getElement(i)) - lambda((T0+T->getElement(i))/2)*(T->getElement(i)-T0);
			  val->setElement(i, c);
		}
		else if (i == N - 1) {
			c = lambda((T1 + T->getElement(i)) / 2)*(T1 - T->getElement(i)) - lambda((T->getElement(i-1)+T->getElement(i))/2)*(T->getElement(i)-T->getElement(i-1));
			val->setElement(i, c);
		}
		else { 
			c = lambda((T->getElement(i + 1) + T->getElement(i)) / 2)*(T->getElement(i + 1) - T->getElement(i)) -  lambda((T->getElement(i-1)+ T->getElement(i))/2)*(T->getElement(i)-T->getElement(i-1));
			val->setElement(i, c);
		}
	}
	
}

void jacobi( double T0, double T1, Vector *T, Vector *val, Vector *v1, Vector *v2, Vector *v3) {
	double dx = 1e-4, c = 0.0;
	int N = T->get(), k = 0;
	Vector Val1(N);
	bool a1, a2, a3;

	for (int i = 0; i < N; i++) {

		//численная производная
		T->setElement(i, T->getElement(i) + dx);
		F(T0, T1, T, &Val1);
		T->setElement(i, T->getElement(i) - dx);
		k = 0;

		for (int j = max(0, i - 1); j < min(i + 2, N); j++) {

			c = (Val1.getElement(j) - val->getElement(j)) / dx;

			a1 = (k == 0 && i > 0);
			a2 = ((k == 1 && i > 0) || (i == 0 && k == 0));
			a3 = (k == 2 || (i == 0 && k == 1));
		
			if (a1) {
				v1->setElement(i - 1, c);
			}
			else if (a2) {
				v2->setElement(i, c);
			}
			else if (a3) {
				v3->setElement(i + 1, c);
			}
			k += 1;
		}
	}
}

void Jacobi(double T0, double T1, Vector *T, Vector *val, Matrix *Val) {
	double dx = 1e-4, c = 0.0;
	int n = T->get(), k = 0;
	Vector Val1(n);

	for (int i = 0; i < n; i++) {


		T->setElement(i, T->getElement(i) + dx);
		F(T0, T1, T, &Val1);
		T->setElement(i, T->getElement(i) - dx);


		for (int j = max(0, i - 1); j < min(i + 2, n); j++) {

			c = (Val1.getElement(j) - val->getElement(j)) / dx;


			for (int j = max(0, i - 1); j < min(i + 2, n); j++) {

				c = (Val1.getElement(j) - val->getElement(j)) / dx;
				Val->setElement(j, i, c);

			}
		}
	}
}


int max(int a, int b) {
	int Max = a;
	if (Max <= b) {
		Max = b;
	}
	return Max;

}

int min(int a, int b) {
	int Min = a;
	if (Min >= b) {
		Min = b;
	}
	return Min;
}

/*
* b - диагональ, лежащая над главной 
* c - главная диагональ матрицы A 
* a - диагональ, лежащая под главной 
* f - правая часть 
* x - решение
*/
void solveMatrix(Vector *x, Vector b, Vector c, Vector a, Vector f){

	int n = x->get(), i;
	Vector p(n);
	Vector q(n);
	double m = 0.0;
	
	p.setElement(0, -b.getElement(0) / c.getElement(0)); 
	q.setElement(0, f.getElement(0) / c.getElement(0));
	 
	for (i = 1; i < n; i++) {
		p.setElement(i, -b.getElement(i) / (c.getElement(i) + a.getElement(i)*p.getElement(i - 1)));
		q.setElement(i, (f.getElement(i) - a.getElement(i) * q.getElement(i - 1) ) / (c.getElement(i) + a.getElement(i) * p.getElement(i - 1)));
	}

	x->setElement(n - 1, q.getElement(n - 1));
	
	for (i = n - 2; i >= 0; i--) {
		m = p.getElement(i) * x->getElement(i + 1) + q.getElement(i);
		x->setElement(i, m);
	}
}