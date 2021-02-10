#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include "vector.h"
class Matrix
{
private:
	
	double **V;
	int a, b;
	int max(int a, int b);
	int min(int a, int b);

public:
	enum MatrixType { M_OK, M_SINGULAR };
	Matrix(int a, int b);					//конструктор  выделяет помять под матрицу, заполняет нулями
	Matrix(const Matrix & v);
	~Matrix();								//диструктор чистит память

	int geta() const;						//узнать высоту матрицы
	int getb() const;						//узнать ширину матрицы

	void setElement(int i,int j, double c);	//присвоить значение элементу
	double getElement(int i, int j) const;	//считать значение элемента
	void print_matrix() const;				//вывести матрицу на экран

	bool Sum_Matrix(Matrix &A, Matrix &B);	//сумма матриц

	void Mult_Matrix(Matrix &A, Matrix &B);	//умножение матриц
	void MultOn_Scolar(double c);			//умножение матрицы на число 

	void Transposed_Matrix(Matrix &A);		//транспанирование
	void getMinor(int i, int j, Matrix &A);	//получение минора

	double Det();							//определитель по методу диагонализации
	double Det1();							//по методу LU разложения

	/*решение СЛАУ через LU разложение
	действует на матрицу коэф, возвращает LU матрицу
	* *out - флаг корректности работы
	* X - вектор решения
	* f - вектор свободных коэф
	*/
	Matrix LU(MatrixType *out, Vector *X, Vector &f);	

	Matrix &operator=(const Matrix & v);
};
#endif // MATRIX_H