#pragma once 
#ifndef VECTOR_H
#define VECTOR_H

class Vector
{
private:
	double * V;
	int N;
public:
	enum VectorType { V_OK, V_SINGULAR, V_IVCORRECT_DIMENSION};
	Vector(int N);								//конструктор	
	Vector(const Vector & v);					//2 конструктор для перегрузок
	~Vector();									//диструктор

	void setElement(int i, double a);			//присваивает значение элементу массива
	int get() const;							//размер массива 
	double getElement(int i) const;				//возвращает элемент массива
	void print_vector() const;					//вывод вектора на экран 

	void Mult_on_scalar(double a, Vector &B);	//умножение на число

	void Sum(Vector &A, Vector &B);				//сумма векторов

	double Norma();								//модуль вектора
	void Ort(Vector &B);						//орт вектора

	double Scolar_mult( Vector &B);				//сколярное произведение
	void Vector_mult(Vector &A, Vector &B);		//векторное произведение
	
	double angle_betweenV(Vector &A);			//угол между векторами
	
	Vector operator+(const Vector & v1);		//cложение
	Vector operator-(const Vector & v1);		//вычетание
	Vector &operator=(const Vector & v);		//присваивание
    double& operator*(const Vector & v1);		//скалярное произведение
	Vector &operator^(const Vector & v);		//векторное произведение
	Vector operator*(double c);
	Vector operator-();
};
#endif