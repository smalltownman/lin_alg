#include "Vector.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

Vector::Vector(int N = 0)//конструктор
{
	this->N = N;
	if (N == 0) {
		V = nullptr;
	}
	else {
		V = (double*)malloc(sizeof(double)*N);
		for (int i = 0; i < N; i++) {
			V[i] = 0.0;
		}
	}
	
	
}

Vector::Vector(const Vector & v) {
	N = v.N;
	V = (double*)malloc(sizeof(double) * N);
	for (int i = 0; i < N; i++) {
		V[i] = v.V[i];
	}
}

Vector::~Vector()//диструктор
{
	free(V);
}

void Vector::setElement(int i, double с) { //присвоение значения ячейке массива
	  V[i] = с;
}
int Vector::get() const {
	return N;
}
double Vector::getElement(int i) const {
	return V[i];
}
void Vector::Mult_on_scalar(double a, Vector &B) { //умножение на число 
	for (int i = 0; i < N; i++) {
		V[i] = B.getElement(i) * a;
	}
} 

void Vector::Sum(Vector &A, Vector &B) { //сумма векторов
	for (int i = 0; i < N; i++) { 
		V[i]= A.getElement(i) + B.getElement(i);
	}
}
double Vector::Norma() {//норма
	double  c = 0.0;
	for (int i = 0; i < N; i++) {
		c += V[i] * V[i];
	}
	return sqrt(c);
}
void Vector::print_vector() const { //вывод вектора на экран
	for (int i = 0; i < N; i++) {
		printf("%8.6lg\t", V[i]);
	}
	printf("\n");
}
double Vector::Scolar_mult(Vector &B) {//сколярное произведение 
	double c = 0.0;
	for (int i = 0; i < N; i++) {
		c += V[i] * B.getElement(i);
	}
	return c;
}
void Vector::Vector_mult(Vector &A, Vector &B) {//векторное произведение 
	if (N == 3) {

		double a_x = A.getElement(0);
		double a_y = A.getElement(1);
		double a_z = A.getElement(2);
		double b_x = B.getElement(0);
		double b_y = B.getElement(1);
		double b_z = B.getElement(2);
		V[0] = a_y * b_z - a_z * b_y;
		V[1] = a_z * b_x - a_x * b_z;
		V[2] = a_x * b_y - a_y * b_x;
	}
		else {
			printf("Can't be multiplied");
		}
}
void Vector::Ort(Vector &B) { //орт
	double c = B.Norma();
	for (int i = 0; i < N; i++) {
		V[i] = B.getElement(i) / c;
	}
}

double Vector::angle_betweenV(Vector &A) { //угол между векторами
	if (N <= 3) {
		double  c = 0.0;
		for (int i = 0; i < N; i++) {
			c += V[i] * V[i];
		}
		c = sqrt(c);
		double c2 = A.Norma();
		double c3 = 0.0;
		for (int i = 0; i < N; i++) {
			c3 += V[i] * A.getElement(i);
		}
		double cos = c3 / (c * c2);
		return cos;
	}
	else {
		printf("incorrect vectors\n");
		return 0;
	}
}

Vector Vector::operator+(const Vector & v1){
	Vector Ans(v1);
	for (int i = 0; i < N; i++){
		Ans.V[i] += V[i];
	}
	return Ans;
}

Vector Vector::operator-(const Vector & v1) {
	Vector Ans(v1);
	for (int i = 0; i < N; i++) {
		Ans.V[i] -= V[i];
	}
	return Ans;
}

Vector & Vector::operator=(const Vector & v){
	for(int i = 0; i < N; i++){
		V[i] = v.V[i];
	}
	return * this;
}

double& Vector::operator*(const Vector & v1) {
	double mult = 0.0;
	Vector Ans(v1);
	for (int i = 0; i < N; i++) {
		mult += Ans.V[i] * V[i];
	}
	return mult;
}


Vector & Vector::operator^(const Vector & v) { 
	if (N == 3 && v.N == 3) {
		Vector Ans(v);
		Ans.V[0] = V[1] * v.V[2] - V[2] * v.V[1];
		Ans.V[1] = V[2] * v.V[0] - V[0] * v.V[2];
		Ans.V[2] = V[0] * v.V[1] - V[1] * v.V[0];
		return Ans;
	}
	else {
		return Vector();
	}
}

Vector Vector::operator*(double c) { 
	Vector ans(N);
	for (int i = 0; i < N; i++) {
		ans.V[i] = c * V[i];
	}
	return ans;
}

Vector Vector::operator-() {
	Vector Ans(N);
	for (int i = 0; i < N; i++) {
		Ans.V[i] = -1.0 * V[i];
	}
	return Ans;
}