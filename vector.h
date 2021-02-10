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
	Vector(int N);								//�����������	
	Vector(const Vector & v);					//2 ����������� ��� ����������
	~Vector();									//����������

	void setElement(int i, double a);			//����������� �������� �������� �������
	int get() const;							//������ ������� 
	double getElement(int i) const;				//���������� ������� �������
	void print_vector() const;					//����� ������� �� ����� 

	void Mult_on_scalar(double a, Vector &B);	//��������� �� �����

	void Sum(Vector &A, Vector &B);				//����� ��������

	double Norma();								//������ �������
	void Ort(Vector &B);						//��� �������

	double Scolar_mult( Vector &B);				//��������� ������������
	void Vector_mult(Vector &A, Vector &B);		//��������� ������������
	
	double angle_betweenV(Vector &A);			//���� ����� ���������
	
	Vector operator+(const Vector & v1);		//c�������
	Vector operator-(const Vector & v1);		//���������
	Vector &operator=(const Vector & v);		//������������
    double& operator*(const Vector & v1);		//��������� ������������
	Vector &operator^(const Vector & v);		//��������� ������������
	Vector operator*(double c);
	Vector operator-();
};
#endif