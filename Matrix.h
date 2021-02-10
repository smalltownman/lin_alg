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
	Matrix(int a, int b);					//�����������  �������� ������ ��� �������, ��������� ������
	Matrix(const Matrix & v);
	~Matrix();								//���������� ������ ������

	int geta() const;						//������ ������ �������
	int getb() const;						//������ ������ �������

	void setElement(int i,int j, double c);	//��������� �������� ��������
	double getElement(int i, int j) const;	//������� �������� ��������
	void print_matrix() const;				//������� ������� �� �����

	bool Sum_Matrix(Matrix &A, Matrix &B);	//����� ������

	void Mult_Matrix(Matrix &A, Matrix &B);	//��������� ������
	void MultOn_Scolar(double c);			//��������� ������� �� ����� 

	void Transposed_Matrix(Matrix &A);		//����������������
	void getMinor(int i, int j, Matrix &A);	//��������� ������

	double Det();							//������������ �� ������ ��������������
	double Det1();							//�� ������ LU ����������

	/*������� ���� ����� LU ����������
	��������� �� ������� ����, ���������� LU �������
	* *out - ���� ������������ ������
	* X - ������ �������
	* f - ������ ��������� ����
	*/
	Matrix LU(MatrixType *out, Vector *X, Vector &f);	

	Matrix &operator=(const Matrix & v);
};
#endif // MATRIX_H