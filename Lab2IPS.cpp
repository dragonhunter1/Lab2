#include <stdio.h>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <chrono>
#include <iostream>

using namespace std::chrono;
using namespace std;

// ���������� ����� � �������� ���������� �������
const int MATRIX_SIZE = 5000;

/// ������� InitMatrix() ��������� ���������� � �������� 
/// ��������� ���������� ������� ���������� ����������
/// matrix - �������� ������� ����
void InitMatrix( double** matrix )
{
	for ( int i = 0; i < MATRIX_SIZE; ++i )
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for ( int i = 0; i < MATRIX_SIZE; ++i )
	{
		for ( int j = 0; j <= MATRIX_SIZE; ++j )
		{
			matrix[i][j] = rand() % 2500 + 1;
		}
	}
}

/// ������� SerialGaussMethod() ������ ���� ������� ������ 
/// matrix - �������� ������� �������������� ���������, �������� � ����,
/// ��������� ������� ������� - �������� ������ ������ ���������
/// rows - ���������� ����� � �������� �������
/// result - ������ ������� ����
double SerialGaussMethod( double **matrix, const int rows, double* result )
{
	int k;
	double koef;
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	// ������ ��� ������ ������
	for ( k = 0; k < rows; ++k )
	{
		//
		for ( int i = k + 1; i < rows; ++i )
		{
			koef = -matrix[i][k] / matrix[k][k];

			for ( int j = k; j <= rows; ++j )
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	// �������� ��� ������ ������
	high_resolution_clock::time_point t5 = high_resolution_clock::now();
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for ( k = rows - 2; k >= 0; --k )
	{
		result[k] = matrix[k][rows];

		//
		for ( int j = k + 1; j < rows; ++j )
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}
	high_resolution_clock::time_point t6 = high_resolution_clock::now();
	duration<double> duration1 = (t4 - t3);
	duration<double> duration2 = (t6 - t5);
	cout << "����� ���������� ������� ���� ������ ������  : " << duration1.count() << " ���" << endl;
	cout << "����� ���������� ��������� ���� ������ ������  : " << duration2.count() << " ���" << endl;
	return duration1.count();
}

double SerialGaussMethodParallel( double **matrix, const int rows, double* result )
{
	int k;
	double *koef= new double[rows];
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	// ������ ��� ������ ������
	for ( k = 0; k < rows; ++k )
	{
		//
		cilk_for ( int i = k + 1; i < rows; ++i )
		{
			koef[i] = -matrix[i][k] / matrix[k][k];

			for ( int j = k; j <= rows; ++j )
			{
				matrix[i][j] += koef[i] * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	// �������� ��� ������ ������
	high_resolution_clock::time_point t5 = high_resolution_clock::now();
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for ( k = rows - 2; k >= 0; --k )
	{
		result[k] = matrix[k][rows];

		cilk::reducer_opadd<int> msum(0);

		cilk_for ( int j = k + 1; j < rows; ++j )
		{
			*msum += matrix[k][j] * result[j];
		}
		result[k]=result[k]-msum.get_value();
		result[k] /= matrix[k][k];
	}
	high_resolution_clock::time_point t6 = high_resolution_clock::now();
	duration<double> duration1 = (t4 - t3);
	duration<double> duration2 = (t6 - t5);
	cout << "����� ���������� ������� ���� ������ ������  : " << duration1.count() << " ���" << endl;
	cout << "����� ���������� ��������� ���� ������ ������  : " << duration2.count() << " ���" << endl;
	//delete[]koef;
	return duration1.count();
}

int main()
{
	srand( (unsigned) time( 0 ) );

	int i;
	__cilkrts_set_param("nworkers", "8");
	setlocale(LC_ALL, "Russian");
	// ���-�� ����� � �������, ���������� � �������� �������
	const int test_matrix_lines = MATRIX_SIZE;

	double **test_matrix = new double*[test_matrix_lines];
	
	// ���� �� �������
	for ( i = 0; i < test_matrix_lines; ++i )
	{
		// (test_matrix_lines + 1)- ���������� �������� � �������� �������,
		// ��������� ������� ������� ������� ��� ������ ����� ���������, �������� � ����
		test_matrix[i] = new double[test_matrix_lines + 1];
	}

	// ������ ������� ����
	double *result = new double[test_matrix_lines];
	// ������������� �������� �������
	InitMatrix(test_matrix);
	//test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	//test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	//test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	//test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;
	cout << "���������������� ���������� : "  << endl;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	double time1 =SerialGaussMethod( test_matrix, test_matrix_lines, result );
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	cout << "������������ ���������� : "  << endl;
	high_resolution_clock::time_point t1s = high_resolution_clock::now();
	double time2 = SerialGaussMethodParallel( test_matrix, test_matrix_lines, result );
	high_resolution_clock::time_point t2s = high_resolution_clock::now();

	duration<double> duration = (t2 - t1);
	
	cout<<"��������� ������ ��������� : "<<time1/time2 << endl;


	for ( i = 0; i < test_matrix_lines; ++i )
	{
		delete[]test_matrix[i];
	}

	/*printf( "Solution:\n" );

	for ( i = 0; i < test_matrix_lines; ++i )
	{
		printf( "x(%d) = %lf\n", i, result[i] );
	}*/
	//cout << "����� ���������� ����� ������ ������ ( ���������������� ���������� ) : " << duration.count() << " ���" << endl;
	delete[] result;

	return 0;
}