#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#define EPS 0.1 //����������� �������� �������� ������� ����
//#define SIZE 1500 //������������ ������ ����

using namespace std;

const int SIZE = 3000;
/*
�������, ������������ ������� �������� n �� n+1
n+1-�� ������� ��������� ����� ����
*/
vector <vector<long double>> createMatrix(int n) 
{
	srand(time(0));
	vector <vector <long double> > matrix(n, vector<long double>(n+1, 0.0));
	long double summ=0.0;
	for (size_t i = 0; i < n; i++)
	{
		summ = 0;
		for (size_t j = 0; j < n+1; j++)
		{
			matrix[i][j] = rand();
			summ += matrix[i][j];
		}
		matrix[i][i] = summ;
	}
	return matrix;
}

/*
��������� ������ � ���� ������� � ���������
������� �� ����� ������� ���� � ��������
�� 3�3 �� 1000�1000.
*/
void writeTime(vector<double> t, string fileName) {
	srand(time(0));

	ofstream fout(fileName, ios_base::out | ios_base::trunc);

	for (int i = 0; i < SIZE-2; i++)
	{
		fout << t[i] << endl;
	}
	fout.close();
}

/*
�������, ����������� ����� ����� (������� ��������) 
��� ���������� ������� ����
*/
void YakobyIterationMethod(vector <vector<long double>> matrix, int isParallel, string fileName) 
{
	int ti=0;
	bool isTheEnd = false;
	time_t t1,t2;
	vector <long double> previousVariableValues;
	vector <long double> currentVariableValues;
	vector <double> timeArray(SIZE - 2, 0.0);
	// ������ � ������������ �������� ���������� ��������	
	for (int size_i = 2; size_i < SIZE; size_i++)
	{
		previousVariableValues = vector <long double>(size_i, 0.0);
		// ����� ��������� ������������ ������� �� ��� ���, 
		// ���� �� ����� ���������� ����������� ��������    
		t1 = clock();
		
		while (!isTheEnd)
		{
			// ������ ������ �������� ����������� �� ������� ����       
			currentVariableValues = vector <long double>(size_i);
			int j;
			// ��������� �������� ����������� �� ������� ��������
			// � ������������ � �������������� ���������
#pragma omp parallel for private (j) if(isParallel)
			for (int i = 0; i < size_i; i++)
			{
				currentVariableValues[i] = matrix[i][size_i];
				for (j = 0; j < size_i; j++)
				{
					if (i != j)
					{
						currentVariableValues[i] -= matrix[i][j] * previousVariableValues[j];
					}
				}
				// ����� �� ����������� ��� i-�� �����������
				currentVariableValues[i] /= matrix[i][i];
			}

			// ��������� ������� ����������� ������������ ���������� ��������
			long double error = 0.0;
			for (int i = 0; i < size_i; i++)
			{
				error += abs(currentVariableValues[i] - previousVariableValues[i]);
			}
			// ���� ����������� �������� ����������, �� ��������� �������
			if (error < EPS)
			{
				t2 = clock();
				timeArray[size_i - 2] = ((double)(t2 - t1)) / (double)CLOCKS_PER_SEC;
				if((size_i % 100) == 0)
				{
					cout << "����� = " << timeArray[size_i - 2];
					cout << ", ������� ������� ��� ���� �������� " << size_i << endl;
					if ((size_i % 1000) == 0) {
						writeTime(timeArray, fileName);
						isTheEnd = true;
					}
				}
				break;
			}
			// ��������� � ��������� ��������, ��� 
			// ��� ������� �������� ����������� 
			// ���������� ���������� �� ���������� ��������
			previousVariableValues = currentVariableValues;
		}
		if (isTheEnd) break;
	}
}

int main()
{
	setlocale(LC_ALL, "ru");
	
	// ����� ������� ������� � �������, ��������� �� 
	// �������� ������������ �����
	// ������� ����� ����� ������ (size) x (size + 1),
	// c ������ ������� ��������� ������
	vector <vector <long double> > matrixExample = createMatrix(SIZE);
	cout << "� ������������������" << endl;
	YakobyIterationMethod(matrixExample, 1, "result_time1.txt");
	cout << "��� �����������������" << endl;
	YakobyIterationMethod(matrixExample, 0, "result_time2.txt");
	cout << endl;
	system("pause");
	return 0;
}