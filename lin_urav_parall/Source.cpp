#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#define EPS 0.1 //необходимая точность искомого решения СЛАУ
//#define SIZE 1500 //максимальный размер СЛАУ

using namespace std;

const int SIZE = 3000;
/*
Функция, возвращающая матрицу размером n на n+1
n+1-ый столбик свободные члены СЛАУ
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
Процедура записи в файл массива с затратами
времени на поиск решения СЛАУ с размером
от 3х3 до 1000х1000.
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
Функция, реализующая метод Якоби (простой итерации) 
для вычисления решения СЛАУ
*/
void YakobyIterationMethod(vector <vector<long double>> matrix, int isParallel, string fileName) 
{
	int ti=0;
	bool isTheEnd = false;
	time_t t1,t2;
	vector <long double> previousVariableValues;
	vector <long double> currentVariableValues;
	vector <double> timeArray(SIZE - 2, 0.0);
	// массив с приближенным решением предыдущей итерации	
	for (int size_i = 2; size_i < SIZE; size_i++)
	{
		previousVariableValues = vector <long double>(size_i, 0.0);
		// Будем выполнять итерационный процесс до тех пор, 
		// пока не будет достигнута необходимая точность    
		t1 = clock();
		
		while (!isTheEnd)
		{
			// Введем вектор значений неизвестных на текущем шаге       
			currentVariableValues = vector <long double>(size_i);
			int j;
			// Посчитаем значения неизвестных на текущей итерации
			// в соответствии с теоретическими формулами
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
				// Делим на коэффициент при i-ой неизвестной
				currentVariableValues[i] /= matrix[i][i];
			}

			// Посчитаем текущую погрешность относительно предыдущей итерации
			long double error = 0.0;
			for (int i = 0; i < size_i; i++)
			{
				error += abs(currentVariableValues[i] - previousVariableValues[i]);
			}
			// Если необходимая точность достигнута, то завершаем процесс
			if (error < EPS)
			{
				t2 = clock();
				timeArray[size_i - 2] = ((double)(t2 - t1)) / (double)CLOCKS_PER_SEC;
				if((size_i % 100) == 0)
				{
					cout << "Время = " << timeArray[size_i - 2];
					cout << ", найдено решение для СЛАУ размером " << size_i << endl;
					if ((size_i % 1000) == 0) {
						writeTime(timeArray, fileName);
						isTheEnd = true;
					}
				}
				break;
			}
			// Переходим к следующей итерации, так 
			// что текущие значения неизвестных 
			// становятся значениями на предыдущей итерации
			previousVariableValues = currentVariableValues;
		}
		if (isTheEnd) break;
	}
}

int main()
{
	setlocale(LC_ALL, "ru");
	
	// Будем хранить матрицу в векторе, состоящем из 
	// векторов вещественных чисел
	// Матрица будет иметь размер (size) x (size + 1),
	// c учетом столбца свободных членов
	vector <vector <long double> > matrixExample = createMatrix(SIZE);
	cout << "С распараллеливанием" << endl;
	YakobyIterationMethod(matrixExample, 1, "result_time1.txt");
	cout << "Без распараллеливания" << endl;
	YakobyIterationMethod(matrixExample, 0, "result_time2.txt");
	cout << endl;
	system("pause");
	return 0;
}