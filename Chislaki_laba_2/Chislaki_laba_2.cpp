#include <iostream>
#include <iomanip>
using namespace std;

const double dx = 1e-9;
const int k_max = 50;//предельное число итераций

double* Solution_by_the_Gaussian_method(double** Matrix_A, double* Vector_b, int rank)
{
    double* Vector_X, max;
    int k, index;
    const double eps = 0.000001;  // точность
    Vector_X = new double[rank];
    k = 0;
    while (k < rank)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(Matrix_A[k][k]);
        index = k;
        for (int i = k + 1; i < rank; i++)
        {
            if (abs(Matrix_A[i][k]) > max)
            {
                max = abs(Matrix_A[i][k]);
                index = i;
            }
        }
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы Якоби" << endl;
            return 0;
        }
        // Перестановка строк
        for (int j = 0; j < rank; j++)
        {
            double temp = Matrix_A[k][j];
            Matrix_A[k][j] = Matrix_A[index][j];
            Matrix_A[index][j] = temp;
        }
        double temp = Vector_b[k];
        Vector_b[k] = Vector_b[index];
        Vector_b[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < rank; i++)
        {
            double temp = Matrix_A[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < rank; j++)
                Matrix_A[i][j] = Matrix_A[i][j] / temp;
            Vector_b[i] = Vector_b[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < rank; j++)
                Matrix_A[i][j] = Matrix_A[i][j] - Matrix_A[k][j];
            Vector_b[i] = Vector_b[i] - Vector_b[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = rank - 1; k >= 0; k--)
    {
        Vector_X[k] = Vector_b[k];
        for (int i = 0; i < k; i++)
            Vector_b[i] = Vector_b[i] - Matrix_A[i][k] * Vector_X[k];
    }
    return Vector_X;
    delete[]Vector_X;
}
//первое уравнение
double function1(double x1, double x2)
{
    return 1.5 * x1 * x1 * x1 - x2 * x2 - 1;
}
//второе уравнение
double function2(double x1, double x2)
{
    return x1*x2*x2*x2 - x2 - 4;
}
//вычисление производных аналитически
//первая функция
double Differential_function_1(double x1, double x2, int number_x)
{
    if (number_x == 0)
    {
        return 3 * 1.5 * x1 * x1;//по х1
    }
    else
        return -2 * x2;//по х2
}
//вторая функция
double Differential_function_2(double x1, double x2, int number_x)
{
    if (number_x == 0)
    {
        return x2*x2*x2;//по х1
    }
    else
        return x1*3*x2*x2 - 1;//по х2
}
//матрица Якоби конечно-разностным методом
double** find_Jacobian_NUMERIC(double x1, double x2)
{
    double** Jacobian = new double* [2];

    double M = 0.001;

    for (int i = 0; i < 1; i++)
    {
        Jacobian[i] = new double[2];

        Jacobian[i][0] = (function1(x1 + M * x1, x2) - function1(x1, x2)) / M * x1;
        Jacobian[i][1] = (function1(x1, x2 + M * x2) - function1(x1, x2)) / M * x2;        
    } 

    for (int i = 1; i < 2; i++)
    {
        Jacobian[i] = new double[2];

        Jacobian[i][0] = (function2(x1 + M * x1, x2) - function2(x1, x2)) / M * x1;
        Jacobian[i][1] = (function2(x1, x2 + M * x2) - function2(x1, x2)) / M * x2;
    }

    return Jacobian;
    for (int i = 0; i < 2; i++)
    {
        delete[]Jacobian[i];
    }
    delete[]Jacobian;
}
// матрица Якоби аналитическим методом
double** find_Jacobian_ANALITIC(double x1, double x2)
{
    double** Jacobian = new double*[2];
    
    for (int i = 0; i < 1; i++)
    {
        Jacobian[i] = new double[2];
        for (int j = 0; j < 2; j++)
        {
            Jacobian[i][j] = Differential_function_1(x1, x2, j);
        }        
    }
    for (int i = 1; i < 2; i++)
    {
        Jacobian[i] = new double[2];
        for (int j = 0; j < 2; j++)
        {
            Jacobian[i][j] = Differential_function_2(x1, x2, j);
        }
    }

    return Jacobian;
    for (int i = 0; i < 2; i++)
    {
        delete[]Jacobian[i];
    }
    delete[]Jacobian;
}

int main()
{
    setlocale(LC_ALL, "rus");
    double x01, x02; //исходное приближение
    cout << "Введите начальное приближение х1(0): ";
    cin >> x01;
    cout << "                            и х2(0): ";
    cin >> x02;
    int counter_iteration = 1;//количество итераций
    double x1, x2;//переменные
    
    x1 = x01;
    x2 = x02;

    double delta_1, delta_2;

    cout << "Номер итерации:            " << "delta_1:            " << "delta_2:" << setw(20) << x1 << setw(20)<< x2<< endl;
        
    //сохранение предыдущих значений вектора невязки
    double x1_pred = x1;
    double x2_pred = x2;

    //вектор невязки
    double* residual_vector = new double[2];
        
    double* vector_amendment = new double[2];//вектор поправки   

    do
    {
        double** Jacobian = find_Jacobian_ANALITIC(x1, x2);
        //double** Jacobian = find_Jacobian_NUMERIC(x1, x2);

        residual_vector[0] = -function1(x1, x2);
        residual_vector[1] = -function2(x1, x2);

        vector_amendment = Solution_by_the_Gaussian_method(Jacobian, residual_vector, 2);
        
        x1_pred = x1;
        x2_pred = x2;

        //изменение решения
        x1 = x1 + vector_amendment[0];
        x2 = x2 + vector_amendment[1];
                
        if (abs(function1(x1, x2)) >= abs(function2(x1, x2)))
            delta_1 = abs(function1(x1, x2));
        else delta_1 = abs(function2(x1, x2));

        if (abs(x1) < 1)
        {
            delta_2 = abs(x1 - x1_pred);
        }
        else
        {
            delta_2 = abs((x1 - x1_pred) / x1);
        }
        double save_delta_2 = delta_2;
        if (abs(x2) < 1)
        {
            delta_2 = abs(x2 - x2_pred);
        }
        else
        {
            delta_2 = abs((x2 - x2_pred) / x2);
        }
        if (delta_2 <= save_delta_2)
        {
            delta_2 = save_delta_2;
        }

        cout << setw(15) << counter_iteration << setw(20) << delta_1 << setw(20) << delta_2 << setw(20) << x1 << setw(20) <<  x2 << endl;

        double eps1 = 1e-9, eps2 = 1e-9;
        if ((delta_1 <= eps1) && (delta_2 <= eps2))
        {
            cout << "Ответ:" << endl
                << "x1 = " << x1 << "  x2 = " << x2 << endl;
            break;
        }
        else counter_iteration++;
    } while(counter_iteration < k_max);

    delete[] residual_vector;
    delete[] vector_amendment;
}