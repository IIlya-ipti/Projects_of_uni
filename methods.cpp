#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include<functional>
#include <cmath>
using namespace std;
typedef long long li;
double f1(double x, double y) {
    return 1.5 - cos(x) - y;
}
double f2(double x, double y) {
    return 1 - 2 * x + sin(y - 0.5);
}
double f11(double x, double y) {
    return sin(x);
}
double f12(double x, double y) {
    return -1;
}
double f21(double x, double y) {
    return -2;
}
double f22(double x, double y) {
    return cos(y - 0.5);
}
double det(vector<vector<double>> matrix) {
    // определитель двухмерной матрицы
    return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
}
vector<vector<double>> matrix_inv(vector<vector<double>>& F) {
    double det_int = det(F);
    // транспонирование
    double k = F[0][1];
    F[0][1] = F[!0][!1];
    F[!0][!1] = k;
    vector<vector<double>> matrix(F.size(), vector<double>(F.size()));
    // обратная матрица 
    for (int i = 0; i < F.size(); i++) {
        for (int j = 0; j < F.size(); j++) {
            matrix[i][j] = (pow(-1, i + j + 2) * F[!i][!j]) / det_int;
        }
    }
    return matrix;
}
vector<vector<double>> inv(vector<double>& x) {
    // матрица Якоби
    // для частного случая (для курсовичка, поэтому i != j)
    int n = x.size();
    vector<vector<double>> matrix(n, vector<double>(n));
    matrix[0][0] = f11(x[0], x[1]);
    matrix[0][1] = f12(x[0], x[1]);
    matrix[1][0] = f21(x[0], x[1]);
    matrix[1][1] = f22(x[0], x[1]);
    matrix = matrix_inv(matrix);
    return matrix;
}
vector<double> numbers(vector<vector<double>> matrix, vector<function<double(double, double)>>& f, vector<double>& x) {
    // вычисление нового x
    vector<double> new_x(x.size());
    new_x[0] = f[0](x[0], x[1]);
    new_x[1] = f[1](x[0], x[1]);
    vector<double> y(x.size(), 0);
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < x.size(); j++) {
            y[i] += matrix[i][j] * new_x[j];
        }
    }
    return { x[0] - y[0], x[1] - y[1] };
}
void method_of_Newton(vector<function<double(double, double)>>& F, vector<double> x0) {
    // метод Ньютона
    vector<double>x(2);
    x = numbers(inv(x0), F, x0);
    double eps = 0.0001;
    int k = 1;
    cout << " x = " << x[0] << " y = " << x[1] << " f1 =  " << f1(x[0], x[1]) << " f2 = " << f2(x[0], x[1]) << " k = " << k << endl;
    while ((abs(x[0] - x0[0]) > eps) || (abs(x[1] - x0[1]) > eps)) {
        if (k > 10000)break;
        x0 = x;
        x = numbers(inv(x0), F, x0);
        ++k;
        cout << " x = " << x[0] << " y = " << x[1] << " f1 =  " << f1(x[0], x[1]) << " f2 = " << f2(x[0], x[1]) << " k = " << k << endl;
    }
    cout << " x = " << x[0] << " y = " << x[1] << " k = " << k << endl;
}
double F1(double x) {
    return 1.5 - cos(x);
}
double F2(double y) {
    return (0.5 + 0.5 * sin(y - 0.5));
}
void prime_iterations(vector<function<double(double)>>& F, vector<double> x0) {
    // метод простых итераций
    vector<double>x(2);
    double eps = 0.0001;
    int k = 1;
    x[0] = F2(x0[1]);
    x[1] = F1(x0[0]);
    cout << " x = " << x[0] << " y = " << x[1] << " f1 =  " << f1(x[0], x[1]) << " f2 = " << f2(x[0], x[1]) << " k = " << k << endl;
    while ((abs(x[0] - x0[0]) > eps) || (abs(x[1] - x0[1]) > eps)) {
        x0 = x;
        x[0] = F2(x0[1]);
        x[1] = F1(x0[0]);
        ++k;
        cout << " x = " << x[0] << " y = " << x[1] << " f1 =  " << f1(x[0], x[1]) << " f2 = " << f2(x[0], x[1]) << " k = " << k << endl;
    }
    cout << " x = " << x[0] << " y = " << x[1] << " f1 =  " << f1(x[0], x[1]) << " f2 = " << f2(x[0], x[1]) << " k = " << k << endl;
    cout << "x - x0 = " << abs(x[0] - x0[0]) << " Y -  Y0 = " << abs(x[1] - x0[1]) << endl;
}
int main()
{
    vector<function<double(double, double)>> F;
    F.push_back(f1);
    F.push_back(f2);
    vector<double>x0(2);
    vector<double>x(2);
    x0[0] = 1.5;
    x0[1] = 1.5;
    method_of_Newton(F, x0);

    return 0;
}

