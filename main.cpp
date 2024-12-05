#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <utility>

using namespace std;

/// 高斯核密度函数 https://zhuanlan.zhihu.com/p/643143096
/// 参数x ： 样本数据点
/// 参数mean ： 所有像素点
/// 参数h ： 带宽 
double gaussian(double x, double mean, double h) {
    double constant_1 =  1.0 / (h * sqrt(2*M_PI));
    double exponent = -0.5 * pow((x - mean) / h, 2);
    return constant_1 * exp(exponent);
}

/// 读取文件内容，将数据存放入数组中
vector<pair<double, double>> read_file(FILE* file) {
    // TODO: 
}

/// 获取核函数自变量取值范围[x_min, x_max]
double get_x_min(vector<pair<double, double>> data) {
    // TODO:
}
double get_x_max(vector<pair<double, double>> data) {
    // TODO:
}

/// 计算样本点与像素点之间欧式距离的平方和
double get_distance_square(vector<pair<double, double>> data) {
    // TODO:
}

// 计算像素点与样本点之间距离的平方
double get_points_distance_square(pair<double, double> p1, pair<double, double> p2) {
    return pow(p1.first - p2.first, 2) + pow(p1.second - p2.second ,2);
}

// 处理结果，将结果放入文件中
void handle_result(char* file, vector<double> result) {
    // TODO:
}

char dataset[100]; //文件路径
FILE* file = NULL;
char result_file[100] = "result.txt";

const double gamma = 1.0 ; //  表示公式中的γ 
const double h = 1.0 ;

int main(int argc, char** argv) 
{
    // 获取参数中的文件路径，并读取文件
    strcpy(dataset, argv[1]);
    file = fopen(dataset, "r");
    if (file == NULL) {
        printf("can not find file named: %s.\n",dataset);
    }
    vector<pair<double, double>> data = read_file(file);
    double n = data.size() ;

    double x_min = get_x_min(data);
    double x_max = get_x_max(data);

    // 二次函数逼近上界 Upper Bound Function
    double a_u = ((x_max-x_min+1)*exp(-x_max) - exp(-x_min)) / pow((x_max-x_min), 2);
    double b_u = ( exp(-x_max) - exp(-x_min)) / ( x_max - x_min ) - a_u*(x_max + x_min);
    double c_u = ( exp(-x_min)*x_max - exp(-x_max)*x_min) / (x_max - x_min) + a_u*x_max*x_min;

    // 二次函数逼近下界 Lower Bound Function
    double t = gamma / n * get_distance_square(data);
    double a_l = (exp(-x_max) + (x_max - 1 - t)*exp(-t)) / pow((x_max -t), 2) ;
    double b_l = -exp(-t) - 2*t*(exp(-x_max)+(x_max-1-t)*exp(-t)) / pow((x_max - t) , 2);
    double c_l = (1+t)*exp(-t) + t*t*(exp(-x_max)+(x_max-1-t)*exp(-t)) / pow((x_max - t), 2);

    // 计算结果 并存储
    vector<double> result; // 结果集合
    vector<pair<double, double>> q;
    // 对于像素集合中的每个像素点（从左到右，从上到下），将其概率密度存放到result中
    double w = 1.0 / (n * h * sqrt(2*M_PI)); 
    for(int i = 0; i< q.size(); i++) {
        double sum = 0;
        for(int j = 0; j < n ;j++) {
            double x = gamma * get_points_distance_square(data[j], q[i]);
            double exp = ((a_u*x*x + b_u*(-x) + c_u) + (a_l*x*x + b_l*(-x) + c_l)) / 2;
            sum += w * exp;
        }
        result.push_back(sum);
    }

    // 存储结果
    handle_result(result_file, result);
}