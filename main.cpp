#include "utils.h"
#include "kdtree.h"
#include "top.h"
using namespace std;

#define SCALA_X 320
#define SCALA_Y 240

// 计算像素点与样本点之间距离的平方
double get_points_distance_square(pair<double, double> p1, pair<double, double> p2) {
    return (p1.first - p2.first) * (p1.first - p2.first) + (p1.second - p2.second) * (p1.second - p2.second);
}

/// 计算样本点与像素点之间欧式距离的平方和
double get_distance_square(vector<pair<double, double>> data, Point p) {
    double s = 0.0;
    for (auto point : data) {
        double i = get_points_distance_square(point, p);
        if (std::isnan(i)) {
            continue;;
        }
        //cout << i << endl;
        s += i;

    }
    //cout << s << endl;
    return s;
}


// 处理结果，将结果放入文件中
void handle_result(char* file, vector<double> result) {
    // TODO:
}

double compute_u_QUAD(double a_u, double b_u, double c_u, double x) {
  	return a_u * x * x + b_u * x + c_u;
}
double compute_l_QUAD(double a_l, double b_l, double c_l, double x) {
  return a_l * x * x + b_l * x + c_l;
}
double compute_gaussian_QUAD(KdTreeNode* node, double x, double y) {
  double sum = 0;
  for (auto point : node->points) {
    double dis_square = get_points_distance_square(point, {x, y});
	sum += gaussian(dis_square, QUAD::h);
  }
  return sum / node->Rcount;
}

void get_params(vector<Point>& data,
                double& x_min, double& x_max, double& y_max, double& y_min,
                double& a_u, double& b_u, double& c_u,
                double& a_l, double& b_l, double& c_l,
                double& x_per, double& y_per) {
  int n = data.size();
  x_min = findRanges(data).first.first;
  x_max = findRanges(data).first.second;
  y_min = findRanges(data).second.first;
  y_max = findRanges(data).second.second;
  double x_range = x_max - x_min;
  double y_range = y_max - y_min;
	x_per = x_range / SCALA_X;
  y_per = y_range / SCALA_Y;

    // 二次函数逼近上界 Upper Bound Function
  a_u = ((x_max-x_min+1)*exp(-x_max) - exp(-x_min)) / (x_range * x_range);
  b_u = ( exp(-x_max) - exp(-x_min)) / ( x_max - x_min ) - a_u*(x_max + x_min);
  c_u = ( exp(-x_min)*x_max - exp(-x_max)*x_min) / (x_max - x_min) + a_u*x_max*x_min;

    // 二次函数逼近下界 Lower Bound Function
  double t = QUAD::gamma * get_distance_square(data, data[4]) / n;
    //cout << t << endl;
  a_l = (exp(-x_max) + (x_max - 1 - t)*exp(-t)) / pow((x_max -t), 2) ;
  b_l = -exp(-t) - 2*t*(exp(-x_max)+(x_max-1-t)*exp(-t)) / pow((x_max - t) , 2);
  c_l = (1+t)*exp(-t) + t*t*(exp(-x_max)+(x_max-1-t)*exp(-t)) / pow((x_max - t), 2);
  cout << a_l << " " << b_l << " " << c_l << endl;
  cout << a_u << " " << b_u << " " << c_u << endl;
}

const std::string name = "HT_Sensor_dataset.dat"; //文件路径

class compare {
	double x,y;
    double a_u, b_u, c_u;
    double a_l, b_l, c_l;
public:
  	compare(double x, double y, double a_u, double b_u, double c_u, double a_l, double b_l, double c_l) {
          this->x = x;
          this->y = y;
          this->a_u = a_u;
          this->b_u = b_u;
          this->c_u = c_u;
          this->a_l = a_l;
          this->b_l = b_l;
          this->c_l = c_l;
  	}
	bool operator() (const KdTreeNode* a, const KdTreeNode* b) {
    	Point q({x,y});
          double min_a_u = QUAD::gamma * minDistanceToRectangle(q, a->boundingBox);
          double min_b_u = QUAD::gamma * minDistanceToRectangle(q, b->boundingBox);
          double max_a_l = QUAD::gamma * maxDistanceToRectangle(q, a->boundingBox);
          double max_b_l = QUAD::gamma * maxDistanceToRectangle(q, b->boundingBox);
          double u_a = 1 * a->Rcount * ( min_a_u * min_a_u * a_u + min_a_u * b_u + c_u); //w = 1
          double l_a = 1 * a->Rcount * ( max_a_l * max_a_l * a_l + max_a_l * b_l + c_l);
          double u_b = 1 * b->Rcount * ( min_b_u * min_b_u * a_u + min_b_u * b_u + c_u);
          double l_b = 1 * b->Rcount * ( max_b_l * max_b_l * b_l + max_b_l * b_l + c_l);
          return (u_a - l_a) < (u_b - l_b);
	}
};


int main(int argc, char** argv)
{

    vector<Point> data = readColumns(name, 11, 12);
    double n = data.size() ;
    double x_min,x_max,y_max,y_min;
    double a_u, b_u, c_u;
    double a_l, b_l, c_l;
    double x_per,  y_per;
    auto root = buildTree(data);
    double result[SCALA_X][SCALA_Y];
    for(int i = 0; i < SCALA_X; i++) {
      for(int j = 0; j < SCALA_Y; j++) {
          //cout << "[" << i << "," << j << "]" << endl;
        double x = x_min + i * x_per;
        double y =  y_min + j * y_per;
        priority_queue<KdTreeNode*, vector<KdTreeNode*>, compare> Pri_queue(compare(x,y,a_u,b_u,c_u,a_l,b_l,c_l));
		Pri_queue.push(root.get());
        get_params(root->points, x_min, x_max, y_max, y_min, a_u, b_u, c_u, a_l, b_l, c_l, x_per, y_per);
        double x_l = QUAD::gamma * maxDistanceToRectangle({x, y}, root->boundingBox);
        double l = 1 * root->Rcount * compute_l_QUAD(a_l,b_l,c_l,x_l);
        double x_u = QUAD::gamma * minDistanceToRectangle({x, y}, root->boundingBox);
        double u = 1 * root->Rcount * compute_u_QUAD(a_u,b_u,c_u,x_u);
        while (u > ( 1 + QUAD::e) * l) {
			KdTreeNode* node = Pri_queue.top();
            Pri_queue.pop();
            l -= 1 * node->Rcount * compute_l_QUAD(a_l,b_l,c_l,x_l);
            u -= 1 * node->Rcount * compute_u_QUAD(a_u,b_u,c_u,x_u);
            if (node->left) {
              Pri_queue.push(node->left.get());
              get_params(node->left->points, x_min, x_max, y_max, y_min, a_u, b_u, c_u, a_l, b_l, c_l, x_per, y_per);
              x_l = QUAD::gamma * maxDistanceToRectangle({x, y}, node->left->boundingBox);
              l += 1 * node->left->Rcount * compute_l_QUAD(a_l,b_l,c_l,x_l);
              x_u = QUAD::gamma * minDistanceToRectangle({x, y}, node->left->boundingBox);
              u += 1 * node->left->Rcount * compute_u_QUAD(a_u,b_u,c_u,x_u);
            }
            if (node->right) {
              Pri_queue.push(node->right.get());
              get_params(node->right->points, x_min, x_max, y_max, y_min, a_u, b_u, c_u, a_l, b_l, c_l, x_per, y_per);
              x_l = QUAD::gamma * maxDistanceToRectangle({x, y}, node->right->boundingBox);
              l += 1 * node->right->Rcount * compute_l_QUAD(a_l,b_l,c_l,x_l);
              x_u = QUAD::gamma * minDistanceToRectangle({x, y}, node->right->boundingBox);
              u += 1 * node->right->Rcount * compute_u_QUAD(a_u,b_u,c_u,x_u);
            }
            if (!node->left && !node->right) {
				l += compute_gaussian_QUAD(node, x, y);
                u += compute_gaussian_QUAD(node, x, y);
            }
        }
        result[i][j] = (u + l) / 2;
          cout << result[i][j] << endl;
      }
    }
}