//
// Created by julius on 25-1-22.
//

#ifndef UTILS_H
#define UTILS_H
#include "top.h"

std::vector<Point> readColumns(const std::string& filename, int col1, int col2);

double gaussian(double x, double h);

std::pair<Point, Point> findRanges(const std::vector<Point>& points);
namespace QUAD {
    const double gamma = 1.0 ; //  表示公式中的γ
    const double h = 1.0 ;
    const double e = 0.01; //  ε
    const double t = 0.1; //τ
}

#endif //UTILS_H

















