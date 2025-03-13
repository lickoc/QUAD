//
// Created by julius on 25-1-22.
//

#ifndef KDTREE_H
#define KDTREE_H
#include "top.h"
// 矩形类
struct Rectangle {
    Point min; // 左下角
    Point max; // 右上角
    //构造函数
    Rectangle(Point min, Point max) : min(min), max(max) {}
};

struct KdTreeNode {
    Rectangle boundingBox;
    int Rcount;
    std::unique_ptr<KdTreeNode> left;
    std::unique_ptr<KdTreeNode> right;
    std::vector<Point> points;
    KdTreeNode(const Rectangle& boundingBox): boundingBox(boundingBox), Rcount(0) {}
};

std::unique_ptr<KdTreeNode> buildTree(std::vector<Point>& points, bool splitX = true);
double minDistanceToRectangle(const Point& point, const Rectangle& rect);
double maxDistanceToRectangle(const Point& point, const Rectangle& rect);

#endif //KDTREE_H
