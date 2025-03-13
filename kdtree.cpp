
#include "kdtree.h"



// 欧几里得距离平方计算
double distanceSquared(const Point& p1, const Point& p2) {
    return std::pow(p1.first - p2.first, 2) + std::pow(p1.second - p2.second, 2);
}

// 计算点到矩形的最小距离
double minDistanceToRectangle(const Point& point, const Rectangle& rect) {
    // 投影点到矩形边界的最近点
    double px = std::max(rect.min.first, std::min(point.first, rect.max.first));
    double py = std::max(rect.min.second, std::min(point.second, rect.max.second));

    // 返回该点到投影点的距离
    return distanceSquared(point, Point(px, py));
}

// 计算点到矩形的最大距离
double maxDistanceToRectangle(const Point& point, const Rectangle& rect) {
    // 四个角点
    Point corners[4] = {
        Point(rect.min.first, rect.min.second),
        Point(rect.min.first, rect.max.second),
        Point(rect.max.first, rect.min.second),
        Point(rect.max.first, rect.max.second)
    };

    // 初始化最大距离
    double maxDist = 0;

    // 遍历所有角点，寻找最大距离
    for (const auto& corner : corners) {
        double dist = distanceSquared(point, corner);
        if (dist > maxDist) {
            maxDist = dist;
        }
    }

    return maxDist;
}

// 函数：计算包围一组点的最小矩形
std::pair<std::pair<double, double>, std::pair<double, double>>
calculateBoundingBox(const std::vector<Point>& points) {
    if (points.empty()) {
        throw std::invalid_argument("Point set cannot be empty.");
    }

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& point : points) {
        minX = std::min(minX, point.first);
        maxX = std::max(maxX, point.first);
        minY = std::min(minY, point.second);
        maxY = std::max(maxY, point.second);
    }

    // 返回对角线上的两个点来表示矩形
    return {{minX, minY}, {maxX, maxY}};
}

// 比较函数，用于根据x或y排序
bool compareX(const Point& a, const Point& b) { return a.first < b.first; }
bool compareY(const Point& a, const Point& b) { return a.second < b.second; }
// 函数：将点集按照指定维度划分成两个点集
void splitPointsByDimension(std::vector<Point>& points, bool splitByX,
                            std::vector<Point>& leftSet, std::vector<Point>& rightSet) {
    if (points.empty()) return;

    // 根据给定维度排序点集
    if (splitByX) {
        std::sort(points.begin(), points.end(), compareX);
    } else {
        std::sort(points.begin(), points.end(), compareY);
    }

    size_t midIndex = points.size() / 2;

    // 如果点集大小为奇数，midIndex指向中点；如果为偶数，它指向上半部分的第一个点
    if (points.size() % 2 == 0) {
        // 对于偶数个元素，我们使用nth_element来找到中位数的位置
        std::nth_element(points.begin(), points.begin() + midIndex - 1, points.end(), splitByX ? compareX : compareY);
        std::nth_element(points.begin() + midIndex, points.end() - 1, points.end(), splitByX ? compareX : compareY);
    } else {
        std::nth_element(points.begin(), points.begin() + midIndex, points.end(), splitByX ? compareX : compareY);
    }

    // 分割点集
    leftSet.assign(points.begin(), points.begin() + midIndex);
    rightSet.assign(points.begin() + midIndex, points.end());
}

std::unique_ptr<KdTreeNode> buildTree(std::vector<Point>& points, bool splitX) {
    if (points.size() < 10) return nullptr;

    // 计算点集边界矩形
    auto [p1, p2] = calculateBoundingBox(points); // Assuming calculateBoundingBox returns a pair of Points
    Rectangle root_r(p1, p2);

    auto root = std::make_unique<KdTreeNode>(root_r);
    root->Rcount = points.size();
    root->points = points;
    // 依次按x和y划分点集为两部分
    std::vector<Point> leftSet, rightSet;
    splitPointsByDimension(points, splitX, leftSet, rightSet);

    root->left = buildTree(leftSet, !splitX);
    root->right = buildTree(rightSet, !splitX);

    return root;
}












