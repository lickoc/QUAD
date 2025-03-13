#include "utils.h"

// Function to read specific columns from each line of a space-separated text file, handle missing data represented by '.', and convert them to double
std::vector<std::pair<double, double>> readColumns(const std::string& filename, int col1, int col2) {
    std::ifstream file(filename);
    std::vector<std::pair<double, double>> columns;
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            std::string cell;
            int currentCol = 0;
            double colData1 = NAN, colData2 = NAN; // Initialize with NaN for missing values

            // Read the line into cells separated by spaces
            while (ss >> cell) {
                ++currentCol;
                if (currentCol == col1 || currentCol == col2) {
                    if (cell == "Temp." || cell == "Humidity") {
                        continue;
                    }
                    if (cell != ".") { // Check for missing data representation
                        try {
                            double value = std::stod(cell);
                            if (currentCol == col1) colData1 = value;
                            else if (currentCol == col2) colData2 = value;
                        } catch (const std::invalid_argument& e) {
                            std::cerr << "Invalid argument for conversion: " << cell << " at column " << currentCol << std::endl;
                            // If you want to set a default value instead of NaN, do it here.
                        }
                    } else {
                        // Handle missing data as needed, currently setting to NaN
                        if (currentCol == col1) colData1 = NAN;
                        else if (currentCol == col2) colData2 = NAN;
                    }
                }
                // Once both target columns are read, break the inner loop to avoid unnecessary reading.
                if (currentCol >= col2) break;
            }

            // Add the pair only if both columns have valid or NaN data
            columns.push_back(std::make_pair(colData1, colData2));
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return columns;
}

/// 高斯核密度函数 https://zhuanlan.zhihu.com/p/643143096
/// 参数x ： 样本数据点
/// 参数mean ： 所有像素点
/// 参数h ： 带宽
double gaussian(double x, double h) {
    double constant_1 =  1.0 / (h * sqrt(2*M_PI));
    double exponent = -0.5 * x * x / (h * h);
    return constant_1 * exp(exponent);
}

std::pair<std::pair<double, double>, std::pair<double, double>> findRanges(const std::vector<std::pair<double, double>>& points) {
    if (points.empty()) {
        throw std::invalid_argument("Point set is empty.");
    }

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& point : points) {
        double x = point.first;
        double y = point.second;

        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }

    return {{minX, maxX}, {minY, maxY}};
}

/*
int main() {
    const std::string filename = "HT_Sensor_dataset.dat"; // Replace with your actual file path.
    int col1 = 11;
    int col2 = 12;

    std::vector<std::pair<double, double>> columns = readColumns(filename, col1, col2);

    for (size_t i = 0; i < columns.size(); ++i) {
        std::cout << "Line " << i + 1 << ": Column " << col1
                  << " - " << (std::isnan(columns[i].first) ? "Missing" : std::to_string(columns[i].first))
                  << ", Column " << col2
                  << " - " << (std::isnan(columns[i].second) ? "Missing" : std::to_string(columns[i].second))
                  << std::endl;
    }

    return 0;
}
*/