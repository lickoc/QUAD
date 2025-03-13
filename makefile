# 编译器
CXX = g++

# 编译选项，比如C++标准版本
CXXFLAGS = -std=c++14 -Wall

# 目标可执行文件名称
TARGET = main

# 源文件和对象文件
SRCS = main.cpp kdtree.cpp utils.cpp
OBJS = $(SRCS:.cpp=.o)

# 默认目标
all: $(TARGET)

# 链接对象文件生成可执行文件
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# 规则：从 .cpp 文件生成 .o 文件
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清理生成的文件
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
