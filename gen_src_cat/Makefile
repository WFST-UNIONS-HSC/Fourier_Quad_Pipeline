CXX = mpicxx

CPPFLAGS = -Iinclude
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Wpedantic -MMD -MP
LDFLAGS =
LDLIBS =

SRCS = main.cpp \
       src/process_extcat/process_extcat.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
TARGET = process_extcat

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

test: $(TARGET)
	python3 tests/test_process_extcat.py ./$(TARGET)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)

.PHONY: all test clean
