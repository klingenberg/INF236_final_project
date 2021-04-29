# the compiler: gcc for C program, define as g++ for C++
CC = clang

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -Xpreprocessor -fopenmp -lomp -lm -O3

# the build target executable:
TARGET = driver

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c

clean:
	$(RM) $(TARGET)
