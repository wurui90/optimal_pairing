SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp = .o)
CC = g++
INCLUDES = 
LIBS = 
CCFLAGS = -g -Wall -std=c++0x
TARGET = optimal_pairing

$(TARGET) : $(TARGET).o
	$(CC) $^ -o $@ $(INCLUDES) $(LIBS)

%.o : %.cpp
	$(CC) -c $< $(CCFLAGS) -o $@

clean:
	rm *.o
	rm $(TARGET)

.PHONY:clean

