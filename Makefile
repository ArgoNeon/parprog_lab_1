CC = gcc
CCFLAGS = -Wall

SRC = main.c
OBJ = $(SRC:.c = .o)

TARGET = build/task

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CCFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CCFLAGS) $^ -o $@

clean:
	rm -f $(TARGET)
