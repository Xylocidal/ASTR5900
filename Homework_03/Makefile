#
# Flags for debugging:
#

CC := gcc
CFLAGS := -g -std=c99 -fsanitize=address


#LINKER FLAGS

# LFLAGS = -Lgsl/lib -lgsl -lgslcblas -lm
# INCLDIRS = -Igsl/include

#Source files:
OBJ := numerical_integration.o
TARGET := numint.exe

$(TARGET): $(OBJ)
# 	$(CC) $(CFLAGS) $(INCLDIRS) -o $@ $^ $(LFLAGS)
	$(CC) $(CFLAGS) -o $@ $^
TARGETS := $(TARGET)
default: $(TARGET)

all: $(TARGETS)

.PHONY: clean

clean:
	\rm -f *.o *.mod $(TARGETS) 
