CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c11 -fPIC
LDFLAGS = -shared
OUTPUT = libscatty.so

SRC = integration.c cross_section.c tools.c rates.c stubs.c cmplx.c cgamma.c clog.c
OBJ = $(SRC:.c=.o)

all: $(OUTPUT)

$(OUTPUT): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $(OBJ) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)

cleanlib: 
	rm -f $(OUTPUT)