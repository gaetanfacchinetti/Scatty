CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c11 -fPIC
LDFLAGS = -shared
OUTPUT = libscatty.so
EXEC = scatty

SRC = integration.c cross_section.c tools.c rates.c stubs.c cmplx.c cgamma.c clog.c
MAIN_SRC = scatty.c
OBJ = $(SRC:.c=.o)
MAIN_OBJ = $(MAIN_SRC:.c=.o)

all: $(OUTPUT) $(EXEC)

$(OUTPUT): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $(OBJ) -lm

$(EXEC): $(OBJ) $(MAIN_OBJ)
	$(CC) -o $@ $^ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(MAIN_OBJ)

cleanlib:
	rm -f $(OUTPUT)

cleanexe:
	rm -f $(EXEC)

distclean: clean cleanlib cleanexe