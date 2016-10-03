CFLAGS = -O3 -march=native -lm -Wall -std=c11
.PHONY: clean test
all: strassen test
strassen: strassen.c strassen.h
	$(CC) $(CFLAGS) $< -o $@
test: strassen
	time ./$< 1000 < test/1000.in > result
release: strassen.c strassen.h test/1000.in Makefile README.md
	zip -9 strassen.zip $^ 
clean:
	rm -rf strassen result
