# Makefile for U(1) project


# Default rule for C compilation
.c.o: ; cc -c $*.c -lm

pureU1: pureU1.c
	cc -o pureU1 pureU1.c -lm
	@echo "made pureU1"

u1dev: u1dev.c
	cc -o u1dev u1dev.c -lm

u1twist: u1twist.o u1utils.o u1update.o measure.o flux.o staples.o
	cc -o u1twist u1twist.o u1utils.o u1update.o measure.o \
	flux.o staples.o -lm
	@echo "Made u1twist"
