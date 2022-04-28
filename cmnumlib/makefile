CC = gcc
BUILD = build/
TARGET = main
SRCS = src/*.c
HDRS = lib/*.h
LIBS = -lm -lrt -lasound -ljack -lpthread -lportaudio
OPTS = -D_DEBUG_ -w

.SILENT: all clear run

all: ${SRCS} ${HDRS} clear
	${CC} ${SRCS} ${OPTS} -o ${BUILD}${TARGET}.o ${LIBS} > /dev/null

run: all
	./${BUILD}${TARGET}.o 

out:
	./${BUILD}${TARGET}.o > $1

clear:
	cd ${BUILD} ; rm -f *.bin *.hex *.o *.elf
