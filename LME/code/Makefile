CC = gcc
LD = gcc
CFLAGS=-lm -O3
LFLAGS=-lm -O3

all: LogisticEmbed LogisticPred DPTEmbed DPTPred

common.o: LogisticEmbed_common.h LogisticEmbed_common.c
	$(CC) -o common.o LogisticEmbed_common.c -c $(CFLAGS)  

ptable.o: PairHashTable.h PairHashTable.c
	$(CC) -o ptable.o PairHashTable.c -c $(CFLAGS) 

ttable.o: TransitionTable.h TransitionTable.c
	$(CC) -o ttable.o TransitionTable.c -c $(CFLAGS) 

gtable.o: GridTable.h GridTable.c
	$(CC) -o gtable.o GridTable.c -c $(CFLAGS) 

ltable.o: LandmarkTable.h LandmarkTable.c
	$(CC) -o ltable.o LandmarkTable.c -c $(CFLAGS) 

LogisticEmbed.o: LogisticEmbed.h LogisticEmbed.c
	$(CC) -o LogisticEmbed.o LogisticEmbed.c -c $(CFLAGS) 
	
LogisticPred.o: LogisticPred.h LogisticPred.c
	$(CC) -o LogisticPred.o LogisticPred.c -c $(CFLAGS) 

LogisticEmbed: LogisticEmbed.o common.o ptable.o ttable.o gtable.o ltable.o
	$(LD) -o LogisticEmbed LogisticEmbed.o common.o ptable.o ttable.o gtable.o ltable.o $(LFLAGS) 

LogisticPred: LogisticPred.o common.o ptable.o ttable.o gtable.o ltable.o
	$(LD) -o LogisticPred LogisticPred.o common.o ptable.o ttable.o gtable.o ltable.o $(LFLAGS)

DPTEmbed.o: DPTEmbed.h DPTEmbed.c
	$(CC) -o DPTEmbed.o DPTEmbed.c -c $(CFLAGS) 

DPTPred.o: DPTPred.h DPTPred.c
	$(CC) -o DPTPred.o DPTPred.c -c $(CFLAGS) 

DPTEmbed: DPTEmbed.o common.o ptable.o ttable.o gtable.o ltable.o
	$(LD) -o DPTEmbedEmbed DPTEmbed.o common.o ptable.o ttable.o gtable.o ltable.o $(LFLAGS) 

DPTPred: DPTPred.o common.o ptable.o ttable.o gtable.o ltable.o
	$(LD) -o DPTPred DPTPred.o common.o ptable.o ttable.o gtable.o ltable.o $(LFLAGS)

clean:
	rm -f *.o LogisticEmbed LogisticPred DPTEmbed DPTPred
