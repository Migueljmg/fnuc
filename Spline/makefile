define CC
	g++
endef
CFLAGS = -Wall
LDFLAGS = -lm
MDFLAGS := `root-config --cflags --glibs`

# exercicio 3 
OBJS = mainInterpolator.o EqSolver.o Vec.o cFCgraphics.o DataInterpolator.o 
SRCS = mainInterpolator.C EqSolver.C Vec.C cFCgraphics.C DataInterpolator.C 
HEAS = mainInterpolator.h EqSolver.h Vec.h cFCgraphics.h DataInterpolator.h

#actions
%.o: %.C $(HEAS)
	$(CC) -c -o $@ $< $(OBJS) $(CFLAGS)
mainInterpolator.exe: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ 
	@echo compiling $?
mainInterpolator: $(OBJS)
	g++ -o $@ $(SRCS) $(CFLAGS) -I /SOFT/root-53406 $(MDFLAGS)

clean:
	rm -f *.o libTrab2.a

Interpolator:libTrab2.a
	g++ -o mainInterpolator.exe mainInterpolator.o DataInterpolator.o Vec.o EqSolver.o -I $(MDFLAGS)
mainInterpolator.o:
	g++ -c mainInterpolator.C  `root-config --cflags --glibs`
DataInterpolator.o:
	g++ -c DataInterpolator.C `root-config --cflags --glibs`
Vec.o:
	g++ -c Vec.C
EqSolver.o: 
	g++ -c EqSolver.C
libTrab2.a:mainInterpolator.o DataInterpolator.o Vec.o EqSolver.o
	ar ruv libTrab2
	ranlib libTrab2


