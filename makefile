

ode: ode.cpp ode.h
	g++ -o ode ode.cpp

femrunner: fem.cpp PyMem.o ode.o timer.h
	g++ fem.cpp PyMem.o ode.o -o femrunner -lrt

PyMem.o: PyMem.cpp
	g++ -o PyMem.o PyMem.cpp -c

ode.o: ode.cpp ode.h
	g++ -o ode.o ode.cpp -c

