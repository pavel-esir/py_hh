CXXFLAGS = -O3 -Wall -fmessage-length=0 -fPIC

OBJS = py_hh.o hh_main.o

LIBS =

TARGET = py_hh.so

#PYPATH = /home/pavel/opt/python2.7.11
#PYPATH = /home/pavel/opt/python27
PYPATH = /usr/

all: $(TARGET)

py_hh.cpp: py_hh.pyx
	cython --cplus $< -o $@

py_hh.o: py_hh.cpp
	$(CXX) -c $(CXXFLAGS) `$(PYPATH)/bin/python2-config --includes` -I"$(PYPATH)/lib/python2.7/site-packages/numpy/core/include/" $< -o $@

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) -shared $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET) py_hh.cpp
