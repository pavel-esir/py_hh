CXXFLAGS =	-O3 -Wall -fmessage-length=0 -fPIC

OBJS =		py_hh.o hh_main.o

LIBS =

TARGET =	py_hh.so

PYPATH = /home/pavel/opt/python2.7.11/bin
#PYPATH = /home/pavel/opt/python27/bin

py_hh.cpp: py_hh.pyx
	cython --cplus $< -o $@

py_hh.o: py_hh.cpp
	$(CXX) -c -fPIC -O3 -Wno-cpp `$(PYPATH)/python2-config --include` -I"/home/pavel/opt/python2.7.11/lib/python2.7/site-packages/numpy/core/include/" $< -o $@
#	$(CXX) -c -g -fPIC -O3 -Wno-cpp `$(PYPATH)/python2-config --include` -I"/home/pavel/opt/python27/lib/python2.7/site-packages/numpy/core/include/" cython_cover.cpp -o $@

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) -shared $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) cython_cover.cpp
