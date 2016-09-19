CXXFLAGS = -O3 -Wno-cpp -fmessage-length=0 -fPIC 

OBJS = py_hh.o hh_main.o hh_main_gpu.o

LIBS =

TARGET = py_hh.so

# PYPATH = /home/pavel/opt/python2.7.11
PYPATH = /home/pavel/opt/python27
# PYPATH = /usr/

all: $(TARGET)

py_hh.cpp: py_hh.pyx
	cython --cplus $< -o $@

py_hh.o: py_hh.cpp
	$(CXX) -c $(CXXFLAGS) `$(PYPATH)/bin/python2-config --includes` -I"$(PYPATH)/lib/python2.7/site-packages/numpy/core/include/" $< -o $@

hh_main.o: hh_main.cpp
	$(CXX) -c $(CXXFLAGS) -I/usr/local/cuda-7.5/include/ $< -o $@

hh_main_gpu.o: hh_main_gpu.cu
#	/home/pavel/opt/cuda-7.5/bin/nvcc -c -O3 -arch sm_21 --use_fast_math -Xcompiler -fPIC $< -o $@
	nvcc -c -O3 -arch sm_21 --use_fast_math -Xcompiler -fPIC $< -o $@
	
$(TARGET): $(OBJS)
#	$(CXX) $(CXXFLAGS) -o $@ -shared $(OBJS) -L/home/pavel/opt/cuda-7.5/lib64/ -lcudart 
	$(CXX) $(CXXFLAGS) -o $@ -shared $(OBJS) -lcudart 

clean:
	rm -f $(OBJS) $(TARGET) py_hh.cpp 
