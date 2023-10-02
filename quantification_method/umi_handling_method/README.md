# To Build this executable
Create a build directory inside this current directory and then run cmake, using the commands below.
```
mkdir build
cd build
cmake ../ -DNO_IPO=T
make
```
- Note `DNO_IPO` is needed for CentOS and REHL 
