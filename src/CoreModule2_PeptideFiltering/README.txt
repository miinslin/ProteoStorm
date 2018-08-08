Compile CoreModule2_PeptideFiltering.exe from source
==========
### Software requirements:
1. ```CMake``` (minimum version 3.6)
2. ```Boost C++ Libraries``` (version 1.60.0)

### Guidelines for installing required software:

**Linux **
sudo apt-get install build-essential
* ```CMake```
	```sh
	$ wget https://cmake.org/files/v3.10/cmake-3.10.2.tar.gz
	$ tar xzf cmake-3.10.2.tar.gz
	$ cd cmake-3.10.2
	$ ./bootstrap
	$ make
	$ sudo make install
	```
	
* ```Boost C++ Libraries```
	```sh
	$ wget http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
	$ tar xzf boost_1_60_0.tar.gz
	$ cd boost_1_60_0
	$ ./bootstrap.sh -prefix=boost/
	$ ./bjam
	$ ./bjam install
	```

**Windows**
* ```CMake``` 
	[version 3.10.2](https://cmake.org/files/v3.10/cmake-3.10.2-win64-x64.msi)

* ```Boost C++ Libraries``` 
	[version 1.60.0](http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz)

	```sh
	$ tar xzf boost_1_60_0.tar.gz
	$ cd path/to/boost_1_60_0
	$ ./bootstrap.sh -prefix=boost/
	$ ./bjam
	$ ./bjam install
	```
	Copy the three new directories (i.e., bin, include, and lib) to ~/cygwin/usr/local/


### Create CoreModule2_PeptideFiltering.exe:
* Extract files from CoreModule2_PeptideFiltering.zip
* In ./CoreModule2_PeptideFiltering/CMakeLists.txt, specify:
	- BOOST_ROOT
	- BOOST_INCLUDEDIR
	- BOOST_LIBRARYDIR
* Create executable
	```sh
	$ cd CoreModule2_PeptideFiltering
	$ mkdir build
	$ cd build
	$ cmake ..
	$ make
	```