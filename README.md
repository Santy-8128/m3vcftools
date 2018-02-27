
## Prerequisites

Automatic installation of m3vcftools requires [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget) and cmake v3.2. These prerequisites can be installed as follows:

Ubuntu 16.04
```
sudo apt-get install cmake python-pip python-dev
pip install cget
```
Ubuntu 14.04
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:george-edison55/cmake-3.x
sudo apt-get update
sudo apt-get install cmake python-pip python-dev
pip install cget
```
MacOS
```
brew install cmake
sudo easy-install pip
pip install --user cget --ignore-installed six
```


## Installation
The easiest way to install m3vcftools is to git clone followed by using cget and cmake directly.
```bash
git clone https://github.com/Santy-8128/m3vcftools
cd m3vcftools
cget install -f ./requirements.txt                      # Install dependencies locally.
mkdir build && cd build                                 # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                    # Build.
```
