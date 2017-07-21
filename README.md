
## Installation
The easiest way to install Minimac4 and its dependencies is to use [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget).
```bash
cget install --prefix <install_prefix> statgen/Minimac4
```

Alternatively, you can setup a dev environment cmake directly.
```bash
cd m3vcftools
cget install -f ./requirements.txt                      # Install dependencies locally.
mkdir build && cd build                                 # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                    # Build.
```
