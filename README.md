
## Installation
The easiest way to install m3vcftools is to setup a dev environment cmake directly.
```bash
cd m3vcftools
cget install -f ./requirements.txt                      # Install dependencies locally.
mkdir build && cd build                                 # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                    # Build.
```
