# hemoMIPs



## Pre-requirements

### Conda

### hemoMIPs Conda environment

### Shed Skin

is an experimental compiler, that can translate pure, but implicitly statically typed Python (2.4-2.6) programs into optimized C++. To fasten the mapping process one one of our python scripts have to be translated to C++ with Shed Skin.

First we need an environment with python v2.6 and the requirements for Shed Skin. Therefore we created the environment file `envs/shedskin.yml`. Be sure that you are in your root hemoMIPs pipeline folder.
```
conda env create -f envs/shedskin.yml -n shedskin
source activate shedskin
```
Then download and install Shed Skin v0.9.4 into the bin directory of the hemoMIPs pipeline.

```bash
# Download Shed Skin 0.9.4
wget https://github.com/shedskin/shedskin/releases/download/v0.9.4/shedskin-0.9.4.tgz
# Create bin folder
mkdir -p bin
# Extract Shed Skin and remove file
tar -xzf shedskin-0.9.4.tgz -C bin
rm shedskin-0.9.4.tgz
# Install Shed Skin
cd bin/shedskin-0.9.4
python setup.py install
```
Now we can test the shedskin intalation:
```bash
shedskin test
make
```
The result should look like
```
*** SHED SKIN Python-to-C++ Compiler 0.9.4 ***
Copyright 2005-2011 Mark Dufour; License GNU GPL version 3 (See LICENSE)

[analyzing types..]
********************************100%
[generating c++ code..]
[elapsed time: 1.29 seconds]
```
If the installation or test fails please have a look a the [Shed Skin Dokumentation](https://shedskin.readthedocs.io/en/latest/). Maybe you need the Boehm garbage collector.

#### Compiling XYZ scriipt

Now we need to compile the XYZ script using Shed Skin:

```bash
# Go back to the root folder
cd ../../
shedskin sc
```
