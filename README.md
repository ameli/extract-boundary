# Extract Boundary

This program extracts all __inner__ and __outer__ boundaries of 2D grids and identifies __zero-velocity__ boundary points. Moreover, __inward__ or __outward__ normal vectors on boundary points are calculated. You can use 2D grids with __any type of cell__ and various file extensions such as `*.vtk` and `*.vtu`.

## User Guide
Please refer to [Documentation page](https://ameli.github.io/extract-boundary) or [doc/UserGuide.pdf](https://github.com/ameli/extract-boundary/raw/main/doc/UserGuide.pdf).

## Quick installation
You may install by building source or getting the ready package from Ubuntu repos.

### 1. Building from source

Just download the [installer](https://raw.github.com/ameli/extract-boundary/main/installer) and do the following:

    $ sudo chmod +x installer
    $ sudo ./installer

### 2. Ubuntu Package

For Ubuntu users, the package has been added to the Ubuntu repository. Type in terminal:

    $ sudo add-apt-repository ppa:ameli/extractboundary
    $ sudo apt-get update
    $ sudo apt-get install extractboundary

## Sample Usage

You may use the code either in command line or as a ParaView plugin.

### 1. Command Line

* Specifying input and output paths within files:

    `$ extractboundary /InputPath/InputData.vtk /OutputPath/OutputData.vtk`

* No output file is specifies. Default name will be used for output file:

    `$ extractboundary InputData.vtk`

* Using options: 
    * outward normals (with -n)
    * search all inner and outer boundaries (with -s)
    * identify zero-velocity points (with -z)

    `$ extractboundary InputData.vtk OutputData.vtk -n outward -s all -z`

### 2. ParaView Plugin

In ParaView, from _tools_ menu open _Manage Plugins_. Then load [`bin/libExtractBoundary.so`](https://github.com/ameli/extract-boundary/raw/main/bin/libExtractBoundary.so) file. Next, in _Filters_ menu go to _Extensions_ and apply _Extract Boundary_ filter to your pipeline.

You may refer to [UserGuide](https://github.com/ameli/extract-boundary/raw/main/doc/UserGuide.pdf) or [wiki](https://github.com/ameli/extract-boundary/wiki/Extract-Boundary) for picture illustrations.

## License

Its free and open source under GNU/zlib license. Please see [LICENSE.txt](https://raw.github.com/ameli/extract-boundary/main/LICENSE.txt) for terms.

## Author

Siavash Ameli  
[Shadden Research Group](http://shaddenlab.berkeley.edu/)  
University of California, Berkeley
