The source code for building the `mkhssp` and `hsspconv`
programs is bundled in the `hssp` project. The mkhssp executable
creates stockholm files with hssp annotations in them. The hsspconv
executable converts stockholm to the original hssp format.

# Development

The provided Dockerfile sets up a development environment. Build the docker
image using the command `docker build -t hssp .` and run the image in a
container, with a local source copy and data files mounted, with the command
`docker run -v /home/jon/projects/hssp:/app -v /mnt/extra:/srv/data -it hssp`.

# Download and installation instructions

Source code archives are available [here][2].

## Compiling hssp programs

### Pre-requisites

Compiler:
* Must support at least the c++ 11 standard.

System libraries:

* libzeep version >= 3.0 (for mkhssp --fetch-dbrefs only)
* libboost version >= 1.48
* libz
* libbz2
* autoconf
* automake
* autotools-dev

### Instructions

Download and uncompress the hssp [source code archive][2] (version >= 2.2.6):

    wget https://github.com/cmbi/hssp/archive/hssp-?.?.?.tar.gz
    tar -zxvf hssp-?.?.?.tar.gz
    cd hssp-?.?.?

Configure and build the hssp executables:

    ./autogen.sh
    ./configure
    make

To build only one executable of the hssp project, e.g. `mkhssp`, type:

    make mkhssp

To test the `mkhssp` executable type:

    ./mkhssp

To add the executables to /usr/local/bin type:

    sudo make install

# Citing hssp

The [reference][3] for the new versions of hssp and other protein structure
bioinformatics [facilities][4] is:

```
A series of PDB-related databanks for everyday needs
Wouter G. Touw, Coos Baakman, Jon Black, Tim A. H. te Beek,
 E. Krieger, Robbie P. Joosten and Gert Vriend.
Nucl. Acids Res. (2015) 43, D364-D368
```

The original reference for DSSP is:

```
Dictionary of protein secondary structure: pattern recognition of
 hydrogen-bonded and geometrical features.
Kabsch W and Sander C, Biopolymers (1983) 22, 2577-2637.
```

The original reference for HSSP is:

```
Database of homology-derived protein structures and the structural
 meaning of sequence alignment.
Sander C and Schneider R, Proteins (1991) 9, 56-68.
```

# Contact

In 2013, maintenance of hssp has been taken over from Maarten Hekkelman by
Coos Baakman. If you want to provide feedback,
either send an e-mail to hssp.cmbi@radboudumc.nl or have a look at
[existing issues][5] (if necessary, create a new issue).


[1]: http://swift.cmbi.umcn.nl/gv/hssp/
[2]: https://github.com/cmbi/hssp/releases
[3]: http://dx.doi.org/10.1093/nar/gku1028
[4]: http://swift.cmbi.umcn.nl/gv/facilities/
[5]: https://github.com/cmbi/hssp/issues
