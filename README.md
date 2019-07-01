The source code for building the `mkdssp`, `mkhssp`, `hsspconv`, and
`hsspsoap` programs is bundled in the `xssp` project. The DSSP executable is
`mkdssp`.

# Development

The provided Dockerfile sets up a development environment. Build the docker
image using the command `docker build -t xssp .` and run the image in a
container, with a local source copy and data files mounted, with the command
`docker run -v /home/jon/projects/xssp:/app -v /mnt/extra:/srv/data -it xssp`.

# Download and installation instructions

Pre-compiled *old* versions of DSSP are available from the
[old repository][1]. New source code archives are available [here][2].

## Compiling xssp programs

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

Download and uncompress the xssp [source code archive][2] (version >= 2.2.6):

    wget https://github.com/cmbi/xssp/archive/xssp-?.?.?.tar.gz
    tar -zxvf xssp-?.?.?.tar.gz
    cd xssp-?.?.?.tar.gz

Configure and build the xssp executables:

    ./autogen.sh
    ./configure
    make

To build only one executable of the xssp project, e.g. `mkdssp`, type:

    make mkdssp

To test the `mkdssp` executable type:

    ./mkdssp

To add the executables to /usr/local/bin type:

    sudo make install

# Citing xssp

The [reference][3] for the new versions of xssp and other protein structure
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

In 2013, maintenance of xssp has been taken over from Maarten Hekkelman by
Coos Baakman, Jon Black, and Wouter Touw. If you want to provide feedback,
either send an e-mail to xssp.cmbi@radboudumc.nl or have a look at
[existing issues][5] (if necessary, create a new issue).


[1]: http://swift.cmbi.umcn.nl/gv/dssp/DSSP_5.html
[2]: https://github.com/cmbi/xssp/releases
[3]: http://dx.doi.org/10.1093/nar/gku1028
[4]: http://swift.cmbi.umcn.nl/gv/facilities/
[5]: https://github.com/cmbi/xssp/issues
