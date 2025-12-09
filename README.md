Ortheus
=======

Evolutionary transducer based method for aligning and reconstructing indel
histories (Ensembl Compara modifications)

This repository holds the necessary changes of [Benedict Paten's version](https://github.com/benedictpaten/ortheus)
to run the latest Ensembl Compara pipeline.
Ensembl is **not** the official maintainer of this software.
Pull-requests can still be submitted, but we will only accept them if they
can provide a benefit to Ensembl.

The main new features are:
* new `-T` option to refer to a file that contains the Newick tree, instead
  of passing the Newick tree itself as a command-line argument
* new `-F` option to refer to a file that contains the paths of the Fasta
  files to align, instead of passing these as command-line arguments.

Significant bugs that are now fixed are:
* wrong branch lengths when rerooting the tree
* off-by-one gaps in the alignments
* _AttributeError_ failure in `calculateProbableRootOfGeneTree`

## Branches and tags

There is a single branch (master) where all the development goes.
The version number stated in the source code (0.5.0) is not maintained.

`ensembl_production_XX` tags are used to refer to the version used for the
production of Ensembl version XX. Due to deployment constraints, these tags
may not include the latest changes of the master branch.
Instead, we provide `ensembl_release_candidate_Y` tags, Y starting from 1, for
the "next" version we will deploy in production.

## Installation

### Dependencies

(1) Download and install Pecan

(2) Download and install Semphy, if tree estimation is required

### Installation

Within a suitable environment, ortheus may be installed as follows:
```
git clone --recurse-submodules https://github.com/Ensembl/ortheus.git
cd ortheus
pip install .
```

If you would like to run unit tests, you might consider specifying the
`test` optional dependencies when calling the `pip install` command:
```
pip install .[test]
```

### Testing

ortheus may be tested from within its repository's root directory by calling
`make test`, or by typing `python3 allTests.py` to run the tests directly.

The tests will not verify the installation of Pecan or Semphy currently, but rather just test the core Ortheus algorithms.

See https://github.com/benedictpaten/ for different projects.

For usage instructions see doc/README
