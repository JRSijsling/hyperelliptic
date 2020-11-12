Description
--

This repository contains Magma code for reconstruction and isomorphisms of hyperelliptic curves.

Prerequisites
--

You need to have Magma installed on your machine to use this code.

Installation
--

You can enable the functionality of this package in Magma by attaching the `hyperelliptic/magma/spec` file with `AttachSpec`. To make this independent of the directory in which you find yourself, and to active this on startup by default, you may want to indicate the relative path in your `~/.magmarc` file, by adding a line like
```
AttachSpec("~/Programs/hyperelliptic/magma/spec");
```

Usage
--

Example that tests the routines in this package can be found in the directory `examples`.

Verbose comments are enabled by
```
SetVerbose("G2Twists", n);
```
or
```
SetVerbose("G3Twists", n);
```
where `n` is either `1` or `2`. A higher value gives more comments.

Citing this code
--

Please cite the following preprint if this code has been helpful in your research:

Reynald Lercier and Christophe Ritzenthaler,
*Hyperelliptic curves and their invariants: Geometric, arithmetic and algorithmic aspects*,
[Journal of Algebra 372 (2012) 595–636](http://dx.doi.org/10.1016/j.jalgebra.2012.07.054)
