# multiz
DNA multiple sequence aligner, official version from Penn State's Miller Lab

This is a nearly functionally equivalent copy of the tarball and documentation
for multiz that was posted at the Miller Lab website January 21, 2009,
`http://www.bx.psu.edu/miller_lab`.

The motivation for this repository is to provide an equivalent of that tarball
that can be installed by modern package managers.

The only functional change that has been made is that this version uses lastz
to generate the underlying pariwise alignments. The 2009 release used blastz
instead. Blastz has now been deprecated for nearly a decade, and lastz is a
drop-in replacement for it.

The only other changes are the addition of this README.md, the MIT license,
and several syntax changes to prevent compiler warnings.
