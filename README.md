# SCALE (Sterge Computer ALgebra Emporium)

A collection of computer algebra packages for Go.

## SFAX (Sterge Fields And eXtensions)
SFAX implements various algebraic structures using Go's generics
* Support for Groups, Rings and Fields
* Implementations for prime fields, rationals, finite extensions, polynomial rings, and fields of rational functions
## SLAB (Sterge Linear AlgeBra)
SLAB implements some linear algebra routines, mostly for $\mathbb{Q}$ and $\mathbb{Z}$.
* Generic matrix multiplication for field elements from SFAX
* Multimodular and multithreaded multiplication for $\mathbb{Z}$
* Cleared denominator multiplication for $\mathbb{Q}$