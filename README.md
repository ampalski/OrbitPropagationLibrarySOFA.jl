# OrbitPropagationLibrarySOFA

[![Build Status](https://github.com/ampalski/OrbitPropagationLibrarySOFA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ampalski/OrbitPropagationLibrarySOFA.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

A pure Julia implementation of the IAU timing and coordinate tranformation functions, as published in the SOFA routines with supplemental material derived from Vallado. 

This package is at minimum viable product, and is intended to be used as the timing and coordinate backbone to an Earth orbit propagation library and not explicitly for astronomy the way SOFA is. As such, some function signatures will be different. As more functionality is added (see below), some function signatures will change.

In order to maintain compatibility with the USSF standards, initial implementation focuses on IAU-76 and other pre-2000 models.

## Future Plans
- Support for TEME
- Support for GCRF (vice the current J2000)
- User defined EOP values, particularly for date values beyond what is included.
- More complete leap second handling
- Refactoring for use with StaticArrays.jl

## License Note

This package uses routines and computations derived from software provided by SOFA under license (see the LICENSE); and does not itself constitute software provided by and/or endorsed by SOFA.
