# OrbitPropagationLibrarySOFA

[![Build Status](https://github.com/ampalski/OrbitPropagationLibrarySOFA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ampalski/OrbitPropagationLibrarySOFA.jl/actions/workflows/CI.yml?query=branch%3Amain)

A pure Julia implementation of the IAU timing and coordinate tranformation functions, as published in the SOFA routines with supplemental material derived from Vallado. 

This package is still in development, and is intended to be used as the timing and coordinate backbone to an Earth orbit propagation library and not explicitly for astronomy the way SOFA is. As such, some function signatures will be different. 

In order to maintain compatibility with the USSF standards, initial implementation will focus on IAU-76 and other pre-2000 models.

This package uses routines and computations derived from software provided by SOFA under license (see the LICENSE); and does not itself constitute software provided by and/or endorsed by SOFA.
