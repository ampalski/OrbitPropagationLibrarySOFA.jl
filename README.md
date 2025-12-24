# OrbitPropagationLibrarySOFA

[![Build Status](https://github.com/ampalski/OrbitPropagationLibrarySOFA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ampalski/OrbitPropagationLibrarySOFA.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

A pure Julia implementation of the IAU timing and coordinate tranformation functions, as published in the SOFA routines with supplemental material derived from Vallado. 

This package is at minimum viable product, and is intended to be used as the timing and coordinate backbone to an Earth orbit propagation library and not explicitly for astronomy the way SOFA is. As such, some function signatures will be different. As more functionality is added (see below), some function signatures will change.

In order to maintain compatibility with the USSF standards, initial implementation focuses on IAU-76 and other pre-2000 models.

## Usage Example

```
using OrbitPropagationLibrarySOFA
using StaticArrays
date = [2014.0, 8, 28, 6, 46, 24.461] #UTC
JD, _ = datevec2jdate(date)
JDUT1 = convert_jd(JD, :UT1)
JDTT = convert_jd(JDUT1, :TT)

rJ2000 = SA[23141.52, 35279.3, -5.05699]
rMOD = j20002mod76(r, JDTT)
rTOD = mod2tod76(rMOD, JDTT)
rPEF = tod2pef76(rTOD, JDUT1)
rITRF = pef2itrf76(rPEF, JD)
```

## Future Plans
- Support for GCRF (vice the current J2000)
- User defined EOP values, particularly for date values beyond what is included.
- More complete leap second handling

## License Note

This package uses routines and computations derived from software provided by SOFA under license (see the LICENSE); and does not itself constitute software provided by and/or endorsed by SOFA.
