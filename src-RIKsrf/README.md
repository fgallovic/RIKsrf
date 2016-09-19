#Source files of individual codes
----------------

The evaluation of slip rate functions starts with generating rupture times using `RIKrandomruptvel`.
Then the main code `RIKsrfRandomruptvel` is to be run afterwords.

Appended script `compile.sh` shows how the codes can be compiled.

Credits:
 - Subroutines from Numerical recipes (FFT, random numbers, spline interpolation).
 - Subroutine XAPIIR (IIR filter design and implementation) by Dave Harris (1990).
