The file dipole.cpp is written with the help of ChatGPT in the Antenna course to teach students that it is possible to wire computational electromagnetics codes without being an expert in the programming language.

Run
g++ dipole.cpp -o dipoleCPP -std=c++20
The run ./dipoleCPP

Dr. Paknys' files can be downloaded here: https://www.wiley.com//legacy/wileychi/paknys/fortran.html?type=SupplementaryMaterial
I have put the file here.
For Dr. Paknys' codes in Fortran, simply run:
gfortran constants_m.f moler_m.f wire.f -o dipoleFortra
Then run ./dipoleFortran

Notice that the input impedance by both methods are around 75.8 - j 2.52 Ohms
