# fast-bve

Compile with `[mpi c++ compiler] -std=c++14 -O3 [lapack linker] driver.cpp -o driver`

On laptop, compile with `mpic++ -O3 -std=c++17 -framework Accelerate driver.cpp -o driver`

On Cheyenne, if using filesystem,
`module load gnu`
then compile with
`mpicxx -O3 -std=c++17 -lstdc++fs driver.cpp -o driver -L/glade/u/apps/ch/opt/netlib/3.9.0/lib64 -llapack -lblas -lgfortran`
otherwise
`mpicxx -O3 -std=c++14 driver.cpp -o driver`

On Derecho,
`mpicxx -O3 -std=c++14 -qmkl driver.cpp -o driver`

Modify things in namelist.txt
