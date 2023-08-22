# fast-bve

Compile with `[mpi c++ compiler] -std=c++14 -O3 [lapack linker] driver.cpp -o driver`

On laptop, compile with `mpic++ -O3 -std=c++14 -framework Accelerate driver.cpp -o driver`

On Derecho,
`mpicxx -O3 -std=c++14 -qmkl driver.cpp -o driver`

To compile with BLTC/BLTDD,
`mpicxx -O3 -framework Accelerate -std=c++14 -L../BVE-BaryTree/Build/lib/ -Wl,-rpath,../BVE-BaryTree/build/lib/ -lBaryTree_cpu bltc_rhs.cpp -o bltc_rhs`
and on Derecho
`mpicxx -O0 -std=c++14 -qmkl -L../BVE-BaryTree/build/lib/ -Wl,-rpath=../BVE-BaryTree/build/lib/ -lBaryTree_cpu bltc_rhs.cpp -o bltc_rhs`

Modify things in namelist.txt
