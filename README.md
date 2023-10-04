# fast-bve

To compile with BLTC/BLTDD,
`mpicxx -O3 -framework Accelerate -std=c++14 -L../BVE-BaryTree/Build/lib/ -Wl,-rpath,../BVE-BaryTree/build/lib/ -lBaryTree_cpu bltc_rhs.cpp -o bltc_rhs`
and on Derecho
`mpicxx -O0 -std=c++14 -qmkl -L../BVE-BaryTree/build/lib/ -Wl,-rpath=../BVE-BaryTree/build/lib/ -lBaryTree_cpu bltc_rhs.cpp -o bltc_rhs`

Modify things in bin/namelist.txt

Clone repo, `mkdir build`, `cd build`, then `cmake ..`, then `make`. 
