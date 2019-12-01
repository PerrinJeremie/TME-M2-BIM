cd GTL-1.2.4/
./configure --prefix=$HOME/GTL/
make
make install
cd ../tv-0.5/ncl-2.0/
make
cd ../../superTree
cp ../tv-0.5/ncl-2.0/src/*.cpp ncl-2.0/src/
cp ../tv-0.5/ncl-2.0/src/.deps/* ncl-2.0/src/.deps/
cp ../tv-0.5/ncl-2.0/src/*.h ncl-2.0/src/
make
