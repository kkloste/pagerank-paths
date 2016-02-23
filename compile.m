function compile
mypath = fileparts(mfilename('fullpath'));
curdir = pwd();
try
    cd(mypath);
    if ismac
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims ppr_fast_grid_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims pprgrow_mex.cc
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims ppr_paths_rho_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims ./util/graphcutsweep_val.c
    else
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims ppr_fast_grid_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims pprgrow_mex.cc
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims ppr_paths_rho_mex.cpp
        mex -O CXXFLAGS="\$CXXFLAGS -Wall -I. -std=c++11" -largeArrayDims ./util/graphcutsweep_val.c
    end
    cd(curdir);
catch me
    cd(curdir)
    rethrow(me)
end
