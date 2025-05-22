if ( NOT BUILD_TYPE )
    message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
    set(BUILD_TYPE BIG)
endif()

set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fno-automatic -fallow-argument-mismatch")

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build
    set(ENV{FC} mpif90)
    set(CMAKE_Fortran_COMPILER mpif90)
    set(USER_Fortran_FLAGS_RELEASE "-fconvert=little-endian -O3 -ffast-math -mtune=native -march=native")
    add_definitions(-DUSE_MPI -DUSE_MPI_IO)
    set(CMAKE_BUILD_TYPE RELEASE)
    
else() # compiler for serial build
    set(ENV{FC} gfortran)
    set(CMAKE_Fortran_COMPILER gfortran)
    
    if    ( ${BUILD_TYPE} STREQUAL "BIG" )
        set(USER_Fortran_FLAGS_RELEASE "-fconvert=big-endian -ffpe-summary=none -O3 -ffast-math -mtune=native -march=native")
        set(CMAKE_BUILD_TYPE RELEASE)
        
    elseif( ${BUILD_TYPE} STREQUAL "LITTLE" )
        set(USER_Fortran_FLAGS_RELEASE "-fconvert=little-endian -ffpe-summary=none -O3 -ffast-math -mtune=native -march=native")
        set(CMAKE_BUILD_TYPE RELEASE)
        
    elseif( ${BUILD_TYPE} STREQUAL "PROFILE" )
        set(USER_Fortran_FLAGS_DEBUG "-fconvert=little-endian -O0 -pg -ggdb -ffpe-summary=none")
        set(CMAKE_BUILD_TYPE DEBUG)

    else()
        set(USER_Fortran_FLAGS_DEBUG "-fconvert=little-endian -O0 -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow")#,underflow,precision,denormal")
        add_definitions(-D_DEBUG)
        set(CMAKE_BUILD_TYPE DEBUG)
        
    endif()
    
endif()

add_definitions(-DUSE_FFTW)
#set(FFTW_INCLUDE_DIR   "/usr/local/include")
#set(FFTW_LIB           "/usr/local/lib/libfftw3.a")
set(FFTW_LIB           "-lfftw3")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

add_definitions(-DUSE_NETCDF)
set(NC_INCLUDE_DIR     "/sw/buster-x64/io/netcdf-c-4.7.4-fortran-4.5.2-cxx4-4.3.1-cxx-4.2-gccsys/include")
set(NC_LIB             "-L/sw/buster-x64/io/netcdf-c-4.7.4-fortran-4.5.2-cxx4-4.3.1-cxx-4.2-gccsys/lib -Wl,-rpath -Wl,/sw/buster-x64/io/netcdf-c-4.7.4-fortran-4.5.2-cxx4-4.3.1-cxx-4.2-gccsys/lib -lnetcdff")
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${NC_INCLUDE_DIR})
set(LIBS ${LIBS} ${NC_LIB})
