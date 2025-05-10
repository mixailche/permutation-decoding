# Permutation decoding of polar subcodes

Program implementation of the permutation-based decoder of polar subcodes and the builder of polar subcodes for permutation decoding with Arikan kernel. Currently supports polar EBCH subcodes and a randomized construction. Supported decoding algorithms are SC, SCL and their affine permutation based versions (not necessary authomorhism ensemble decoding).

## Build

Can be built using either Visual Studio or CMake.

#### Build via CMake

In the project root:
```
> mkdir build
> cd build
> cmake ..
> cmake --build .
```

## Run

There are two executables produced: `Builder.exe` and `Decoder.exe`. Both accept space-separeted list of CLI arguments, each one according to the format: `-name value`. Some arguments are optional - if they are not specified, the defualt value is used.\
`Builder.exe` constructs the polar subcode and writes its specification to the file. `Decoder.exe` consumes the specification file and runs the simulation of the given algorithm for AWGN channel and BPSK modulation, returning the obtained FER and the average number of operations with real numbers.

### Builder

Usage: `Builder.exe <arguments>`. Expected arguments:
+ `-out` - path to the output file
+ `-code` - type of the code; allowed values are:
  + `ebch` - construct a polar subcode of an EBCH code
  + `rand` - construct a randomized polar subcode
  + `rand-perm` - construct a randomized polar subcode for permutation decoding
+ `-len` - code length; integer, must be a power of 2
+ `-dim` - code dimension; integer
+ Parameter for estimating the virtual bit subchannels error probabilities. Exactly one of arguments must be specified:
  + `-BEC` - target erasure probability in the binary erasing channel; real number between 0 and 1
  + `-gauss` - target SNR in AWGN channel (gaussian approximation for density evolution is used); real number
+ For `-code rand` or `-code rand-perm`:
  + `-a` - number of type-A dynamic frozen symbols; integer
  + `-b` - number of type-B dynamic frozen symbols; integer
+ For `-code ebch`:
  + `-d` - design distance of the parent EBCH code
+ For `-code rand-perm`:
  + `-blocks` - block structure for affine authomorphisms of base code; comma separated sequence of integers, their sum must equal $log2(N)$ where $N$ is the code length

#### Example
```
Builder.exe \
  -out my/spec/file \
  -code rand-perm \
  -len 256 -dim 128 \
  -gauss 2.0 \
  -a 8 -b 32 \
  -blocks 1,1,1,1,1,3
```

### Decoder

Usage: `Decoder.exe <arguments>`. Expected arguments:
+ `-spec` - path to the specification file
+ `-alg` - decoding algorithm; allowed values are:
  + `sc` - successive cancellation (SC) decoding
  + `scl` - successive cancellation list (SCL) decoding
  + `perm-sc` - permutation-based SC decoding
  + `perm-scl` - permutation-based SCL decoding
+ `-snr` - SNR for the noise simulation; real number
+ `-iter` - max number of iterations; integer; **optional - default value is `1000000`**
+ `-errors` - max number errors while simulating; integer; **optional - default value is `100`**
+ `period` - number of iterations between intermediate result writing; **optional - default value is `100`**
+ For `-alg scl` and `-alg perm-scl`:
  + `-L` - max number of paths in the SCL decoder; integer
+ For `perm-sc` and `perm-scl`:
  + `-nperms` - number of permutations; integer
  + `-blocks` - block structure for affine permutations; same format as for `-blocks` argument of `Builder.exe`
  + Parameters for permutation list building. **Both arguments are optional having default value `0`**. They stem from the decomposition of every equivalence class representative matrix: $A=PU$ where $P$ is a block permutation and $U$ is an upper-triangular matrix.
    + `-dp` - minimum distance between $P$-matrices in the decompositions of distinct permutations used by decoder; integer 
    + `-du` - minimum distance between $U$-matrices; integer

If the value of the `-nperms` argument is greater than the total number of the permutations equivalence classes, program stucks in an infinite loop.

#### Example
```
Decoder.exe \
  -spec project_root/examples/spec_rand-perm_256_128_8_0 \
  -alg perm-scl \
  -snr 2.5 -iter 100000000 \
  -L 8 -nperms 8 \
  -blocks 1,1,1,1,1,3
```

### Specification file formal

First line consists of the code length and dimension. The others describe freezing constraints. The line corresponding to the equation
$$u_{i_1}+u_{i_2}+\dots+u_{i_{k-1}}=u_{i_k}$$
has the following format:
```
k i_1 i_2 ... i_k
```

