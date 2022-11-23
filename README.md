# PARALLEL INTEGRATION WITH TRAPEZOIDS/SIMPSON 

MPI is library for message passing interface standard designed to function on parallel computing architectures.


## Compilation



```bash
mpicc programName.c -o outputProgramName  -lm
```

## Usage

```c
mpirun -np num_proc  ./program [trapezoid/interval numbers] a b [integration function val:0,1,2]


```

## Contributing

Written by De Fazio Francesco for the parallel algorithms and distributed systems exam.
