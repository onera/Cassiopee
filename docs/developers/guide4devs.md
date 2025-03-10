For developers:
===============

Editor:
-------
- indent with spaces (2 or 4 depending on file complexity). Dont use tabs.
- use utf8/lf encoding
- respect typo for commas: `t1, t2, a3`.
- respect type for two points: `a:`

Python Functions:
-----------------
- no snake (my_function) in function names or arguments. Use Camel (myFunction).
- short comment string in function header
- function must have a test before commit
- in place function starts with `_` and return None
- internal functions ends with `__`
- function must have a copyRef counterpart that calls in place function
- no IO in function (must work on input t and return t or a copy of t)
- try to unify argument names with existing functions
- complexifying an existing function, ie adding an argument or modifying the meaning of an argument, must be discussed
- check that a function performs correctly on `FlowSolutionNodes`, `FlowSolution#Centers`, `ZoneBC` and `ZoneGridConnectivity`
- if a function is fully operational, write the doc
- always run validCassiopee in full against the global test database before commit

Tests:
------
- no snake in test file name.
- first line of test shoud be `# - functionName (pyTree) -`
- seq tests finishes by _t1, _t2.
- parallel test finishes by _m1, _m2 and are run on 2 procs.
- in tests, dont use input files, create the test case in the script.
- a test must run in less than 10 seconds (ideally 1 sec).

OpenMP:
-------
- for computationally intensive loops, use omp parallelisation wherever possible. For example:

```c
#pragma omp parallel
{
    E_Int localInt;
    E_Float localFloat;
    #pragma omp for
    for (E_Int i = 0; i < ncells; i++)
    {
        localInt = i;
        localFloat = 12./i;
       /* ... */
    }
}
```

- the `pragma omp parallel` clause creates threads and has a slight cost: group `pragma omp for` clauses within its scope if possible.
- the `pragma omp for` clause splits the for loop into equal chunks
- in more complex cases, you can get the thread id and number of threads as:

```c
#pragma omp parallel
{
  E_Int ithread = __CURRENT_THREAD__;
  E_Int numThreads = __NUMTHREADS__; 
  /* ... */
}
```
