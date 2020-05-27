# General Info

Approximate and exact algorithms for the k-clique densest subgraph problem.
These algorithms are described in "KClist++: A Simple Algorithm for Finding k-Clique Densest Subgraphs in Large Graphs" (Bintao Sun, Maximilien Danisch, T-H. Hubert Chan, Mauro Sozio; PVLDB 2020).

This repository includes the implementations of the following algorithms:
- Seq-kClist++: kCDensestNoAtom.c
- Sim-kClist++: kCDensestCmpSync.c
- Seq-Sampling++: kCDensestSamplingkClistpp.c
- Our exact algorithm: kCDensestMem.c
- The competitor MaxFlow: kCDensestG.cpp
- The competitor MaxFlow-Sampling: kCDensestSamplingMaxFlow.c

These codes are adopted from the codes in https://github.com/maxdan94/kClist.

# Note

Most of the programs use OpenMP for parallelization. Feel free to remove the relevant codes if you do not need them.

Some comments of the codes refer to the sequential update as "Frank-Wolfe" (which might be rather confusing and misleading). The sequential update using the "++" operator can be viewed as a sequential, asynchronous version of the Frank-Wolfe algorithm. See the paper for more details.

# Initial Contributor

Bintao Sun (btsun@connect.hku.hk)
