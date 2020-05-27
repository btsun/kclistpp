/*
Info:
  This program corresponds to the exact algorithm in the PVLDB 2020 paper.
  Feel free to use these lines as you wish.
  This program enumerates all k-cliques, store them in main memory, and apply
the "++" operator repeatedly to find the k-clique densest subgraph, until the
suspected k-clique densest subgraph passes the optimality test based on either
the improved Goldberg's condition or a max-flow.
  This program can handle both the case k = 2 (where all edges are treated as
the k-cliques) and the case k >= 3 (where the subroutine to list all k-cliques,
kClist, is executed once). Note again that all k-cliques are stored in main
memory, consuming super-linear space. One advantage, however, is that we can
shuffle all the cliques to prevent the cliques containing the same node from
coming in batch.

To compile:
  "g++ kCDensestMem.c BinaryHeap.c Graph.c MaxFlow.cpp -O3 -o kCDensestMem -lm -fopenmp"

To execute:
  "./kCDensestMem p k edgeListFileName tag"
  p is the number of threads.
  k is the size of a clique considered as in "k-clique".
  edgeListFileName is the name of the file that contains the graph. Each line of
the file contains one edge represented by two integers separated by a space.
  tag is a string specifying the dataset (e.g., "dblp"), which is used to
generate the output file name.

Output:
  A series of suspected k-clique densest subgraphs. One record per line,
containing
  - number of iterations of sequential updates run so far (always a power of 2);
  - the number of nodes in the suspected k-clique densest subset;
  - the k-clique density of the suspected k-clique densest subset;
  - the time elapsed since the beginning of the execution.
  When the exact solution is eventually found, the program additionally prints
  - the number of edges in the k-clique densest subgraph;
  - the optimality test that is passed ("Goldberg" or "Max Flow");
  - the number of max-flow calls.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include "Graph.h"
#include "MaxFlow.hpp"

unsigned MAX_CLIQUES = 100000000; // Maximum number of cliques for memory allocation; will increase if needed

static int UnsignedCmp(const void *a, const void *b) {
  return (long long)*(unsigned *)a - (long long)*(unsigned *)b;
}

Subgraph *AllocSubgraph(Graph *g, unsigned char k) {
  Subgraph *sg = (Subgraph *)malloc(sizeof(Subgraph));
  sg->n = (unsigned *)calloc(k, sizeof(unsigned));
  sg->d = (unsigned **)malloc(k * sizeof(unsigned *));
  sg->adj = (unsigned *)malloc(g->core * g->core * sizeof(unsigned));
  sg->label = (unsigned char *)calloc(g->core, sizeof(unsigned char));
  sg->nodes = (unsigned **)malloc(k * sizeof(unsigned *));
  sg->core = g->core;
  for (unsigned i = 1; i < k; ++i){
    sg->d[i] = (unsigned *)malloc(g->core * sizeof(unsigned));
    sg->nodes[i] = (unsigned *)malloc(g->core * sizeof(unsigned));
  }
  return sg;
}

static unsigned *id_sg2g = NULL, *id_g2sg = NULL; // to improve (???)
#pragma omp threadprivate(id_g2sg, id_sg2g)

void MakeSubgraph(Graph *g, unsigned u, unsigned v, Subgraph *sg, unsigned char k) {
  if (id_sg2g == NULL){
    id_g2sg = (unsigned *)malloc(g->n * sizeof(unsigned));
    id_sg2g = (unsigned *)malloc(g->core * sizeof(unsigned));
    for (unsigned i = 0; i < g->n; ++i) {
      id_g2sg[i] = UINT_MAX;
    }
  }

  for (unsigned i = 0; i < sg->n[k - 1]; ++i) {
    sg->label[i] = 0;
  }

  for (unsigned i = g->cd[v]; i < g->cd[v + 1]; ++i) { // For each out-neighbor of v
    id_g2sg[g->adj[i]] = UINT_MAX - 1;
  }

  unsigned j = 0;
  for (unsigned i = g->cd[u]; i < g->cd[u + 1]; ++i) { // For each out-neighbor of u
    unsigned x = g->adj[i];
    if (id_g2sg[x] == UINT_MAX - 1) {
      id_g2sg[x] = j;
      id_sg2g[j] = x;
      sg->label[j] = k - 2;
      sg->nodes[k - 2][j] = j;
      sg->d[k - 2][j] = 0; // New degrees
      ++j;
    }
  }
  sg->n[k - 2] = j;

  for (unsigned i = 0; i < sg->n[k - 2]; ++i) { // Reorder adjacency list and compute new degrees
    unsigned x = id_sg2g[i];
    for (unsigned l = g->cd[x]; l < g->cd[x + 1]; ++l) {
      unsigned y = g->adj[l];
      j = id_g2sg[y];
      if (j < UINT_MAX - 1) {
        sg->adj[sg->core * i + sg->d[k - 2][i]++] = j;
      }
    }
  }

  for (unsigned i = g->cd[v]; i < g->cd[v + 1]; ++i) {
    id_g2sg[g->adj[i]] = -1;
  }
}

// Clique-density-friendly decomposition

unsigned *cknodes; // Nodes of a clique
#pragma omp threadprivate(cknodes)
double *rho;
double *alpha;
double *rho_tentative;
unsigned *level;
unsigned *reordered;
unsigned *ck; // List of all cliques
unsigned *p_ckend; // Pointer to the end of ck[]
unsigned long long cnt_clique;
typedef enum {FRANK_WOLFE = 2,
              PAVA_PREPROCESS = 3} task_t;

void AllocCdf(Graph *g, unsigned k) {
  rho = (double *)malloc(g->n * sizeof(double));
  rho_tentative = (double *)malloc(g->n * sizeof(double));
  level = (unsigned *)malloc(g->n * sizeof(unsigned));
  reordered = (unsigned *)malloc(g->n * sizeof(unsigned));
}

inline int CDF_RerunFrankWolfeCmp(const unsigned u, const unsigned v) {
  if (level[u] < level[v]) return -1;
  if (level[u] > level[v]) return 1;
  if (rho[u] > rho[v]) return -1; // A node with larger rho value is "smaller"!
  if (rho[u] < rho[v]) return 1;
  return 0;
}

void CDF_FrankWolfeUpdateRates(int clique_size, unsigned *p_cknodes, double *p_alpha) {
  // Water-filling
  /*for (unsigned i = clique_size; i > 0; --i)
    for (unsigned j = 0; j + 1 < i; ++j)
      if (rho[cknodes[j]] > rho[cknodes[j + 1]]) {
        unsigned tmp = cknodes[j];
        cknodes[j] = cknodes[j + 1];
        cknodes[j + 1] = tmp;
      }
  double budget = 1.0;
  for (unsigned i = 0; i < clique_size; ++i) {
    double val = budget / (i + 1);
    if (i + 1 < clique_size && (rho[cknodes[i + 1]] - rho[cknodes[i]]) * (i + 1) < budget)
      val = rho[cknodes[i + 1]] - rho[cknodes[i]];
    for (unsigned j = 0; j <= i; ++j) {
      #pragma omp atomic
      rho[cknodes[j]] += val;
    }
    budget -= val * (i + 1);
  }*/
  unsigned node_index = 0;
  for (unsigned i = 1; i < clique_size; ++i) {
    if (rho[p_cknodes[node_index]] > rho[p_cknodes[i]])
      node_index = i;
  }
  #pragma omp atomic
  rho[p_cknodes[node_index]] += 1.0;
  #pragma omp atomic
  p_alpha[node_index] += 1.0;
}

void CDF_PavaPreprocessUpdateRates(int clique_size, unsigned *cknodes) {
  unsigned node_getting_weight = cknodes[0];
  for (unsigned l = 1; l < clique_size; ++l)
    if (level[cknodes[l]] > level[node_getting_weight])
      node_getting_weight = cknodes[l];
  #pragma omp atomic
  rho_tentative[level[node_getting_weight]] += 1.0;
}

void CDF_CliqueScan(Graph *g, unsigned char k, task_t task) {
  #pragma omp parallel for
  for (unsigned long long i = 0; i < cnt_clique; ++i) {
  //  for (unsigned j = 0; j < k; ++j)
  //    cknodes[j] = ck[i * k + j];
    switch (task) {
      case FRANK_WOLFE: {
        CDF_FrankWolfeUpdateRates(k, ck + i * k, alpha + i * k);
        break;
      }
      case PAVA_PREPROCESS: {
        CDF_PavaPreprocessUpdateRates(k, ck + i * k);
        break;
      }
    }
  }
}

void CDF_CliqueEnumThread(Subgraph *sg,
                          unsigned char clique_size,
                          unsigned char l) {
  if (clique_size == 3) {
    for (unsigned i = 0; i < sg->n[1]; ++i) {
      unsigned u = sg->nodes[1][i];
      cknodes[0] = id_sg2g[u];
      #pragma omp critical
      {
        if (cnt_clique >= MAX_CLIQUES) {
          MAX_CLIQUES *= 2;
          ck = (unsigned *)realloc(ck, MAX_CLIQUES * clique_size * sizeof(unsigned));
          p_ckend = ck + cnt_clique * clique_size;
        }
        for (unsigned j = 0; j < clique_size; ++j)
          *(p_ckend++) = cknodes[j];
        ++cnt_clique;
      }
    }
    return;
  }
  if (l == 2) {
    for (unsigned i = 0; i < sg->n[2]; ++i) {
      unsigned u = sg->nodes[2][i];
      cknodes[1] = id_sg2g[u];
      for (unsigned j = u * sg->core, end = u * sg->core + sg->d[2][u]; j < end; ++j) {
        unsigned v = sg->adj[j];
        cknodes[0] = id_sg2g[v];
        #pragma omp critical
        {
          if (cnt_clique >= MAX_CLIQUES) {
            MAX_CLIQUES *= 2;
            ck = (unsigned *)realloc(ck, MAX_CLIQUES * clique_size * sizeof(unsigned));
            p_ckend = ck + cnt_clique * clique_size;
          }
          for (unsigned k = 0; k < clique_size; ++k)
            *(p_ckend++) = cknodes[k];
          ++cnt_clique;
        }
      }
    }
    return;
  }

  for (unsigned i = 0; i < sg->n[l]; ++i) { // Enumerate in reverse order. Very confusing! "++i" is actually the reverse order.
    unsigned u = sg->nodes[l][i];
    cknodes[l - 1] = id_sg2g[u];

    sg->n[l - 1] = 0;
    unsigned end = u * sg->core + sg->d[l][u];
    for (unsigned j = u * sg->core; j < end; ++j) { // Relabel nodes and forming U'.
      unsigned v = sg->adj[j];
      if (sg->label[v] == l) {
        sg->label[v] = l - 1;
        sg->nodes[l - 1][sg->n[l - 1]++] = v;
        sg->d[l - 1][v] = 0; // New degrees
      }
    }
    for (unsigned j = 0; j < sg->n[l - 1]; ++j) { // Reorder adjacency list and compute new degrees
      unsigned v = sg->nodes[l - 1][j];
      for (unsigned k = sg->core * v, end = sg->core * v + sg->d[l][v]; k < end; ++k) {
        unsigned w = sg->adj[k];
        if (sg->label[w] == l - 1) {
          ++sg->d[l - 1][v];
        }
        else{
          sg->adj[k--] = sg->adj[--end];
          sg->adj[end] = w;
        }
      }
      qsort(sg->adj + sg->core * v, sg->d[l - 1][v], sizeof(unsigned), UnsignedCmp); // Sort the nodes in reverse order
    }

    CDF_CliqueEnumThread(sg, clique_size, l - 1);

    for (unsigned j = 0; j < sg->n[l - 1]; ++j) { // Restore labels
      unsigned v = sg->nodes[l - 1][j];
      sg->label[v] = l;
    }
  }
}

void CDF_CliqueEnum(Graph *g, unsigned char k) {
  Subgraph *sg;
  cnt_clique = 0;
  p_ckend = ck = (unsigned *)malloc(MAX_CLIQUES * k * sizeof(unsigned));
  #pragma omp parallel private(sg) reduction(+: cnt_clique)
  {
    cknodes = (unsigned *)malloc(k * sizeof(unsigned));
    sg = AllocSubgraph(g, k);
    #pragma omp for schedule(dynamic, 1) nowait
    for(unsigned i = 0; i < g->e; ++i) {
      cknodes[k - 1] = g->edges[i].s;
      cknodes[k - 2] = g->edges[i].t;
      MakeSubgraph(g, g->edges[i].s, g->edges[i].t, sg, k);
      CDF_CliqueEnumThread(sg, k, k - 2);
    }
    FreeSubgraph(sg, k);
  }
  ck = (unsigned *)realloc(ck, cnt_clique * k * sizeof(unsigned));
  alpha = (double *)malloc(cnt_clique * k * sizeof(double));
}

static int CDF_NodeCmp(const void *a, const void *b) {
  double d = rho[*(const unsigned *)a] - rho[*(const unsigned *)b];
  if (d > 0) return -1;
  if (d < 0) return 1;
  return 0;
}

typedef struct {
	unsigned n; // Total number of aggregated points
	unsigned *nag; // nag[i]: number of points aggregated in i
	double *val; // val[i]: value of the aggregated points
	// double *ub;
} IsotonicRegression;

// Pool Adjacent Violators Algorithm. Values to fit are stored in vect and n is the size of vect.
IsotonicRegression *CDF_Pava(double *vect, unsigned n) {
	IsotonicRegression *fit = (IsotonicRegression *)malloc(sizeof(IsotonicRegression));
	unsigned *nag = (unsigned *)malloc(n * sizeof(unsigned));
	double *val = (double *)malloc(n * sizeof(double));
	nag[0] = 1;
	val[0] = vect[0];
	unsigned j = 0;
	for (unsigned i = 1; i < n; ++i) {
		j += 1;
		val[j] = vect[i];
		nag[j] = 1;
		while (j > 0 && val[j] >= val[j - 1] * 0.999999) {
			val[j - 1] = (nag[j] * val[j] + nag[j - 1] * val[j - 1]) / (nag[j] + nag[j - 1]);
			nag[j - 1] += nag[j];
			--j;
		}
	}
	fit->n = j + 1;
	fit->nag = nag;
	fit->val = val;
	return fit;
}

IsotonicRegression *CDF_PavaPreprocess(Graph *g, unsigned char k) {
  for (unsigned i = 0; i < g->n; ++i)
    reordered[i] = i;
  qsort(reordered, g->n, sizeof(unsigned), CDF_NodeCmp); // Reorder the nodes by decreasing rho values
  for (unsigned i = 0; i < g->n; ++i)
    level[reordered[i]] = i;
  CDF_CliqueScan(g, k, PAVA_PREPROCESS);
  IsotonicRegression *partition = CDF_Pava(rho_tentative, g->n);
	for (unsigned j = 0, i = 0; j < partition->n; ++j)
		for (unsigned l = 0; l < partition->nag[j]; ++l, ++i)
			level[reordered[i]] = j;
  return partition;
}

bool CDF_CheckStability(Graph *g, unsigned subset_size, const unsigned char k) {
  for (unsigned i = 0; i < cnt_clique; ++i) {
    unsigned max_level = 0, max_level_cnt = 0;
    for (unsigned j = 0; j < k; ++j) {
      if (level[ck[i * k + j]] > max_level) {
        max_level = level[ck[i * k + j]];
        max_level_cnt = 1;
      } else if (level[ck[i * k + j]] == max_level) {
        ++max_level_cnt;
      }
    }
    double sum = 0;
    for (unsigned j = 0; j < k; ++j) {
      if (level[ck[i * k + j]] < max_level) {
        sum += alpha[i * k + j];
        #pragma omp atomic
        rho[ck[i * k + j]] -= alpha[i * k + j];
        alpha[i * k + j] = 0;
      }
    }
    for (unsigned j = 0; j < k; ++j) {
      if (level[ck[i * k + j]] == max_level) {
        #pragma omp atomic
        rho[ck[i * k + j]] += sum / max_level_cnt;
        alpha[i * k + j] += sum / max_level_cnt;
      }
    }
  }
  double prefix_min_rho = rho[reordered[0]];
  double suffix_max_rho = -1;
  for (unsigned i = 1; i < subset_size; ++i)
    if (prefix_min_rho > rho[reordered[i]])
      prefix_min_rho = rho[reordered[i]];
  for (unsigned i = g->n - 1; i >= subset_size; --i)
    if (suffix_max_rho < rho[reordered[i]])
      suffix_max_rho = rho[reordered[i]];
  return prefix_min_rho * 0.999999 > suffix_max_rho;
}

bool CDF_CheckDensestGoldberg(Graph *g, const unsigned n, const unsigned char clique_size) {
  qsort(reordered, n, sizeof(unsigned), CDF_NodeCmp); // Reorder the nodes by decreasing rho values
  unsigned long long m = 0;
  for (unsigned j = 0; j < n; ++j)
    m += (unsigned long long)(rho_tentative[j] + 0.5);
  double sum_rho = 0;
  double jck = 0; // j choose k
  bool skip = true;
  for (unsigned j = 1; j < n; ++j) {
    sum_rho += rho[reordered[j - 1]];
    if (skip) {
      if (j == clique_size)
        jck = 1;
      else if (j > clique_size)
        jck = (jck * j) / (j - clique_size);
      if (jck / j > (double)m / (double)n)
        skip = false, fprintf(stderr, "Jump to j = %u\n", j);
      else
        continue;
    }
    double ub = sum_rho / j;
    if (ub - (double)m / (double)n >= 1.0 / n / j &&
      ub - (double)m / (double)n >= (ceil((double)m * j / n) - (double)m * j / n) / j) {
      return false;
    }
  }
  return true;
}

bool CDF_CheckDensestMaxFlow(Graph *g, const unsigned n, const unsigned char clique_size, unsigned *p_cnt_max_flow) {
  ++(*p_cnt_max_flow);
  Network network;
  unsigned *id_in_network = (unsigned *)malloc(g->n * sizeof(unsigned));
  unsigned long long m = 0;
  vector<Network::Vertex> R;
  Network::Vertex s = network.AddVertex(), t = network.AddVertex();
  for (unsigned i = 0; i < g->n; ++i)
    id_in_network[i] = n;
  for (unsigned i = 0; i < n; ++i) {
    id_in_network[reordered[i]] = i;
    R.push_back(network.AddVertex());
  }
  for (unsigned i = 0; i < cnt_clique; ++i) {
    bool flag = true;
    for (unsigned j = 0; j < clique_size; ++j) {
      if (id_in_network[ck[i * clique_size + j]] >= n) {
        flag = false;
        break;
      }
    }
    if (flag) {
      ++m;
      Network::Vertex v = network.AddVertex();
      for (unsigned j = 0; j < clique_size; ++j)
        network.AddEdge(v, R[id_in_network[ck[i * clique_size + j]]], n);
      network.AddEdge(s, v, n);
    }
  }
  for (unsigned i = 0; i < n; ++i)
    network.AddEdge(R[i], t, m);
  free(id_in_network);
  return network.MaxFlow(s, t) >= m * n;
}

void ShuffleCliques(const unsigned k) {
  for (unsigned i = 1; i < cnt_clique; ++i) {
    unsigned rand_index = rand() % (i + 1);
    for (unsigned j = 0; j < k; ++j) {
      unsigned temp = ck[i * k + j];
      ck[i * k + j] = ck[rand_index * k + j];
      ck[rand_index * k + j] = temp;
    }
  }
}

void CDF_Main(const unsigned char k, Graph *g, FILE *ofp, time_t t0) {
  unsigned cnt_max_flow = 0;
  fprintf(ofp, "[Number of Iterations]\t[Number of Nodes]\t[k-Clique Density]\t[Time (seconds)]\t[Info]\n");
  AllocCdf(g, k);
  if (k >= 3) {
    CDF_CliqueEnum(g, k); // Collect all k-cliques
  } else {
    ck = (unsigned *)malloc((unsigned long long)(g->e) * k * sizeof(unsigned));
    alpha = (double *)malloc(g->e * k * sizeof(double));
    cnt_clique = g->e;
    for (unsigned long long i = 0; i < (unsigned long long)(g->e); ++i) {
      ck[i << 1] = g->edges[i].s;
      ck[(i << 1) + 1] = g->edges[i].t;
    }
  }
  ShuffleCliques(k);
  for (unsigned i = 0; i < g->n; ++i)
    rho[i] = 0;
  for (unsigned long long i = 0; i < cnt_clique * k; ++i)
    alpha[i] = 0;
  for (unsigned num_iter = 1; ; num_iter <<= 1) {
    fprintf(stderr, "Start: number of iterations = %u.\n", num_iter);
    // Step 1: run the Frank-Wolfe based algorithm for num_iter rounds
    for (unsigned t = num_iter / 2 + 1; t <= num_iter; ++t) {
      if (t % 10 == 0)
        fprintf(stderr, "Run round %u...\n", t);
      CDF_CliqueScan(g, k, FRANK_WOLFE);
    }
    // Step 2: give a tentative decomposition
    for (unsigned i = 0; i < g->n; ++i)
      rho_tentative[i] = 0;
    IsotonicRegression *partition = CDF_PavaPreprocess(g, k);
    fprintf(stderr, "Approximate densest subgraph: %u nodes, density = %f.\n", partition->nag[0], partition->val[0]);
  /*  FILE *ofp = fopen("rates.txt", "w");
    for (unsigned i = 0; i < g->n; ++i)
      fprintf(ofp, "r[%u] = %.12f\n", reordered[i], rho[reordered[i]]);
    fclose(ofp);*/
    // Step 3: Check stability and optimality
    if (CDF_CheckStability(g, partition->nag[0], k)) {
      fprintf(stderr, "The potential densest set is stable!\n");
      if (CDF_CheckDensestGoldberg(g, partition->nag[0], k)) {
        fprintf(stderr, "The first %u nodes forms a densest subgraph by criteria A!\n", partition->nag[0]);
        fprintf(ofp, "[Number of Iterations]\t[Stopping Condition]\t[Number of Nodes]\t[Number of Edges]\t[k-Clique Density]\t[Number of Max-Flow Calls]\t[Time (seconds)]\n");
        fprintf(ofp, "%u\tGoldberg\t%u\t%u\t%.12f\t%u\t", num_iter, partition->nag[0], CountEdges(g, partition->nag[0], reordered), partition->val[0], cnt_max_flow);
        break;
      }
      else if (CDF_CheckDensestMaxFlow(g, partition->nag[0], k, &cnt_max_flow)) {
        fprintf(stderr, "The first %u nodes forms a densest subgraph by criteria B!\n", partition->nag[0]);
        fprintf(ofp, "[Number of Iterations]\t[Stopping Condition]\t[Number of Nodes]\t[Number of Edges]\t[k-Clique Density]\t[Number of Max-Flow Calls]\t[Time (seconds)]\n");
        fprintf(ofp, "%u\tMax Flow\t%u\t%u\t%.12f\t%u\t", num_iter, partition->nag[0], CountEdges(g, partition->nag[0], reordered), partition->val[0], cnt_max_flow);
        break;
      }
      else {
        fprintf(stderr, "Cannot guarantee it is densest by either criteria A or criteria B.\n");
        fprintf(ofp, "%u\t%u\t%.12f\t%ld\tSTABLE BUT NOT DENSEST\n", num_iter, partition->nag[0], partition->val[0], time(NULL) - t0);
      }
    } else {
      fprintf(stderr, "The potential densest subset is not stable!\n");
      fprintf(ofp, "%u\t%u\t%.12f\t%ld\tNOT STABLE\n", num_iter, partition->nag[0], partition->val[0], time(NULL) - t0);
    }
  /*  ofp = fopen("rates_rerun.txt", "w");
    for (int i = 0; i < partition->nag[0]; ++i)
      fprintf(ofp, "r[%u] = %.12f\n", reordered[i], rho[reordered[i]]);
    fclose(ofp);*/
  }
}

int main(int argc, char **argv) {
  EdgeList *el;
  Graph *g;
  unsigned char k = atoi(argv[2]);
  char *file_name = argv[3];
  unsigned num_threads = atoi(argv[1]);
  omp_set_num_threads(num_threads);

  time_t t0, t1, t2;
  t0 = t1 = time(NULL);

  printf("Reading edgelist from file %s\n", file_name);
  el = ReadEdgeList(file_name);
  printf("Number of nodes = %u\n", el->n);
  printf("Number of edges = %u\n", el->e);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n",(t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  printf("Building the graph structure\n");
  SortByCore(el); // Do core decomposition and render degeneracy ordering to the nodes
  Relabel(el);
  g = MakeGraph(el);
  printf("Number of nodes (degree > 0) = %u\n", g->n);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  printf("Iterate over all cliques\n");
  char output_file_name[100] = "stat_exact_";
  strcat(output_file_name, argv[4]);
  strcat(output_file_name, "_");
  strcat(output_file_name, argv[1]);
  strcat(output_file_name, "_");
  strcat(output_file_name, argv[2]);
  strcat(output_file_name, ".txt");
  FILE *ofp = fopen(output_file_name, "w");
  try {
    CDF_Main(k, g, ofp, t0);
  } catch(std::exception &e) {
    fprintf(ofp, "%s\n", e.what());
  }

  printf("Number of %u-cliques: %llu\n", k, cnt_clique);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  FreeGraph(g);
  printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));
  fprintf(ofp, "%ld\n", t2 - t0);
  fclose(ofp);
  return 0;
}
