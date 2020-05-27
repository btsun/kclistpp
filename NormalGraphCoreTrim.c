/*
Info:
  A preprocessing tool that trims an input graph while guaranteeting that the
maximal densest subgraph is untouched.
  Feel free to use these lines as you wish.
  This program does the same trimming as in kClistCoreTrim.c, but for the case
k = 2. In this case, we do not need to rely on the kClist subprocedure. See
kClistCoreTrim.c for more details.

To compile:
  "gcc NormalGraphCoreTrim.c BinaryHeap.c Graph.c -O3 -o NormalGraphCoreTrim".

To execute:
  "./NormalGraphCoreTrim edgelist.txt tag".
  - edgelist.txt is the name of the file that contains the graph. Each line of
the file contains one edge represented by two integers separated by a space.
  - tag is a string specifying the dataset (e.g., "dblp"), which is used to
generate the output file name.

Output:
  The program will output
  - the trimmed graph in trimmed_[tag]_2.txt, in the same format as
edgelist.txt.
  - the statistics of the trimmed graph and running time in
stat_trim_[tag]_2.txt.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "BinaryHeap.h"
#include "Graph.h"

// Compute degeneracy ordering.
// The larger the node's core number, the higher (smaller) the node's rank.
void TrimByCoreDecomp(EdgeList *el, char *trimg, char *stat, time_t t0) {
  unsigned n = el->n, e = el->e;
  unsigned *d = (unsigned *)calloc(n, sizeof(unsigned));
  unsigned *cd = (unsigned *)malloc((n + 1) * sizeof(unsigned));
  unsigned *adj = (unsigned *)malloc(2 * e * sizeof(unsigned));
  unsigned *deg_when_del = (unsigned *)malloc(n * sizeof(unsigned));
  unsigned *node_to_del = (unsigned *)malloc(n * sizeof(unsigned));
  bool *alive = (bool *)malloc(n * sizeof(bool));
  unsigned cnt_remaining_edges = e;
  for (unsigned i = 0; i < e; ++i) {
    ++d[el->edges[i].s];
    ++d[el->edges[i].t];
  }
  cd[0] = 0;
  for (unsigned i = 1; i < n + 1; ++i) {
    cd[i] = cd[i - 1] + d[i - 1];
    d[i - 1] = 0;
  }
  for (unsigned i = 0; i < e; ++i) {
    adj[cd[el->edges[i].s] + d[el->edges[i].s]++] = el->edges[i].t;
    adj[cd[el->edges[i].t] + d[el->edges[i].t]++] = el->edges[i].s;
  }

  BinaryHeap *heap = BH_MakeHeap(n, d);

  double max_rho = -1;
  unsigned r = 0;
  for (unsigned i = 0; i < n; ++i) {
    if (max_rho < (double)cnt_remaining_edges / (n - i))
      max_rho = (double)cnt_remaining_edges / (n - i);
    KeyValuePair kv = BH_PopMin(heap);
    node_to_del[i] = kv.key;
    deg_when_del[i] = kv.value;
    cnt_remaining_edges -= kv.value;
    for (unsigned j = cd[kv.key]; j < cd[kv.key + 1]; ++j) {
      BH_Update(heap, adj[j]);
    }
  }

  FILE *ofp = fopen(stat, "w");
  cnt_remaining_edges = e;
  for (unsigned i = 0; i < n; ++i)
    alive[i] = true;
  for (unsigned i = 0; i < n; ++i) {
    if (deg_when_del[i] < max_rho) {
      alive[node_to_del[i]] = false;
      cnt_remaining_edges -= deg_when_del[i];
    }
    else {
      fprintf(ofp, "[Number of Nodes]\t[Number of Edges]\t[Number of k-Cliques]\t[Time (seconds)]\n");
      fprintf(ofp, "%u\t%u\t%u\t", n - i, cnt_remaining_edges, cnt_remaining_edges);
      break;
    }
  }
  fclose(ofp);
  ofp = fopen(trimg, "w");
  for (unsigned i = 0; i < e; ++i) {
    if (alive[el->edges[i].s] && alive[el->edges[i].t])
      --cnt_remaining_edges, fprintf(ofp, "%u %u\n", el->edges[i].s, el->edges[i].t);
  }
  fclose(ofp);

  assert(cnt_remaining_edges == 0);

  BH_FreeHeap(heap);
  free(d);
  free(cd);
  free(adj);
  free(deg_when_del);
  free(node_to_del);
}

int main(int argc, char **argv) {
  EdgeList *el;
  char *file_name = argv[1];

  time_t t0, t1, t2;
  t0 = t1 = time(NULL);

  printf("Reading edgelist from file %s\n", file_name);
  el = ReadEdgeList(file_name);
  printf("Number of nodes = %u\n", el->n);
  printf("Number of edges = %u\n", el->e);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n",(t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  char trim_file_name[100] = "trimmed_";
  strcat(trim_file_name, argv[2]);
  strcat(trim_file_name, "_2.txt");
  char stat_file_name[100] = "stat_trim_";
  strcat(stat_file_name, argv[2]);
  strcat(stat_file_name, "_2.txt");
  printf("Trimming the graph by core decomposition\n");
  TrimByCoreDecomp(el, trim_file_name, stat_file_name, t0); // Do core decomposition trim the graph

  t2 = time(NULL);
  printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));
  FILE *ofp = fopen(stat_file_name, "a");
  fprintf(ofp, "%ld\n", t2 - t0);
  fclose(ofp);

  return 0;
}
