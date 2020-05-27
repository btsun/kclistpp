#ifndef _BINARY_HEAP_H_
#define _BINARY_HEAP_H_

typedef struct {
  unsigned key;
  unsigned value;
} KeyValuePair;

typedef struct {
  unsigned n_max; // Max number of nodes
  unsigned n; // Number of nodes
  unsigned *p; // Pointers to nodes
  KeyValuePair *kv; // Nodes
} BinaryHeap;

void BH_Update(BinaryHeap *heap, unsigned key);

KeyValuePair BH_PopMin(BinaryHeap *heap);

BinaryHeap *BH_MakeHeap(unsigned n, unsigned *v);

void BH_FreeHeap(BinaryHeap *heap);

#endif
