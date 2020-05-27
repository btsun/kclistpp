#include <stdlib.h>
#include "BinaryHeap.h"

// C implementation of a min-heap (value of root is no greater than value of sons)

BinaryHeap *BH_Construct(unsigned n_max) {
  BinaryHeap *heap = (BinaryHeap *)malloc(sizeof(BinaryHeap));
  heap->n_max = n_max;
  heap->n = 0;
  heap->p = (unsigned *)malloc(n_max * sizeof(unsigned));
  for (unsigned i = 0; i < n_max; ++i) heap->p[i] = -1;
  heap->kv = (KeyValuePair *)malloc(n_max * sizeof(KeyValuePair));
  return heap;
}

void BH_Swap(BinaryHeap *heap, unsigned i, unsigned j) {
  KeyValuePair kv_tmp = heap->kv[i];
  unsigned p_tmp = heap->p[kv_tmp.key];
  heap->p[heap->kv[i].key] = heap->p[heap->kv[j].key];
  heap->kv[i] = heap->kv[j];
  heap->p[heap->kv[j].key] = p_tmp;
  heap->kv[j] = kv_tmp;
}

void BH_BubbleUp(BinaryHeap *heap, unsigned i) {
  unsigned j = (i - 1) >> 1;
  while (i > 0) {
    if (heap->kv[j].value > heap->kv[i].value) {
      BH_Swap(heap, i, j);
      i = j;
      j = (i - 1) >> 1;
    }
    else break;
  }
}

void BH_BubbleDown(BinaryHeap *heap) {
  unsigned i = 0, j1 = 1, j2 = 2, j;
  while (j1 < heap->n) {
    j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
    if (heap->kv[j].value < heap->kv[i].value) {
      BH_Swap(heap, i, j);
      i = j;
      j1 = (i << 1) + 1;
      j2 = j1 + 1;
      continue;
    }
    break;
  }
}

void BH_Insert(BinaryHeap *heap, KeyValuePair kv) {
  heap->p[kv.key] = (heap->n)++;
  heap->kv[heap->n - 1] = kv;
  BH_BubbleUp(heap, heap->n - 1);
}

void BH_Update(BinaryHeap *heap, unsigned key) { // Decrease the value of a node by 1
  unsigned i = heap->p[key];
  if (i != -1) {
    --(heap->kv[i].value);
    BH_BubbleUp(heap, i);
  }
}

KeyValuePair BH_PopMin(BinaryHeap *heap) {
  KeyValuePair kv = heap->kv[0];
  heap->p[kv.key] = -1;
  heap->kv[0] = heap->kv[--(heap->n)];
  heap->p[heap->kv[0].key] = 0;
  BH_BubbleDown(heap);
  return kv;
}

BinaryHeap *BH_MakeHeap(unsigned n, unsigned *v) {
  KeyValuePair kv;
  BinaryHeap *heap = BH_Construct(n);
  for (unsigned i = 0; i < n; ++i) {
    kv.key = i;
    kv.value = v[i];
    BH_Insert(heap, kv);
  }
  return heap;
}

void BH_FreeHeap(BinaryHeap *heap) {
  free(heap->p);
  free(heap->kv);
  free(heap);
}
