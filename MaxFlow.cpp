#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
#include "MaxFlow.hpp"
using namespace std;
using namespace boost;

Network::Network(): g(), rev(get(edge_reverse, g)) {}

Network::Vertex Network::AddVertex() {
  return add_vertex(g);
}

void Network::AddEdge(Vertex &v1, Vertex &v2, const long capacity) {
  Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
  Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
  put(edge_capacity, g, e1, capacity);
  rev[e1] = e2;
  rev[e2] = e1;
}

unsigned long long Network::MaxFlow(Vertex &s, Vertex &t) {
  return push_relabel_max_flow(g, s, t);
}
