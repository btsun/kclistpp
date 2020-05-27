#define BOOST_DISABLE_ASSERTS
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include "MaxFlowDouble.hpp"
using namespace std;
using namespace boost;

Network::Network(): g(), rev(get(edge_reverse, g)), capacity(get(edge_capacity, g)), residual_capacity(get(edge_residual_capacity, g)) {}

Network::Vertex Network::AddVertex() {
  return add_vertex(g);
}

Network::Edge Network::AddEdge(Vertex &v1, Vertex &v2, const double capacity) {
  Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
  Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
  put(edge_capacity, g, e1, capacity);
  rev[e1] = e2;
  rev[e2] = e1;
  return e1;
}

double Network::MaxFlow(Vertex &s, Vertex &t) {
  return push_relabel_max_flow(g, s, t); // Boost library also provides boykov_kolmogorov_max_flow (needs "#include <boost/graph/boykov_kolmogorov_max_flow.hpp>")
}

struct Network::NonSaturatedEdges {
  NonSaturatedEdges() {}
  NonSaturatedEdges(property_map<Graph, edge_residual_capacity_t>::type residual_capacity):
    residual_capacity(residual_capacity) {}
  property_map<Graph, edge_residual_capacity_t>::type residual_capacity;
  bool operator ()(const Edge &e) const {
    return residual_capacity[e] > 1e-9;
  }
};

void Network::BfsOnResidualGraph(Vertex &s) {
  NonSaturatedEdges filter(get(edge_residual_capacity, g));
  filtered_graph<Graph, NonSaturatedEdges> fg(g, filter);
  color = get(vertex_color, g);
  boost::queue<Vertex> Q;
  default_bfs_visitor vis;
  breadth_first_search(fg, s, Q, vis, color);
}
