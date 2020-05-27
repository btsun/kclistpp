#ifndef _MAX_FLOW_DOUBLE_HPP_
#define _MAX_FLOW_DOUBLE_HPP_

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
using namespace std;
using namespace boost;

class Network {
public:
  typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
  typedef adjacency_list < vecS, vecS, directedS,
    property < vertex_name_t, std::string,
    property < vertex_index_t, long,
    property < vertex_color_t, boost::default_color_type,
    property < vertex_distance_t, long,
    property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,

    property < edge_capacity_t, double,
    property < edge_residual_capacity_t, double,
    property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;
  typedef Traits::vertex_descriptor Vertex;
  typedef Traits::edge_descriptor Edge;
  Network();
  property_map<Graph, edge_capacity_t>::type capacity;
  property_map<Graph, edge_residual_capacity_t>::type residual_capacity;
  property_map<Graph, vertex_color_t>::type color;
  Vertex AddVertex();
  Edge AddEdge(Vertex &v1, Vertex &v2, const double capacity);
  double MaxFlow(Vertex &s, Vertex &t);
  void BfsOnResidualGraph(Vertex &s);

private:
  Graph g;
  property_map < Graph, edge_reverse_t >::type rev;
  struct NonSaturatedEdges;
};

#endif
