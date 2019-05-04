#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/properties.hpp>
#include <vector>

using namespace boost;
using namespace std;
    
class phase_t {
    public :
        unsigned int bench_id;
        unsigned int alloc;

        phase_t();
        phase_t(unsigned int bench_id,unsigned int alloc);
};
phase_t::phase_t(unsigned int bench_id,unsigned int alloc) : bench_id{bench_id}, alloc{alloc} {}
phase_t::phase_t() : bench_id{0}, alloc{0} {}
enum  RulesEdgeType {RULE1,RULE2,RULE3,RULE4};
typedef vector<phase_t>  alloc2_t;

// typedef struct NodeType {
//     string key;
//     string value;
// };
// typedef struct EdgeType {
//     string key;
//     string value;
// };

typedef adjacency_list<vecS,vecS,directedS,alloc2_t,RulesEdgeType> Graph;
// typedef adjacency_list<vecS,vecS,directedS,NodeType,EdgeType> Graph;
    
int main() {
    Graph g;
    // Then fill the graph
    add_edge(
        add_vertex(alloc2_t{phase_t(0,1),phase_t(0,2),phase_t(0,4)}, g),
        add_vertex(alloc2_t{phase_t(0,3),phase_t(0,6),phase_t(0,7)}, g),
        RULE1, g
    );
    // add_edge(
    //     add_vertex({"A","Apple"},g),
    //     add_vertex({"B","Ball"},g),
    //     {"prec","relation"},g
    // );

    // property_map<Graph, vertex_index_t>::type vim(get(vertex_index, g));
    // dynamic_properties dp;
    // dp.property("color",get(&NodeType::key,g));
    // dp.property("node_id",get(vertex_index,g));
    // write_graphviz_dp(std::cout, g, dp);
    write_graphviz(cout,g);
}