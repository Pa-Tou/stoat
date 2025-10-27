import bdsg
import bdsg_helper

graph_base="data/binary/pg.full.pg"
graph = bdsg.bdsg.PackedGraph()
graph.deserialize(graph_base)

n = graph.get_handle(4225)

samples_ins = bdsg_helper.get_paths_on_node(graph, n, True)

h1 = [x  for x in samples_ins if x.endswith("_h1")]
h0 = [x  for x in samples_ins if x.endswith("_h0")]

print(len(h1))
print(len(h0))