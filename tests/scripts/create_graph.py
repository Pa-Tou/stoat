from bdsg.bdsg import HashGraph
import subprocess
import argparse
import os

def vg_process(filename) :
    pg_filename = filename.replace(".hg", ".pg")
    dist_filename = filename.replace(".hg", ".dist")
    subprocess.run(["vg", "convert", filename, "-p"], stdout=open(pg_filename, "wb"))
    subprocess.run(["vg", "index", pg_filename, "-j", dist_filename])

def vg_view(filename, max_node) :
    pg_filename = filename.replace(".hg", ".pg")
    string_vg = ["vg", "find", "-x", pg_filename, "-r", f"1:{str(max_node)}", "-c", "10", "|", "vg", "view", "-dp", "-", "|", "dot", "-Tsvg", "-o", filename.replace(".hg", ".svg")]
    string_vg = " ".join(string_vg)
    subprocess.run(string_vg, shell=True)

def create_simple_snp_graph(filename="simple_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "C", "G", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 4, 5]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 4, 5]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 6)
    print("simple snp graph created")

def create_3th_snp_graph(filename="3th_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "C", "T", "G", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[1], nodes[4])
    gr.create_edge(nodes[2], nodes[5])
    gr.create_edge(nodes[3], nodes[5])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[6])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 5, 6]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 5, 6]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 7)
    print("3th snp graph created")

def create_insert_snp_graph(filename="insert_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "C", "GTT", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 4, 5]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 4, 5]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 6)
    print("insert snp graph created")

def create_deletion_snp_graph(filename="deletion_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "C", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 4]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 5)
    print("deletion_snp created")

def create_insert_deletion_graph(filename="insert_deletion.hg") :
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "GTT", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 4]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 5)
    print("insert deletion graph created")

def create_loop_simple_graph(filename="loop_simple.hg") :

    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "A", "CG", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[2], nodes[2]) # LOOP
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])

    path = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 4, 5]:
        gr.append_step(path, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 4, 5]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 6)
    print("simple loop graph created")

def create_loop_double_graph(filename="loop_double.hg") :

    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "CG", "TC", "GG", "GA", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[5])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[2]) # 1er LOOP
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[2]) # 2eme LOOP
    gr.create_edge(nodes[4], nodes[6])
    gr.create_edge(nodes[5], nodes[6])
    gr.create_edge(nodes[6], nodes[7])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 5, 6, 7]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 3, 4, 6, 7]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 8)
    print("double loop graph created")

def create_loop_graph(filename="loop.hg") :

    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "A", "CG", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[2], nodes[1]) # LOOP
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])

    path = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 4, 5]:
        gr.append_step(path, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 4, 5]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 6)
    print("loop graph created")

def create_loop_plus_graph(filename="loop_plus.hg") :

    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "G", "AT", "C", "A", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[6])
    gr.create_edge(nodes[2], nodes[3]) # NESTED
    gr.create_edge(nodes[2], nodes[4]) # NESTED
    gr.create_edge(nodes[3], nodes[5]) # NESTED
    gr.create_edge(nodes[4], nodes[5]) # NESTED
    gr.create_edge(nodes[5], nodes[2]) # NESTED + LOOP
    gr.create_edge(nodes[5], nodes[7])
    gr.create_edge(nodes[6], nodes[7])
    gr.create_edge(nodes[7], nodes[8])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 6, 7, 8]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 3, 5, 7, 8]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 9)
    print("loop plus graph created")

def create_linear_path(filename="linear.hg"):
    gr = HashGraph()
    seqs = ["AAA", "T", "GGG", "C"]
    nodes = [gr.create_handle(s) for s in seqs]

    for i in range(len(nodes) - 1):
        gr.create_edge(nodes[i], nodes[i + 1])

    path = gr.create_path_handle("ref")
    for node in nodes:
        gr.append_step(path, node)

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 4)
    print("linear path graph created")

def create_snp_and_nested_snp_graph(filename="snp_and_nested_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "G", "AT", "C", "T", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[6])
    gr.create_edge(nodes[2], nodes[3]) # NESTED
    gr.create_edge(nodes[2], nodes[4]) # NESTED
    gr.create_edge(nodes[3], nodes[5]) # NESTED
    gr.create_edge(nodes[4], nodes[5]) # NESTED
    gr.create_edge(nodes[5], nodes[7])
    gr.create_edge(nodes[6], nodes[7])
    gr.create_edge(nodes[7], nodes[8])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 6, 7, 8]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 3, 5, 7, 8]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 9)
    print("snp and nested snp graph created")

def create_complex_ins_graph(filename="complex_ins.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "G", "A", "C", "T", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[6])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[1], nodes[7]) # BIG DELETION
    gr.create_edge(nodes[3], nodes[5])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[4], nodes[6])
    gr.create_edge(nodes[5], nodes[7])
    gr.create_edge(nodes[6], nodes[7])
    gr.create_edge(nodes[7], nodes[8])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 6, 7, 8]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 3, 5, 7, 8]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 9)
    print("snp_and_nested_plus created")

def create_4th_graph(filename="4th.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "CA", "TCG", "CTAT", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[2], nodes[5])
    gr.create_edge(nodes[3], nodes[5])
    gr.create_edge(nodes[4], nodes[6])
    gr.create_edge(nodes[5], nodes[6])
    gr.create_edge(nodes[6], nodes[7])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 4, 6, 7]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 5, 6, 7]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 8)
    print("4th graph created")

def create_repetition_graph(filename="repetition.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "GTT", "CAA", "TGG", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[6])
    gr.create_edge(nodes[1], nodes[5])
    gr.create_edge(nodes[2], nodes[5])
    gr.create_edge(nodes[3], nodes[5])
    gr.create_edge(nodes[4], nodes[5])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 5, 6]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 3, 4, 5, 6]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 7)
    print("repetition graph created")

def create_large_del_graph(filename="large_del.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "GTT", "A", "T", "CCC", "G", "TT", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[8])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[3], nodes[5])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[6])
    gr.create_edge(nodes[5], nodes[7])
    gr.create_edge(nodes[6], nodes[7])
    gr.create_edge(nodes[7], nodes[8])
    gr.create_edge(nodes[8], nodes[9])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 8, 9]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 4, 5, 6, 7, 8, 9]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 10)
    print("large deletion graph created")

def create_inversion_graph(filename="inversion.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "AT", "CAT", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[2], gr.flip(nodes[3])) # begin 3 -> end 4
    gr.create_edge(nodes[1], nodes[5]) # DEL
    gr.create_edge(gr.flip(nodes[3]), nodes[4]) # begin 4 -> begin 5
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[6])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 4, 5, 6]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 7)
    print("inversion graph created")

def create_nested_plus_graph(filename="nested_plus.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "A", "TT", "AC", "T", "C", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[7])
    gr.create_edge(nodes[2], nodes[3]) # nested
    gr.create_edge(nodes[2], nodes[4]) # nested
    gr.create_edge(nodes[3], nodes[5])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[6])
    gr.create_edge(nodes[6], nodes[7])
    gr.create_edge(nodes[5], nodes[7])
    gr.create_edge(nodes[7], nodes[8])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 5, 7, 8]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 9)
    print("nested_plus graph created")

def create_jean_graph(filename="jean.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "AT", "CAT", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[3], nodes[2]) # begin 4 -> end 3
    gr.create_edge(nodes[1], nodes[5]) # DEL
    gr.create_edge(gr.flip(nodes[3]), nodes[4]) # begin 4 -> begin 5
    gr.create_edge(nodes[4], gr.flip(nodes[4])) # loop 5 end -> 5 end
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[6])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 4, 5, 6]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 7)
    print("jean graph created")

def create_multicomponent_chain_graph(filename="multicomponent_chain.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "C", "CCCC", "TA", "CT", "TAG", "T", "C", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[10])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], gr.flip(nodes[6])) # end 4 -> end 5
    gr.create_edge(nodes[6], nodes[7])
    gr.create_edge(nodes[7], nodes[8])
    gr.create_edge(nodes[7], nodes[9])
    gr.create_edge(nodes[8], nodes[10])
    gr.create_edge(nodes[9], nodes[10])
    gr.create_edge(nodes[10], nodes[11])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 10, 11]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 12)
    print("multicomponent_chain graph created")

def create_looping_chain_graph(filename="looping_chain.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "TAG", "T", "C", "CAT", "AATT", "AT", "A", "T", "TA", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[3], nodes[5])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[6])
    gr.create_edge(nodes[6], nodes[7])
    gr.create_edge(nodes[7], nodes[8])
    gr.create_edge(nodes[7], nodes[9])
    gr.create_edge(nodes[8], nodes[10])
    gr.create_edge(nodes[9], nodes[10])
    gr.create_edge(nodes[10], nodes[11])
    gr.create_edge(gr.flip(nodes[11]), gr.flip(nodes[1])) # looping chain
    gr.create_edge(nodes[11], nodes[12])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)
    vg_process(filename)
    vg_view(filename, 13)
    print("looping_chain graph created")

# Example call to generate all
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Graph Generator")
    parser.add_argument("-o", "--output", type=str, default="graph_test", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    def wrap_filename(name):
        return os.path.join(args.output, name)

    # Call desired graph creation functions using wrap_filename()
    create_simple_snp_graph(wrap_filename("simple_snp.hg"))
    create_3th_snp_graph(wrap_filename("3th_snp.hg"))
    create_insert_snp_graph(wrap_filename("insert_snp.hg"))
    create_deletion_snp_graph(wrap_filename("deletion_snp.hg"))
    create_insert_deletion_graph(wrap_filename("insert_deletion.hg"))
    create_loop_simple_graph(wrap_filename("loop_simple.hg"))
    create_loop_graph(wrap_filename("loop.hg"))
    create_loop_double_graph(wrap_filename("loop_double.hg"))
    create_loop_plus_graph(wrap_filename("loop_plus.hg"))
    create_linear_path(wrap_filename("linear.hg"))
    create_snp_and_nested_snp_graph(wrap_filename("snp_and_nested_snp.hg"))
    create_complex_ins_graph(wrap_filename("complex_ins.hg"))
    create_4th_graph(wrap_filename("4th.hg"))
    create_repetition_graph(wrap_filename("repetition.hg"))
    create_large_del_graph(wrap_filename("large_del.hg"))
    create_inversion_graph(wrap_filename("inversion.hg"))
    create_nested_plus_graph(wrap_filename("nested_plus.hg"))
    create_jean_graph(wrap_filename("jean.hg"))
    create_multicomponent_chain_graph(wrap_filename("multicomponent_chain.hg"))
    create_looping_chain_graph(wrap_filename("looping_chain.hg"))

# python3 scripts/create_graph.py

# vg convert simple_snp.hg > simple_snp.pg
# vg index simple_snp.pg -j simple_snp.dist

# vg find -x 4th.pg -r 1:6 -c 10 | vg view -dp - | dot -Tsvg -o 4th.svg
