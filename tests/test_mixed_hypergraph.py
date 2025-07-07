from os import remove

from halp.undirected_hypergraph import UndirectedHypergraph
from halp.directed_hypergraph import DirectedHypergraph
from halp.mixed_hypergraph import MixedHypergraph

# Since most of the tests for basic, single-direction usage of hypergraphs
# reutilize the tests in test_undirected_hypergraph and test_directed_hypergraph,
# the tests here specifically focus on features pertaining to mixed hypergraphs.
 
def test_direction_preservation():
    graph = MixedHypergraph()
    graph.add_nodes(['A', 'B', 'C', 'D', 'E', 'F'])
    graph.add_undirected_hyperedge(['A', 'B', 'C'])
    graph.add_directed_hyperedge(['C', 'F'], ['D', 'E'])

    directed_edge = graph.get_directed_hyperedge_id(['C', 'F'], ['D', 'E'])
    assert graph.is_hyperedge_id_directed(directed_edge)

    undirected_edge = graph.get_undirected_hyperedge_id(['A', 'B', 'C'])
    assert graph.is_hyperedge_id_undirected(undirected_edge)

def test_underlying_graph():
    graph = MixedHypergraph()
    graph.add_nodes(['A', 'B', 'C', 'D'])
    graph.add_undirected_hyperedge(['A', 'B', 'C', 'D'])
    graph.add_directed_hyperedge(['A', 'B'], ['C'])
    graph.add_directed_hyperedge(['B'], ['C', 'D'])

    undirected = graph.underlying_undirected_hypergraph()
    assert undirected.has_hyperedge(['A', 'B', 'C', 'D'])
    assert len(undirected.get_hyperedge_id_set()) == 1

    directed = graph.underlying_directed_hypergraph()
    assert directed.has_hyperedge(['A', 'B'], ['C'])
    assert directed.has_hyperedge(['B'], ['C', 'D'])
    assert len(directed.get_hyperedge_id_set()) == 2

def test_extend():
    H = MixedHypergraph()

    H_u = UndirectedHypergraph()
    H_u.add_nodes(['A', 'B', 'C'])
    H_u.add_hyperedge(['A', 'B'])
    H_u.add_hyperedge(['B', 'C'])
    H.extend_undirected_hypergraph(H_u)
    assert len(H.get_hyperedge_id_set()) == 2
    assert len(H.get_node_set()) == 3

    H_d = DirectedHypergraph()
    H_d.add_nodes(['D', 'E', 'F'])
    H_d.add_hyperedge(['D'], ['E', 'F'])
    H_d.add_hyperedge(['D', 'E'], ['F'])
    H_d.add_hyperedge(['A'], ['D'])
    H.extend_directed_hypergraph(H_d)
    assert len(H.get_hyperedge_id_set()) == 5
    assert len(H.get_node_set()) == 6

    H.get_directed_hyperedge_id(['A'], ['D'])

def test_read_invalid_basic():
    invalid_H = MixedHypergraph()
    try:
        invalid_H.read("tests/data/invalid_mixed_hypergraph.txt")
        assert False
    except IOError:
        pass
    except BaseException as e:
        assert False, e
    
def test_read_invalid_direction():
    invalid_H = MixedHypergraph()
    try:
        invalid_H.read("tests/data/invalid_mixed_hypergraph_direction.txt")
        assert False
    except IOError:
        pass
    except BaseException as e:
        assert False, e

def test_overlapping():
    H = MixedHypergraph()
    H.add_undirected_hyperedge(['A', 'B', 'C'])
    H.add_undirected_hyperedge(['A', 'C', 'B'])
    H.add_directed_hyperedges([(['A', 'B'], ['C']),
                               (['A'], ['B', 'C']),
                               (['A'], ['C', 'B'])])
    
    assert len(H.get_hyperedge_id_set()) == 3

def test_get_hyperedge_nodes():
    H = MixedHypergraph()
    id = H.add_directed_hyperedge(['A', 'C'], ['B'])
    assert H.get_hyperedge_nodes(id) == set(['A', 'B', 'C'])

def test_trim_node():
    H = MixedHypergraph()
    H.add_undirected_hyperedge(['A', 'B', 'C'])
    H.add_directed_hyperedge(['A'], ['B', 'C'])
    H.trim_nodes(['C'])

    assert len(H.get_node_set()) == 2
    assert len(H.get_hyperedge_id_set()) == 2

    H.get_undirected_hyperedge_id(['A', 'B'])
    H.get_directed_hyperedge_id(['A'], ['B'])

    H.trim_node('B')

    assert len(H.get_node_set()) == 1
    assert len(H.get_hyperedge_id_set()) == 1

    H.get_undirected_hyperedge_id(['A'])

def test_read_copy():
    H = MixedHypergraph()
    H.read("tests/data/basic_mixed_hypergraph.txt")

    def correct(G):
        # correctness checks
        assert len(H.get_hyperedge_id_set()) == 8
        assert len(H.get_node_set()) == 8

        G.get_undirected_hyperedge_id(['s', 'y', 'x'])
        G.get_undirected_hyperedge_id(['x', 's'])
        G.get_undirected_hyperedge_id(['a', 'u', 't'])

        G.get_directed_hyperedge_id(['s'], ['x'])
    
    correct(H)
    correct(H.copy())

def test_symmetric_image():
    H = MixedHypergraph()
    H.add_undirected_hyperedge(['A', 'B', 'C'])
    H.add_directed_hyperedge(['A'], ['B', 'C'])

    assert len(H.get_hyperedge_id_set()) == 2

    new_H = H.get_symmetric_image()
    assert len(new_H.get_hyperedge_id_set()) == 2

    new_H.get_undirected_hyperedge_id(['A', 'B', 'C'])
    new_H.get_directed_hyperedge_id(['B', 'C'], ['A'])

def test_write():
    H = MixedHypergraph()
    H.add_nodes(['A', 'B', 'C'])

    H.add_directed_hyperedge(['A'], ['B', 'C'])
    H.add_directed_hyperedge(['A', 'B'], ['C'])
    H.add_undirected_hyperedge(['A', 'B'])
    H.add_undirected_hyperedge(['B', 'C'])

    H.write('test_mixed_read_and_write.txt')

    # read the just written file
    new_H = MixedHypergraph()
    new_H.read('test_mixed_read_and_write.txt')

    new_H.get_directed_hyperedge_id(['A'], ['B', 'C'])
    new_H.get_directed_hyperedge_id(['A', 'B'], ['C'])
    new_H.get_undirected_hyperedge_id(['A', 'B'])
    new_H.get_undirected_hyperedge_id(['B', 'C'])

    remove('test_mixed_read_and_write.txt')
