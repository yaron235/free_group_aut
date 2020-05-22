from sage.groups.free_group import is_FreeGroup
from networkx import MultiDiGraph, DiGraph, Graph
from networkx.algorithms import maximum_spanning_tree, shortest_path
from collections import defaultdict
from copy import deepcopy
from itertools import combinations


class CoreGraph(object):
    """
    A class for Stallings' core graphs.
    Object in this class are supposed to be immutable.
    Two core graphs are the same if they have the same generators for the subgroup
    (as obtained by the code's computation).
    """

    def __init__(self, F, arg):
        """
        If arg is a list:
        Generate the core graph corresponding to the group generated by group_gens.
        If arg is a graph:
        Check for validity and do all folds.
        """
        assert is_FreeGroup(F), "F must be a free group"
        self.F = F
        self.F_rank = F.rank()
        self.letter_types = range(1, self.F_rank + 1)
        self.letters = list(range(-self.F_rank, 0)) + list(range(1, self.F_rank + 1))
        # -r, ..., -1, 1, ..., r

        if isinstance(arg, list):
            group_gens = arg
            assert all([gen in F for gen in group_gens]), "The generators must be elements of F."
            self.group_gens = group_gens

            G = MultiDiGraph()
            G.add_node((0,))  # the marked vertex (id)
            for i, gen in enumerate(self.group_gens):
                word = gen.Tietze()
                word_len = len(word)
                # new_nodes = [(i, j) for j in range(1, word_len)]
                # G.add_nodes_from(new_nodes)
                get_node = lambda j: (0,) if (j % word_len == 0) else (i, j)
                for j in range(word_len):
                    G.add_edge(get_node(j), get_node(j + 1), label=word[j])
                    G.add_edge(get_node(j + 1), get_node(j), label=-word[j])

        elif isinstance(arg, MultiDiGraph):
            # We are going to copy the graph, add reverse edges when needed,
            # and sort the edges.
            # The reason we sort the edges is to get a "canonical" version
            # of the object, so subgroup_gens would be the same in different
            # objects with the same graph.
            G = MultiDiGraph()
            G.add_nodes_from(arg.nodes())
            edges = arg.edges(data='label')
            G_edges = [e for e in edges]
            assert len(edges) == len(arg.edges()), "Every edge has to be labelled."
            for src, dst, letter in edges:
                assert letter in self.letters, \
                    f"The edge betwen {src} and {dst} has an invalid label"
                if (dst, src, -letter) not in G.edges(data='label'):
                    G_edges.append((dst, src, -letter))
            G.add_weighted_edges_from(sorted(G_edges), weight='label')

        else:
            raise ValueError("arg must be a list of words or a MultiDiGraph.")

        self.G = do_all_folds(G)

        # The subgraph of positive edges
        G_pos = MultiDiGraph()
        G_pos.add_edges_from([e for e in self.G.edges(data=True) if e[2]['label'] > 0])
        self.G_pos = G_pos

        self.subgroup_gens = tuple(sorted(self.get_subgroup()))

    def get_subgroup(self):
        """
        Generate the subgroup from the graph as pi_1^X(G).
        """
        G_pos = self.G_pos
        T = maximum_spanning_tree(G_pos.to_undirected())
        # Because in list(T.edges), src and dst may be messed up,
        # we have to check each one
        T_edges = []
        for src, dst, letter in T.edges(data='label'):
            for key, eattr in self.G[src][dst].items():
                if eattr['label'] == letter:
                    T_edges.append((src, dst, letter))
                    break
                elif eattr['label'] == -letter:
                    T_edges.append((dst, src, letter))
                    break

        G_pos_edges = G_pos.edges(data='label')

        gens = []

        edges_to_add = [edge for edge in G_pos_edges if edge not in T_edges]

        origin = (0, )
        for src, dst, letter in edges_to_add:
            gen = self.F(1)
            gen *= self.read_path(shortest_path(T, (origin), src), T_edges)
            gen *= self.F([letter])
            gen *= self.read_path(shortest_path(T, dst, origin), T_edges)
            gens.append(gen)

        assert len(gens) == self.rank(), "Something is wrong."
        return gens

    def read_path(self, path, edge_set):
        """
        Go along path and read the emerging word.
        Each time we choose the edge in edge_set. If there's more than one,
        choose arbitrarily.
        """
        result = []
        for i in range(len(path) - 1):
            src = path[i]
            dst = path[i + 1]
            found_edge = False
            for key, value in self.G[src][dst].items():
                letter = value['label']
                if (src, dst, letter) in edge_set or (dst, src, -letter) in edge_set:
                    result.append(letter)
                    found_edge = True
                    break
            if not found_edge:
                raise ValueError(f"No edge between {src} and {dst} is found in the edge set.")
        return self.F(result)

    def rank(self):
        """
        Return the rank of the subgroup corresponding to the core graph.
        """
        return 1 - len(self.G.nodes()) + len(self.G.edges()) / 2

    def get_whitehead_hypergraph(self):
        """
        Return a representation of the Whitehead hypergraph
        of the core graph - a bipartite graph, with one side
        representing the vertices and the other side the
        hyperedges.
        Recall that the Whitehead hypergraph is constructed
        as follows:
        Its vertices are the letters of F and their inverses
        (with respect to a fixed basis), and to every vertex
        in the core graph corresponds an hyperedge of all
        the labels of edges ending on it.
        """
        W = Graph()
        W.add_nodes_from(self.letters)
        for v in self.G.nodes():
            new_edges = [(letter, v) for src, dst, letter in self.G.in_edges(v, data='label')]
            W.add_edges_from(new_edges)
        return W

    def __repr__(self):
        return f"Core graph of {self.subgroup_gens}"

    def __str__(self):
        return f"Core graph on nodes {list(self.G.nodes())}, with edges \n" + \
            "\n".join([f"{o}->{i} : {w}" for o, i, w in self.G_pos.edges(data='label')])

    def __eq__(self, other):
        if not isinstance(other, CoreGraph):
            raise TypeError("other is not a CoreGraph object")

        return self.get_subgroup == other.get_subgroup


def remove_edge(G, node1, node2, letter):
    """
    Remove an edge from {node1} to {node2} labelled {letter} in {G}.
    If there is more than one option, choice is arbitrary.
    """
    is_found = False
    for key, eattr in G[node1][node2].items():
        if eattr['label'] == letter:
            G.remove_edge(node1, node2, key)
            is_found = True
            break
    if not is_found:
        raise ValueError(f"No edge labelled {letter} from {node1} to {node2}")


def remove_edge_and_inverse(G, node1, node2, letter):
    """
    Remove an edge from {node1} to {node2} labelled {letter}, and the opposite edge,
    labelled {-letter}.
    """
    remove_edge(G, node1, node2, letter)
    remove_edge(G, node2, node1, -letter)


def glue(G, node1, node2):
    """
    Glue the nodes node1, node2.
    We assume that the only edge data is 'label', and that if an edge exists,
    then its inverse exists but with -label.
    Thus we go over all edges going out of node2, and add them and their reverses to node1.
    """
    assert isinstance(G, MultiDiGraph), "G must be a graph"
    assert all([v in G.nodes() for v in [node1, node2]]), "The nodes must be in the graph"

    G_new = deepcopy(G)

    for tmp, nbr, letter in G.edges(node2, data='label'):
        if nbr != node2:  # Not a self loop
            G_new.add_edge(node1, nbr, label=letter)
            G_new.add_edge(nbr, node1, label=-letter)
        elif letter > 0:  # To avoid adding self-loops twice
            G_new.add_edge(node1, node1, label=letter)
            G_new.add_edge(node1, node1, label=-letter)

    G_new.remove_node(node2)
    return G_new


def do_all_folds(G, debug=False):
    """
    Apply all possible folds in a graph. A possible fold is a set of two edges labelled the same
    from the same vertix, and the actual folding is gluing the edges.
    If debug is True, we print some debugging information.
    """
    if debug:
        print("======== Entering function")
        out(G)
    is_needed = False
    G_new = deepcopy(G)

    for node in G.nodes():
        if debug:
            print(f"--- Trying node {node}")
        letters_to_edges = defaultdict(list)  # We will collect edges with the same label here
        for tmp, nbr, letter in G.out_edges(node, data='label',):
            letters_to_edges[letter].append(nbr)

        if debug:
            print(list(letters_to_edges.items()))
        for letter, nbrs in letters_to_edges.items():
            if len(nbrs) > 1:  # There is some gluing to do here
                nbrs = sorted(nbrs)
                if debug:
                    print(f"Working on letter {letter}")
                    print(f"Neighbots are {nbrs}")
                nbr0 = nbrs[0]
                for nbr in set(nbrs[1:]):  # We don't need to go over the same nbr twice
                    if nbr != nbr0:
                        G_new = glue(G_new, nbr0, nbr)
                    if debug:
                        print(f"Glued {nbr0} and {nbr}, {nbr} was removed")
                # We're left with len(nbrs) - 1 edges from node to nbr0
                # However, if node is a neighbor of itself, and node > nbr0,
                # we will have only self-loops around nbr0.
                if debug:
                    print("Going to eliminate self-loops/node-nbr0 edges")
                    out(G_new)
                if node in nbrs:
                    for i in range(len(nbrs) - 1):
                        remove_edge_and_inverse(G_new, nbr0, nbr0, letter)
                else:
                    for i in range(len(nbrs) - 1):
                        remove_edge_and_inverse(G_new, node, nbr0, letter)

                return do_all_folds(G_new, debug)
                # We finished with this node, so we start all over

    # When there are no more folding to be done, we will finally exit the loop
    if is_needed:
        return do_all_folds(G_new, debug)
    else:
        return G_new


def create_covers_DAG(cg, verbose=False, debug=False):
    """
    Given a CoreGraph cg, create its DAG of X-covers.
    That is, create the graph of all possible gluings
    of cg, with an arrow from G to H if H is obtained by
    gluing exactly 2 vertices in G (and folding if necessary).
    In such a case, H is called an immediate quotient of G.

    This is a naive implempentation - we go over leaves
    and add immediate quotients, until there's only one
    vertex left.
    This can be improved by remembering all the partitios
    we get while doing the foldings.
    """
    if verbose:
        print("Creating DAG...")
    F = cg.F
    DAG = DiGraph()
    max_weight = cg.G.number_of_nodes()
    DAG.add_node(cg.subgroup_gens, weight=max_weight)
    gens_to_coregraph = {cg.subgroup_gens: cg}
    for curr_weight in range(max_weight, 1, -1):
        curr_nodes = [(gens, weight) for gens, weight in DAG.nodes(data='weight')
                      if weight == curr_weight]
        if verbose:
            print(f"Working on weight {curr_weight}: {len(curr_nodes)} nodes to process")
        for gens in [gens for gens, weight in curr_nodes]:
            ch = gens_to_coregraph[gens]
            for u, v in combinations(ch.G.nodes(), 2):
                if u > v:
                    u, v = v, u
                new_ch = CoreGraph(F, glue(ch.G, u, v))
                new_gens = new_ch.subgroup_gens
                if new_gens not in gens_to_coregraph:
                    gens_to_coregraph[new_gens] = new_ch
                    DAG.add_node(new_gens, weight=new_ch.G.number_of_nodes())
                DAG.add_edge(gens, new_gens, gluing=(u, v))
    if debug:
        return DAG, gens_to_coregraph
    else:
        return DAG


def out(G):
    """
    Pretty output for a labelled MultiDiGraph.
    """
    edges = list(G.edges(data='label'))
    for i, o, w in edges:
        if (w < 0 and (o, i, -w) not in edges) or (w > 0):
            print(f"{i}->{o} : {w}")
