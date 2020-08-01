from sage.groups.free_group import is_FreeGroup
from sage.rings.infinity import PlusInfinity
from itertools import product
from free_group_aut.algorithms.whitehead_moves import get_whitehead_move, get_whitehead_move_of_cut
from networkx import Graph
from networkx.algorithms.flow import minimum_cut

infinity = PlusInfinity()


def minimize_naive(F, word):
    """
    Find a minimal automorphic representative for a word in F.
    Uses Whitehead minimization.

    This is a naive implementation, going over all possible moves
    instead of finding the most efficient one.
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()
    assert word in F, "word must be in F"

    letters = list(range(1, r + 1)) + list(range(-r, 0))

    curr_word = word
    curr_len = len(word.Tietze())
    for v in letters:
        for choices in product([(0, 0), (0, 1), (1, 0), (1, 1)], repeat=r - 1):
            phi = get_whitehead_move(F, v, choices)
            new_word = phi(word)
            new_len = len(new_word.Tietze())
            if new_len < curr_len:
                curr_word = new_word
                curr_len = new_len
    if curr_len < len(word.Tietze()):
        return minimize_naive(F, curr_word)
    else:
        return word


def minimize(F, word):
    """
    Find the minimal form of a word.
    A minimal form is a representative of (Aut(F)w)
    of shortest length.

    This is a smart implementation - for each letter
    v, we find the best v-cuts.
    """
    assert is_FreeGroup(F), "F must be a free group"
    assert word in F, "word must be in F"

    curr_word = word
    new_word = minimize_one_step(F, curr_word)
    while len(new_word.Tietze()) < len(curr_word.Tietze()):
        curr_word, new_word = new_word, minimize_one_step(F, new_word)

    return curr_word


def minimize_one_step(F, word):
    """
    Find the shortest word that can be
    obtained by applying a single Whitehead
    move to a word.
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()
    assert word in F, "word must be in F"

    word = shortest_cyclic_shift(F, word)

    letters = list(range(1, r + 1)) + list(range(-r, 0))

    G = get_whitehead_graph(F, word)

    best_change = 0
    best_move = tuple()
    for v in letters:
        cap, cut = minimum_cut(G, v, -v)
        change = cap - sum([G[v][u]['capacity'] for u in G[v]])
        if change < best_change:
            best_change = change
            best_move = (v, cut[0])

    assert best_change <= 0, "Something is wrong"

    if best_change == 0:
        return word

    v, cut = best_move
    phi = get_whitehead_move_of_cut(F, v, cut)
    new_word = shortest_cyclic_shift(F, phi(word))
    assert len(new_word.Tietze()) == len(word.Tietze()) + best_change, "Something is wrong"

    return new_word


def get_whitehead_graph(F, word):
    """
    Construct the Whitehead graph of the word.
    This is the graph whose vertices are the
    letters of F, and for any (xy) appearing in
    the word there's an edge between x and y^-1.
    """
    assert is_FreeGroup(F), "F must be a free group"
    assert word in F, "word must be in F"
    r = F.rank()
    G = Graph()
    G.add_nodes_from(list(range(1, r + 1)) + list(range(-r, 0)))
    word_repr = word.Tietze()
    edges = [(word_repr[i], -word_repr[i + 1]) for i in range(-1, len(word_repr) - 1)]
    for u, v in edges:
        if not G.has_edge(u, v):
            G.add_edge(u, v, capacity=1)
        else:
            G[u][v]['capacity'] += 1
    return G


def shortest_cyclic_shift(F, word):
    """
    Find the shortest cyclic shift of a word.
    """
    word_rep = word.Tietze()
    word_len = len(word_rep)
    if word_len == 0:
        return word
    i = 0
    while word_rep[i] == -word_rep[-1 - i]:
        i += 1
    if word_len - i < i - 1:  # word is x^{-1} x, shouldn't happen
        return F(1)
    else:
        return F(word_rep[i:word_len - i])
