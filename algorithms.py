from sage.groups.free_group import is_FreeGroup
from sage.rings.infinity import PlusInfinity
from itertools import product
from free_groups.core_graph import CoreGraph, create_covers_DAG
from free_groups.free_group_automorphism import get_whitehead_move, get_whitehead_move_of_cut, FreeGroupAutomorphism
from networkx import Graph, set_edge_attributes
from networkx.algorithms.flow import minimum_cut

infinity = PlusInfinity()

# === Word Minimization ===


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

    word = shortest_cyclic_shift(F, word)
    new_word = get_best_whitehead_nbrs(F, word)[0]

    if len(new_word.Tietze()) < len(word.Tietze()):
        return minimize(F, new_word)
    else:
        return word


def get_best_whitehead_nbrs(F, word):
    """
    Find the shortest words that can be
    obtained by applying a single Whitehead
    move to a word.
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()
    assert word in F, "word must be in F"

    letters = list(range(1, r + 1)) + list(range(-r, 0))

    G = get_whitehead_graph(F, word)

    best_change = 0
    best_moves = []
    for v in letters:
        cap, cut = minimum_cut(G, v, -v)
        change = cap - sum([G[v][u]['capacity'] for u in G[v]])
        if change < best_change:
            best_change = change
            best_moves = [(v, cut[0])]
        elif change == best_change:
            best_moves.append((v, cut[0]))

    assert best_change <= 0, "Something is wrong"

    new_words = []
    for v, cut in best_moves:
        phi = get_whitehead_move_of_cut(F, v, cut)
        new_words.append(shortest_cyclic_shift(F, phi(word)))

    assert all([len(new_word.Tietze()) == len(word.Tietze()) + best_change
                for new_word in new_words]), "Something is wrong"

    return new_words


def canonical_letter_permute_form(F, word):
    """
    Apply a permutation of the letters of F so
    that the different letters appearing in the
    word appear in order (1,2,3,...).
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()
    assert word in F, "word must be in F"

    word_rep = word.Tietze()
    letter_order = []
    sgns = []
    for letter in word_rep:
        if abs(letter) not in letter_order:
            letter_order.append(abs(letter))
            sgns.append((1 if letter > 0 else -1))
    gen_imgs = [i for i in range(1, r + 1)]
    for i, letter in enumerate(letter_order):
        gen_imgs[letter - 1] = (sgns[i] * (i + 1),)
    for i, v in enumerate([v for v in range(1, r + 1) if v not in letter_order]):
        gen_imgs[v - 1] = (len(letter_order) + i + 1,)
    phi = FreeGroupAutomorphism(F, gen_imgs)
    return phi(word)

# === Primitivity Rank and Algebraic Extensions ===


def get_algebraic_extensions(F, group_gens, verbose=False):
    """
    Get the list of all algebraic extensions of G.
    (The algebraic extension are exactly the X-covers of G
    which are not an immediate quotient of any cover of
    smaller rank.)
    """
    cg = CoreGraph(F, group_gens)
    DAG = create_covers_DAG(cg, verbose)
    alg_exts = []
    for ch in DAG:
        if ch == cg.subgroup_gens:
            continue
        if not any([len(ch) > len(v) for v in DAG.predecessors(ch)]):
            alg_exts.append((ch, len(ch)))
    return alg_exts


def get_critical_groups_of_subgroup(F, group_gens, verbose=False):
    """
    Return the primitivity rank and the critical groups of
    group_gens.
    """
    alg_exts = get_algebraic_extensions(F, group_gens, verbose)
    if len(alg_exts) == 0:
        return infinity, []
    min_rank = min([rk for gens, rk in alg_exts])
    critical_groups = [gens for gens, rk in alg_exts if rk == min_rank]
    assert min_rank <= F.rank(), "Something is wrong."
    return min_rank, critical_groups


def get_critical_groups(F, word, verbose=False):
    return get_critical_groups_of_subgroup(F, [minimize(F, word)], verbose)