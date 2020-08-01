from sage.groups.free_group import is_FreeGroup
from itertools import product
from free_group_aut.algorithms.minimize import minimize
from free_group_aut.algorithms.whitehead_moves import get_whitehead_move
from networkx import Graph, has_path


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
    phi = F.hom([F(img) for img in gen_imgs])
    return phi(word)


def get_minword_wh_nbrs(F, word):
    """
    Assuming word is (AutF)-minimial, find
    all words that can be obtained from
    word by exactly one Whitehead move.

    The assumption on minimality is not
    necessary, but we just find whitehead
    move that don't change the length.
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()
    assert word in F, "word must be in F"

    letters = list(range(1, r + 1)) + list(range(-r, 0))

    word_len = len(word.Tietze())
    nbrs = set()
    for v in letters:
        for choices in product([(0, 0), (0, 1), (1, 0), (1, 1)], repeat=r - 1):
            phi = get_whitehead_move(F, v, choices)
            nbr = phi(word)
            if len(nbr.Tietze()) == word_len:
                nbrs.add(nbr)

    return list(nbrs)


def are_automorphic(F, w1, w2):
    assert is_FreeGroup(F), "F must be a free group"
    assert w1 in F and w2 in F, "Both words must be in F"
    w1 = canonical_letter_permute_form(F, minimize(F, w1))
    w2 = canonical_letter_permute_form(F, minimize(F, w2))
    if w1 == w2:
        return True
    if len(w1.Tietze()) != len(w2.Tietze()):
        return False
    G = Graph()
    G.add_nodes_from([w1, w2])
    tested_words = set()
    while len(tested_words) < len(G.nodes()):
        edges_to_add = []
        for word in G.nodes():
            if word in tested_words:
                continue
            nbrs = [canonical_letter_permute_form(F, nbr) for nbr in get_minword_wh_nbrs(F, word)]
            edges_to_add += [(word, nbr) for nbr in nbrs]
            tested_words.add(word)
        G.add_edges_from(edges_to_add)
    return has_path(G, w1, w2)
