from sage.groups.free_group import is_FreeGroup
from sage.rings.infinity import PlusInfinity
from itertools import product
import os
import pickle
from tqdm import tqdm
from free_groups.core_graph import CoreGraph, create_covers_DAG
from free_groups.free_group_automorphism import get_whitehead_move, get_whitehead_move_of_cut, FreeGroupAutomorphism
from networkx import Graph, connected_components
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

# === Generating All Automrophism Classes ===


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


def get_all_aut_classes(F, length, dirname="aut_classes_cache/", verbose=True):
    """
    Get all automorphism classes of words
    in F_r with bounded length.
    Caches the result to dirname.
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()

    cache_dir = os.fsencode(dirname)
    cache_file = os.fsencode(f"r{r}-len{length}.pkl")
    if os.path.exists(cache_dir + cache_file):
        aut_classes = pickle.load(open(cache_dir + cache_file, 'rb'))
        return [set([F(w.Tietze()) for w in cls]) for cls in aut_classes]
    # Maybe we computed something bigger before
    for file in os.listdir(cache_dir):
        filename = os.fsdecode(file)
        r_str, len_str = filename.split(".")[0].split("-")
        r_cached = int(r_str[1:])
        len_cached = int(len_str[3:])
        if r_cached > r and len_cached > length:
            cached_aut_classes = pickle.load(open(cache_dir + file, 'rb'))
            aut_classes = []
            for cls in cached_aut_classes:
                word = cls.pop()
                word_rep = word.Tietze()
                if len(word_rep) == 0:
                    aut_classes.append(set([F(1)]))
                elif len(word_rep) < length and max(set([abs(x) for x in word_rep])) < r:
                    cls.add(word)
                    aut_classes.append(set([F(w.Tietze()) for w in cls]))
            pickle.dump(aut_classes, open(f"aut_classes_cache/r{r}-len{length}.pkl", 'wb'))
            return aut_classes

    letters = list(range(1, r + 1)) + list(range(-r, 0))

    minimal_words = set()
    all_words = set()  # To avoid stuff like (1,-1,2,3) and (2,-2,2,3)

    # We only consider words that are in "canonical order"
    # We also assume we only check tuple starting with a
    # (They can still be a*a^-1*b*...)
    tuples = product(letters, repeat=length - 1)
    if verbose:
        tuples = tqdm(tuples)
    for tup in tuples:
        tup = [1] + list(tup)
        word = F(tup)
        if word not in all_words and word == canonical_letter_permute_form(F, word):
            all_words.add(word)
            minimal_words.add(canonical_letter_permute_form(F, minimize(F, word)))
    if length > 0:
        # Due to cancellations (e.g. (1,-1, ...)), we only
        # need to check words of length N, N - 1.
        tuples = product(letters, repeat=length - 1)
        if verbose:
            tuples = tqdm(tuples)
        for tup in tuples:
            word = F(tup)
            if word == canonical_letter_permute_form(F, word):  # We only consider canonized orders
                minimal_words.add(canonical_letter_permute_form(F, minimize(F, word)))
    if verbose:
        print(f"Finished minimizing letters, found {len(minimal_words)} minimal words.")
        print("Creating the Whitehead moves graph on minimal words.")

    G = Graph()
    G.add_nodes_from(minimal_words)
    for word in minimal_words:
        nbrs = get_minword_wh_nbrs(F, word)
        assert(len(nbrs[0].Tietze()) == len(word.Tietze())), "Something's wrong"
        for nbr in nbrs:
            nbr = canonical_letter_permute_form(F, nbr)
            assert nbr in minimal_words, f"Found a word ({nbr}) not in minimal_words"
            G.add_edge(word, nbr)
    aut_classes = list(connected_components(G))
    pickle.dump(aut_classes, open(f"aut_classes_cache/r{r}-len{length}.pkl", 'wb'))
    return aut_classes

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


def get_primitivity_rank_of_subgroup(F, group_gens, verbose=False):
    return get_critical_groups_of_subgroup(F, group_gens, verbose)[0]


def get_primitivity_rank(F, word, verbose=False):
    return get_critical_groups(F, word, verbose)[0]
