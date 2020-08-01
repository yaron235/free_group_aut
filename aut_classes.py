from sage.groups.free_group import is_FreeGroup
from itertools import product
import os
import pickle
from tqdm import tqdm
from free_group_aut.algorithms.minimize import minimize
from free_group_aut.algorithms.whitehead_algorithm import \
        get_minword_wh_nbrs, canonical_letter_permute_form
from networkx import Graph, connected_components


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
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
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
        print(f"Minimizing all words in {F} of length <={length}")
        print(f"{2*(len(letters)) ** (length - 1)} words to minimize.")
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
