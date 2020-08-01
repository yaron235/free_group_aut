from sage.rings.infinity import PlusInfinity
from free_group_aut.algorithms.core_graph import CoreGraph, create_covers_DAG
from free_group_aut.algorithms.minimize import minimize

infinity = PlusInfinity()


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
