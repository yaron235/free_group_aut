from sage.groups.free_group import is_FreeGroup


class FreeGroupAutomorphism(object):
    """
    A class for representing automoprhism of free groups.
    """

    def __init__(self, F, gen_imgs):
        """
        F is the free groups.
        gen_imgs are the images of the generators x_1,...,x_r.
        """
        assert is_FreeGroup(F), "F must be a free group!"
        self.F = F
        self.rank = F.rank()
        assert len(gen_imgs) == self.rank, "The number of gen_imgs must be equal to the rank of F"
        try:
            self.gen_imgs = [F(img) for img in gen_imgs]
        except (ValueError, TypeError):
            raise ValueError("Could not convert gen_imgs to elements of F")
        self.inverse_gen_imgs = [img ** (-1) for img in self.gen_imgs]

    def __call__(self, word):
        """
        Apply the autuomorphism on a word.
        """
        assert word in self.F, "word must be in F."
        result = self.F(1)
        for letter in word.Tietze():
            if letter > 0:
                result *= self.gen_imgs[letter - 1]
            else:
                result *= self.inverse_gen_imgs[- letter - 1]
        return result

    def __repr__(self):
        return f"Automorphism defined by {self.gen_imgs}"


def get_whitehead_move(F, v, choices):
    """
    Generate the whitehead move corresponding
    to v with the given choices of (s,t).
    Recall that the whitehead move is defined by
    v->v, x_i-> v^s * x * v^-t

    F - a free group.
    v - a letter of F.
    choices - a list of tuples of pair of binary numbers, of length
    rank(F) - 1.
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()
    assert v in range(1, r + 1) or v in range(-r, 0), "v must be a valid letter of F \
                                            (a number in [-r,...,-1,1,...,r])"
    assert len(choices) == r - 1, "choices must be of length rank(F) - 1"

    gen_imgs = []
    for x in range(1, abs(v)):
        s, t = choices[x - 1]
        gen_imgs.append([c for c in [s * v, x, -t * v] if c])
    gen_imgs.append([abs(v)])
    for x in range(abs(v) + 1, r + 1):
        s, t = choices[x - 2]
        gen_imgs.append([c for c in [s * v, x, -t * v] if c])
    return FreeGroupAutomorphism(F, gen_imgs)


def get_whitehead_move_of_cut(F, v, Y):
    """
    Generate the whitehead move corresponding
    to v with and the cut Y.
    Recall that the whitehead move is defined by
    x -> v^(-1 if -x in Y else 0) * x * v^(1 if x in Y else 0)
    v -> v

    F - a free group.
    v - a letter of F.
    Y - a set containing v but not -v.
    """
    assert is_FreeGroup(F), "F must be a free group"
    r = F.rank()
    assert v in range(1, r + 1) or v in range(-r, 0), "v must be a valid letter of F \
                                            (a number in [-r,...,-1,1,...,r])"
    assert v in Y and -v not in Y, "The cut must seperate v and -v"

    gen_imgs = []
    for x in range(1, r + 1):
        img = []
        if -x in Y:
            img.append(-v)
        img.append(x)
        if x in Y:
            img.append(v)
        gen_imgs.append(img)
    gen_imgs[abs(v) - 1] = [abs(v)]

    return FreeGroupAutomorphism(F, gen_imgs)
