# free_group_aut
Implementation of various algorithms on free groups. 

Contains:

- Whitehead's minimization algorithm
- Whitehead's algorithm for checking if two words are automorphic
- Computation of primitivity rank
- Computation of algebraic extensions
- ...

### Requirements

This module requires sage and the python package networkx (can be installed using ```pip install networkx```).

### Installation

Copy the directory to your project, and import using

```python
import free_group_aut
```

### Usage Examples

```python
>> import free_group_aut
>> F.<a, b, c> = FreeGroup(3)
>> # Minimization - finding minimimal representative
>> free_group_aut.minimize(F, a)
a
>> free_group_aut.minimize(F, a*b^50)
a
>> free_group_aut.minimize(F, a*b*c*a^-1*b^-1*c^-1)
b*c*b^-1*c^-1
>> # Check if two words are automorphic
>> free_group_aut.are_automorphic(F, a^2*b^2, a*b*a^-1*b^-1)
False
>> # Calculate primitivity rank
>> free_group_aut.get_primitivity_rank(F, a^2*b^2)
2
>> free_group_aut.get_critical_group(F, a^2*b^2)
(2, [(a, b)])
>> free_group_aut.get_algebraic_extensions(F, [a^2, b^3])
[((a, b^3), 2), ((b, a^2), 2), ((a, b), 2)]
>> # it returned a list of all the algebraic extensions, given as (gens, rank).
>> #
>> # Get all automorphism classes of bounded length
>> free_groups.get_all_aut_classes(F, 4)
[{a^2*b^2, a*b*a^-1*b, a*b*a*b^-1, a*b^2*a}, {a^3}, {a}, {a^4}, {a*b*a^-1*b^-1}, {1}, {a^2}]
>> # returned a list of sets of minimal representativesfor automorphism classes.
```



