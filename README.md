# Partially-Balanced-Latin-Rectangles

This was an undergraduate research project conducted under Professor Sergey Bereg at the University of Texas at Dallas. The project attempts to expand on the work of other authors and identify more applicable balance criteria for Latin Rectangles. It also attempts to find Latin Rectangles with minimal imbalance. 

In a k x n Latin rectangle, each cell has a number in {1, 2, 3… n} and no number is repeated in every row and in every column. Consider the following example:

1, 2, 3, 4
2, 3, 4, 1
3, 4, 1, 2

Latin rectangles have unique properties, making them extremely applicable in many fields that rely on experimental design. They have been found to be particularly useful when time is a confounding or limiting factor in an experiment. Latin squares have proven popular enough that MEDLINE, a bibliographic database, found over 4,000 publications with the phrase “Latin square” in important locations, such as the title or abstract [1]. 

Arbitrary Latin squares may be limited in usability, but additional restraints may make them more suited to particular purposes. The nature of these restrictions varies substantially across fields. 

### Imbalance Functions

Many have attempted to use various criteria in order to define some standard of ‘balance’ for Latin rectangles, as well as to define an imbalance function if the rectangle is not itself balanced. 

However, the purpose of balance criteria varies drastically between authors, so different criteria are used. We have compiled various imbalance functions and proposed some of our own. 

Define d_i(u, v) to be the positive distance between u and v in row i.

For some arbitrary function g(x), let s_g(u, v) be the sum of all g(d_i(u, v)) across all rows. 

Define s(u, v) to be the special case of s_g(u, v) where g(x)=x. We call this the sum-imbalance function, introduced by Diaz et al. [2].

A Latin rectangle is g(x)-balanced if for all distinct u, v ∈ {1, 2, 3… n}, s_g(u, v) is constant. 

The g(x)-imbalance function returns the sum of the absolute differences between s_g(u, v) and the average value of s_g(u, v) for all distinct choices of u and v. The same follows for other various imbalance functions (which are too long to list here, can be described by just their name alone, with implementations found in main.py).

Finally, a Latin rectangle is order-balanced if for all u, v ∈ {1, 2, 3… n} such that u≠v, |x-y|≤1, where x is the number of rows where u precedes v and y is the number of rows where v precedes u.

### Algorithms

A modified version of the Random-restart hill climbing algorithm proposed by Diaz et al. is used as the primary algorithm to find Latin rectangles with low imbalance [2]. Given a Latin rectangle of size (k, n), the algorithm attempts to minimize a given imbalance function across the search space of Latin rectangles. Since the number of Latin squares grows at a substantially larger rate than an exponential function, a brute-force search is impossible.

Therefore, a greedy optimization algorithm is utilized, which randomly selects a pair of columns in a Latin square to swap. Note that swapping preserves the Latin square property. The cyclic Latin square is the starting point for the algorithm, and the top k rows are used to form the Latin rectangle. If swapping columns results in a decrease in imbalance, the swap occurs; otherwise, it doesn't.
 
After a predetermined number of swaps, a new Latin square is created using a modified Jacobson-Matthews algorithm. However, only a few iterations of Jacobson-Matthews are used. The Jacobson-Matthews approach is crucial to confirm that every iteration produces a Latin square and that every possible Latin square can be generated.

Once enough iterations are completed, the process is restarted from the initial cyclic Latin square. This process is repeated for every imbalance function and Latin rectangle size, and the best results are saved in .csv files.

Most imbalance functions examine every pair of elements in each row. The time complexity for computing the imbalance of a k x n Latin rectangle is then O(k\*n^2). Switching two columns can be accomplished in O(1) time. Although the Jacobson-Matthews algorithm can run in O(1) time, I opted for a simpler O(n) implementation since the overall time complexity would still be O(k*n^2).

The pseudocode for computing imbalance and generating low-imbalance Latin rectangles are presented below:

##### Sum imbalance function, using the function g(x) as a parameter

for all (u, v):
	s_g(u, v) = 0
for all (i): # loops through all rows
s_g(u, v) += g(d_i(u, v))
     total_imbalance += abs(s_g(u, v)-average_imbalance)
return total_imbalance

##### Generating low-imbalance LRs

for all (k, n):
	bestImbalance = 1e10 (or some other high constant)
	repeat [a] times:
	L = cyclicLS(n)
repeat [b] times:
	L = JacobsonMatthews(L)
	repeat [c] times:
		shuffled = L.shuffleColumns()
		if imbalance(shuffled) < imbalance(L):
			L = shuffled
		if imbalance(L) < bestImbalance:
			bestImbalance = imbalance(L)
			bestRectangle = L
return bestRectangle

### Results

For each proposed imbalance function, all Latin rectangle sizes up to n=k=12 were run through the earlier algorithm. The lowest _found_ sum imbalances were recorded as follows:

| SUM | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
| - | - | - | - | - | - | - | - | - | - | - | - | 
| 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3 | 1.33 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 4 | 2.67 | 4 | 5.33 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5 | 8 | 6 | 8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 6 | 16 | 12 | 13.33 | 16 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 7 | 28 | 24 | 25.33 | 24.67 | 26 | 18.67 | 0 | 0 | 0 | 0 | 0 |
| 8 | 40 | 40 | 42 | 46 | 46 | 44 | 40 | 0 | 0 | 0 | 0 |
| 9 | 65.33 | 62 | 70.67 | 73.33 | 66 | 86.67 | 77.33 | 36 | 0 | 0 | 0 |
| 10 | 94 | 96 | 99.33 | 110.67 | 112 | 117.33 | 122.67 | 88 | 128 | 0 | 0 |
| 11 | 126 | 130 | 140 | 156 | 164 | 152 | 160 | 190 | 182 | 132 | 0 |
| 12 | 170.67 | 180 | 200 | 201.33 | 210 | 221.33 | 229.33 | 226 | 228 | 186 | 180 |

(Other tables can be found by adjusting main.py accordingly.)

There are a few interesting features of the results. Noting that each table is only an upper bound, it appears that if we were to fix n and observe k, the imbalance function is strictly increasing for any choice of the function g(x). In particular, the sum-imbalance function appears to model a quadratic function. 

As sum-imbalanced Latin rectangles were found that had imbalances equal or close to to prior results, it is likely that the other functions also found useful bounds [2].

### Applications

Latin squares have many applications in experimental design, including controlling for the impact of individual traits on experiment results. For example, in a study testing four prescriptions on four different individuals, a Latin square design can help control for the negative impact of interactions between prescriptions if given in the same order, as each combination is covered and randomized. Balance functions would ensure that every pair of treatments is separated by close-to-equal amounts. The specific balance function used depends on the situation.

The sum-imbalance function, as proposed by Diaz et al., does have merits and practical uses (such as agronomic experiments). However, the imbalance function places too large of a weight on distance, although many phenomena (like contamination or drug residue) have much stronger effects at closer proximity.

The inverse-imbalance and inverse-square-imbalance functions appear to do a better job at modeling contamination emitted from a source, as they assign higher values to closer distances. The inverse-square-imbalance function may return higher quality results in three dimensions due to the inverse-square-law.

The exponential imbalance function can be used for any situation in which the quantity of an item tested shows exponential decay. This could be useful for drug trials, as many substances, such as alcohol, have known half-lives, or other biological or chemical phenomena. 

Order-balanced Latin squares also have substantial use, especially when it comes to cause-and-effect relationships where distance has no effect. Because their existence has been proven for every Latin square size, models requiring any Latin square size can utilize them effectively. 

### References

[1] Richardson, John T. E. (2018). The use of Latin-square designs in educational and psychological research. Educational Research Review, 24 pp. 84–97.
[2] Diaz, Mateo, et al. In Search of Balance: The Challenge of Generating Balanced Latin Rectangles. Center for Applied Mathematics, Cornell University, 2017, https://www.cs.cornell.edu/gomes/pdf/2017_diaz_cpaior_balance.pdf. 
