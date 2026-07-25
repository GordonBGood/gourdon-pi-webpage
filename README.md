# Gourdon Prime Counting Webpage

## Description of Provided Webpage

This repository implements a slightly simplified version of the Gourdon Prime Counting Algorithm as a web page that can be used to count the number of primes to a limit between zero and 9007199254740991 (`2**53-1`).  This upper range is limited to the highest integer number that can be used efficiently in JavaScript as the largest whole integer number that fits into the 53 bit mantissa of JavaScript's 64-bit IEEE floating point numbers since all numeric quantities are represented as such.  In order to represent higher value numbers, one would need to use representations such a "BigInt" that represent 64-bit integers as two 32-bit ones (or more), and manipulating such would be quite slow for the required multiplications and divisions as to required cross 32-bit portion operations.  JavaScript counting operations are already quite slow as JavaScript only supports representing integers limited to 32-bits for many bit operations such as bit shifting, bit and'ing, bit or'ing, etc., and thus JavaScript requires many more operations to emulate 64-bit operations; as well, current versions of JavaScript/ECMASCRIPT do not have access to hardware CPU "pop counting" instructions (including SIMD AVX512 "pop counting" on newer CPU's), so the "pop counting" emulation as must be used is many times slower.  Thus the JavaScript code as used for this webpage is three to five times slower than if the code were translated to a native compiling languages such as C/C++, where multi-threading could also be used to reduce the execution time by a factor of the number of CPU cores.

## The Live GitHub Pages Webpage Link

If you just want to use the prime counting facility provided by this webpage, the live GitHub Pages webpage link is:  https://gordonbgood.github.io/gourdon-pi-webpage/.

## History of Prime Counting Algorithms

Up until about 1808, the best way of counting the number of primes up to a limit was the ancient Sieve of Eratosthenes (SoE) which may have been adapted even by Eratosthenes to treat only odd numbers since "two" is the only even prime.  Carl Friedrich Gauss actually used "page-segmentation" in counting the number of primes within blocks each with a range of 1000 numbers (he also very likely used limited "wheel-factorization" by the only even prime of "two") to calculate the number of primes up to a million (`10**6`), which results were found after his death in 1855 to be entirely accurate at 78,498.

### Legendre's Prime Counting Algorithm

During this time, the first combinational algorithm to determine the number of primes to a limit **without actually considering each of the values within the range up to that limit** using the "inclusion/exclusion" principle (considering combinations of how many primes in the range can be considered to be possible primes versus the number of primes in that range that are impossible to be primes) by Adrien-Marie Legendre was published and in 1808 he published his hand calculated number of primes to a million as 78,526, in error compared to the correct answer given above due to him not using quite enough primes in his calculation.  His Legendre combinational prime counting formula is as follows:

$$
\begin{align*}
\pi(x) &= \phi(x,a) + a - 1 \\
a &= \pi(\sqrt{x}) \\
\phi(x,a) &= {x} - \sum_{i}^{a} {\lfloor \frac{x}{p_i} \rfloor}
                 + \sum_{j>i}^{a} {\lfloor \frac{x}{p_i p_j} \rfloor}
                 - \sum_{k>j>i}^{a} {\lfloor \frac{x}{p_i p_j p_k} \rfloor}
                     \dots + (-1)^a {\lfloor \frac{x}{p_1 p_2 p_3 \dots p_a} \rfloor}  \\
&= \sum_{d \mid P_a} \mu(d) {\lfloor \frac{x}{d} \rfloor} \\
&= \phi(x,a-1) - \phi(x/p_a, a-1) \\
P_a &= \prod_{i=1}^{a} p_i \\
\end{align*} 
$$

Understand that $$P_a$$ is the product of all of the primes up to the square root of the counting limit.  $\mu(d)$ is the Moebius Function for ${d}$ which normally has values of plus or minus one and zero but as zero represents any value that has a factor of a prime squared and $$P_a$$ contains only unique prime factors by definition, a zero result in this case is impossible.  A plus one Moebius value means that the argument has an odd number of prime factors and a minus one value means that the argument has an even number of prime factors, so application here means that the summation adds the floor of the quotients when the number of prime factors is even and subtracts when odd.  Legendre would not have known this function by name as the name is relatively new, but it is clear from the orginal $\phi(x,a)$ expansion that odd numbers of factors are subtracted and even numbers are added.  The notation ${d \mid P_a}$ means that ${d}$ is in the set of all numbers that can be evenly divided into $P_a$, meaning they are all the combinations of products of unique primes up to the square root of ${x}$.  This Moebius Function expressed summation for the summation of the inclusion/exclusion primciple above it is what Legendre would have used to hand calculate the number of primes to a million for a couple of reasons:  1) the idea of function recursion hadn't been invented yet, nor was there a concept of a last-in/first-out stack as required by recursive functions, and 2) this Moebius Function summation can have a much better asymptotic complexity of `O(x**(2/3)/((log x)**2))` than the `O(x/((log x)**2))` that the recursive formulation if combined with "Partial Sieving" or doing processing after each culling pass by each individual base prime. The reason that "Partial Sieving" can exponentially reduce the required number of operations is due to that this formulation initializes to all the values that are products of the primes and gradually reduces them by eliminating the ones that are factors of the base primes in turn starting at the least as compared with starting with all the individual base primes and building up all the combinations of the products (which requires many more recursions) for the recursive formulation.  Seeing practical implementations for each will make this clear.

It is not clear whether Legendre actually used "Partial Sieving" in his hand calculation of the number of primes to a million, although he certainly used some optimizations to reduce them.  If the above formula is used without some terminating conditions, it can't be used practically:  the number of operations implied by the above formula without a terminating condition is immense, even counting just to a million as it is the factorial of the number of primes to the square root of the limit, or 168 for the square root of a million of one thousand (167 for odds-only), which is a huge number with about 300 zeros - more operations than there are atoms in the observable universe!  However, many of these combinations are eliminated just by stopping expanding when the ${d}$ product is greater than the counting limit, in which case the contribution is zero, or even better by stopping when ${x}/{d} <= {lpf}$ where ${lpf}$ is the least prime factor of ${p_i}$ in the expansion expression; for the first condition, the number of division operations are reduced to just `95,733` for the first termination condition and `15,321` for the second, both when using the non "Partial Sieving" implementation of Phi finding the number of primes to a million.  With even the second strong termination, this would have taken Legendre a few years to hand calculate; with "Partial Sieving", the number of operations is reduced to just `2,408` operations to find the number of primes to a million that could have been hand calculation in a few months.  It isn't clear that Legendre actaully used "Partial Sieving", in which case his math error in so many hand calculations could easily have occurred as even with the strong termination condition, there are about 720 combinations of the smallest primes he would have had to process along with the loop over all the multiple by all the primes up to the terminating condition, but all major prime counting algorithms since then, including Meissel's hand calculation of the count of primes to a billion (1,000,000,000), would have used it, with the exception of the first computer adaptation of the algorithms by DH Lehmer, who only used bottom-up recursion with the strong termination condition in calculating Phi.

Rather than confuse the description of the implementation of the algorithms with more words, here are a Python-as-pseudocode implementations of some implementing code for bottom-up recursive Phi calculation; where a "bottom-up" implementation is the superior choice because it doesn't require memoization/caching of intermediate values to avoid exponentially slower execution with increasing counting range.

The recursive formula for Phi is as follows:

$$
\begin{align*}
\phi(x,a) &= \phi(x,a-1) - \phi(x/p_a, a-1) \\
\end{align*} 
$$

A Python implementation of that recursion is as follows:

```python
# given a `prms` list of all the odd primes up to the square root of `x`...
def phi(x, a):
    ```find the Phi value when `x` is the counting limit and\r\n
          `a` is the number of primes to the square root of x.```
    def tinyphi(n): return (x - 1) >> 1
    def lvl(lpfi, lpfisz, m):
        phi = 0
        for pi in range(lpfi, lpfisz):
            p = prms[pi]; nm = m * p
#                if nm > sqrtlmt: return phi # weaker termination condition!
            if p * nm >= limit: return phi + lpfisz - pi # strong termination!
            phi += tinyphi(lmt // nm)
            phi -= lvl(0, pi, nm)
        return phi
    return tinyphi(x) - lvl(0, a - 1, 1) # one less than `a` because odds-only!
```

However, as mentioned, recursion wasn't really a concept in the days of Legendre so he would have implemented it with looping procedurally with what would eventually be a last-in/first-out "stack" replaced by a few memory locations (lines on a piece of paper).  a Python-as-pseudocode implementation with strong termination of looping replacing recursion is as follows - not that this has just as many operations, just doesn't need functions and slacks:

```
TO BE PROVIDED LATER!!!
```

The recursive definition and implementation of Phi is so simple that computer programmers fall into the trap of using it, ignoring the poor execution complexity.  Even in this bottom-up implementation, this expression has an asymptotic complexity of `O(x/((log x)**2))`.  With increasing counting range, this builds much much faster than the non-recursive formula; these formulas use memory of `O(x/(log x))` for the array of primes but require memory use of`O(x**(1/2))` due to the need to SoE sieve to the square root of the counting range in order to find the base primes (although it might be smaller by a constant factor if the sieve is implemented to use just one bit per odd number representation) and increased amount by a constant factor for the non-recursive implementations that require other arrays of this size.

For Legendre's case in using his formula to find the number of primes to a million, the square root of a million is just a thousand, his array(s) would only have needed 500 elements for odds-only from which the 167 odd primes could have been determined, and his small culling base primes would have only gone up to 31 with the values of 3, 5, 7, 11, 13, 17, 19, 23, 29, and 31.  One can play with these algorithms in Python for to see the difference as [the non-recursive Legendre algorithm](https://wandbox.org/permlink/B40zkHgGgn3YBvtT) and [the recursive Legendre algorithm](https://wandbox.org/permlink/Afl7y4II6Uf3QRqk).  Note that for ease of computation, this Python non-recursive program uses another method based on k-rough numbers and additional tables of the same size as that described containing count intermediate count values so that the Moebius function doesn't need to be evaluated.  Legendre would not have used this implementation as k-rough numbers had not been investigated yet in his time, but an equivalent implementation using a "factors" array of that same size containing values for the least prime factors and the Moebius function values for each represented odd number can be used

A Python-as-pseudocode implementation of the Legendre algorithm using "Partial Sieving" is as follows, written as Legendre might have done it if he used "Partial Sieving":

```
TO BE PROVIDED LATER!!!
```

### Meissel's Prime Counting Algorithm

In 1870, Ernst Meissel published his Prime Counting Algorithm and in 1886 he published the results from using his algorithm to hand calculate the number of primes to a billion (1,000,000,000) for which he made a couple of math errors for a slightly wrong answer.  He must have used Partial Sieving and not recursion to calculate Phi as to do otherwise would have resulted in more division operations than he could have computed in his lifetime, and no prime counting algorithm should really be implemented with the recursive formula if computation time is a concern.  His prime counting formula is very similar to that of Legendre with an important difference:  his `a` value is just the number of primes until the cube root of the counting range instead of the square root and he needs to subtract an easy to compute `P2` term to account for the difference.  However, the "easy-to-compute" term requires that one sieves to `x**(2/3)` instead of `x**(1/2)` (the square root) so the time saved in greatly reduce operations is generally spent sieving so has an asymptotic performance of `O(x**(2/3)(log log x))` although memory use is greatly reduced to `O(x**(1/3))` when implemented with a page-segmented (blocked) sieve.  The Meissel prime counting formula is as follows:

$$
\begin{align*}
\pi(x) &= \phi(x,a) + a - 1 - P2(x) \\
a &= \pi(\sqrt[3]{x}) \\
\phi(x,a) &= {x} - \sum_{i}^{a} {\lfloor \frac{x}{p_i} \rfloor}
                 + \sum_{j>i}^{a} {\lfloor \frac{x}{p_i p_j} \rfloor}
                 - \sum_{k>j>i}^{a} {\lfloor \frac{x}{p_i p_j p_k} \rfloor}
                     \dots + (-1)^a {\lfloor \frac{x}{p_1 p_2 p_3 \dots p_a} \rfloor}  \\
&= \sum_{d \mid P_a} \mu(d) {\lfloor \frac{x}{d} \rfloor} \\
&= \phi(x,a-1) - \phi(x/p_a, a-1) \\
P_a &= \prod_{i=1}^{a} p_i \\
P2(x) &= \sum_{i>a}^{b} \pi({\lfloor \frac{x}{p_i} \rfloor}) - \pi(p_i) + 1 \\
b &= \pi(\sqrt{x}) \\
\end{align*} 
$$

It is very easy to adapt the non-recursive implementations of the Legendre algorithm to use in an implementation of the Meissel algorithm in that the only difference is that the sieving range has to increase from the square root of the counting limit to the counting limit to the two thirds power, the "factors" or "k-roughs" table are reduced from the square root of the counting limit in size to the cube root of the counting limit in size, and the fairly simple calculation of the `P2` correction term using the results from the sieving.

### Lehmer's Prime Counting Algorithm

The first adaptation of prime counting algorithms to computers was made by DH Lehmer in 1959 where he was able to calculate the number of primes to ten billion (10,000,000,00) with an error of just one.  Lehmer was greatly constrained by the limited RAM memory of the mainframe computer to which he had access so his algorithm was a step beyond that of Meissel so as to reduce memory required in that it required the computation of a `P3(x)` term slightly more complex than the `P2(x)` term used in the Meissel algorithm but while it reduced memory requirements to `O(x**(1/4))` it required sieving to `x**(3/4)` and thus would be slower than Meissel's algorithm, made worse in that he chose to use a version of the recursive Phi formula although I believe he had enough RAM memory to be able to implement the non-recursive Phi algorithm due to his lower size of the required "factors" array of `O**(1/4)` using the method I described for the Legendre algorithm.  Other than being the first computer program implementation, I don't believe that the Lehmer algorithm has much bearing on further developments of prime counting algorithms other that his was perhaps the first use of a "Tiny Phi" Look Up Table based shortcut for the values of $$\phi(x,c)$$ where `c` doesn't actually improve performance much for large modern counting ranges but is a useful computational shortcut when one is using wheel factorized sieving where those lower prime values have been eliminated due to the wheel-factorization.

### The Lagarias, Miller, and Odlyzko (LMO) Prime Counting Algorithm

We now have two primary classes of prime counting algorithm:  that of Legendre which doesn't require sieving to a high limit but uses a lot of memory and has a high number of computations, and that of Meissel or Lehmer (Meissel plus) which greatly reduces the amount of memory required and the number of computation operations at the cost of sieving to higher limits (although less than doing a full sieve to the counting limit) and with somewhat higher program complexity, especially when using a page-segmented wheel-factorized sieve to reduce the sieving time as much as possible.  In 1985 LMO introduced some major new innovative ideas:  1) that one could tune the algorithm so as to be between the Legendre and Meissel limits, balancing sieving time against the number of counting operations to some optimum level to make a gain in asymptotic complexity while still keeping memory use fairly low, and 2) they adapted their algorithm to multi-threading use on their two core mainframe computer.  They did use a minimal amount of recursive Phi but not for the main Phi calculation where it would make a difference; they divided the computation into different classes as the trivial "S1" class where recursive computations could be used as they take very little time, and "S2" special leaves along with a slightly modified "P2" calculation.  When using a tuning factor of "one" it turns into a Meissel algorithm where there are no easy special leaves and only hard special leaves, and when a tuning factor of `x**(1/6)` is used, the computation becomes the same as for the Legendre algorithm with the easy special leaves computation becoming much larger than the hard special leaves computation.  The sieve size varies between the limits as for between the Legendre algorithm of to the square root of the counting range to that of the Meissel, but being page-segmented does not require much memory, and the best tuning limit is generally some sizable fraction less than as required for the Legendre algorith.  Thus, memory use is about `O(x**(1/3))` while the tuned asymptotic complexity is a little better than the Meissel algorithm at `O(x**(1/3)/(log x))`.

### The Deleglise and Rivat Improvements to the LMO Algorithm

Deleglise and Rivat published some refinements to the LMO algorithm in 1996 with some further subclassification of the "leaves" into the "S2" easy special leaves and the "S2" hard special leaves with different optimizations for each and were thus able to improve the asymptotic complexity to `O(x**(1/3)/((log x)**2))` while still using about the same amount of memory.  With these refinements the tuning factor can be increased above that used for LMO and the required sieving range is reduced proportionally.

### Gourdon's Prime Counting Algorithm

Finally, Xavier Gourdon published his [ultimate prime counting algorithm in 2001](https://t5k.org/nthprime/) which made some further refinements on top of the Deleglise-Rivat version to not change the asymptotic complexity/performance other than by generally a constant reduction factor due to recognizing some redundant computations that were symetrical and therefore dould be reduced by half and some computations that could be divided into those requiring partial sieving and those that could be computed faster by not requiring partial sieving.  Gourdon's algorithm is broken down into a couple of trivial computations of "Phi0" and "Sigma" (which include the LMO "S1" computation) and the major computations of the "A", "B", "C", and "D" parts with the "B" computation being almost the same as Meissel's "P2" computation, the "A" computaion being the easy special leaves compution of Deleglise-Rivat, and the "C" and "D" computations being the hard special leaves of Deleglise-Rivat with the split being made that "C" computations do not need partial sieving whereas "D" computations do.

Gourdon introduces a further tuning factor on page 7 of his paper linked above in the "alpha-z" tuniny factor he calls `d` where `z = d * y`; with this tuning factor causing a slight decrease in the number of "C" and more especially the hard to compute "D" calculations (of which there aren't all that many) but adds some extra complexity to the algorithm.  In light of that this tuning factor only reduces execution time by less than ten percent at its best and in order to somewhat simplify the implementation, I do not use this tuning factor and set it to one which makes `z = y` in the formulas on page 7 of the Gourdon paper.

There are a couple of errors in the Gordon paper:  1) the "Sigma6" formula of page 6 uses a summation over terms expressed as ${\sqrt{x} / \sqrt{p}}$ which is mathematically identical to $\sqrt{x/p}$ for floating point/real numbers but not the same with the truncation of integer operations so only the latter formula gives the correct answers, and 2) the "Phi0" formula of page 7 has a range for `p` of up to `y` where it is correctly given on page 3 as up to `z`, but for my implementation with a "alpha-z" tuning factor `d` of one, they are one and the same.

## Improvements to the Implementation of the Gourdon Algorithm

In 2006, Tomas Oliveira e Silva published a paper with some implementation improvements to the Deleglise-Rivat algorithm that can also be applied to the Gourdon algorithm:  he implemented some clustering and caching improvements to the computation of the "S2" easy special leaves that improve performance for large counting ranges above about `10**21`.  However, for my much more limited counting range I do not implement this improvements as they can cost performance for these much smaller maximum counting range.

From about 2003 there have been implementations including something like the table calculation that I describe in the non-recursive implementation of calculating "Phi" for the Legendre algorithm.  My implementation here does incorporate those improvements as I implement the table pretty well as described there.

Until much after the paper published in 2006, most prime counting algorithm implementations used a indexed tree/Fenwick Tree implementation to organize the counting of the "hard special leaves" contributions; however, counting with this tree structure has several costs that make it not ideal as in:  1) a high memory use, especially in multi-threading as every thread needs its own version, 2) an additional computational complexity factor of about `log x`, and 3) slow computation due to the poor cache associativity due to the random dual array look-ups required.  Kim Walisch realized these problems since about 2015 in his C++ implementation of "primecount" and introduced the idea of counting the changes in "Pi"/"Phi" due to each partial sieve culling pass by subtracting the number of culling operations that actually cause a bit toggle, with these counts tracked in blocks of sieve buffer words to reduce the time required to do a counting scan across the sieve buffer(s), as the cost of this extra culling complexity, while large, only occurs for a fraction of the total culls required so overall cost isn't all that high and much less than the cost without it, and this cost is directly proportial to the overall computational complexity so it doesn't affect that.  Kim Walich uses a more complex method of keeping track of the counting scans across the sieve buffers but, at least for the smallish prime counting limits used by my implementations, I simplify this to just making the cost of block counting scans track the asymptotic complexity by just increasing the size of the blocks just enough with increasing range to keep the time proportional to the asymptotic complexity.  For the very limited counting ranges in my implementations, I don't even do that but just set a fixed block size large enough so that counting scans have a fairly negligible execution cost across the range used although the "D" counting will be slightly slower than as compared to a smaller block size for smaller counting ranges.  My block size is currently eight 64-bit words or sixteen 32-bit words which is the maximum integer size in JavaScript.

## The SoE Implementation

My SoE implementation uses page-segmentation by fixed 512 KiloByte sieve buffer pages and 2/3/5 wheel-factorization (a wheel circumference of 30 values with 8 residuals that fit into one byte) so the range of one sieve buffer is thirty times 512 Kilobytes or 15,728,640.  Even in JavaScript, this is reasonably fast (about four or five times faster than not using page-segmentation due to much better cache associativity) and gains another factor of about three due to the wheel-factorization.  An implementation in a native compiling language such as C/C++ is at least three times faster due to certain extreme optimiztions that can't be implemented in JavaScript.  This is fast enough so the time spent sieving, even with the special culling/counting for partial sieving, is only a small fraction of the overall computation time, with the remaining time mostly spent counting.

## Problems of JavaScript Counting and Dividing

Current JavaScript/ECMASCRIPT does not have access to CPU hardware "pop counting" instructions so these need to be emulated in software at a fairly large cost in execution time; as well, bit manipulations such as bit shifting, bit and'ing, bit or'ing, etc., can only be performed on 32-bit integers in JavaScript instead of 64-bit integers in native compiling languages.  Further, JavaScript division operations are 64-bit floating point operations instead of what would be native 64-bit divisions for native compiled languages at a cost of about double for even newer CPU's.  Overall, division and counting operations are likely to be about three to five times slower than if translated to a native compiling language such as C/C++.

## Realizations Made While Implementing

The "A" and "C" computations only require sieving up to about the square root of the prime counting limit, while the "B" and "D" computation require sieving from about the square root limit up to the limit of about `x**(2/3)/alpha`.  This implies that the page-segmentation could be broken into these two parts with different optmizations for each.  Having two parts would be especially important if one wwere implementing multi-threading, as the small lower part contains most of the counting point values so is very dense while the upper part has a much lower density of counting points although it does increase from low density to quite high density as pages increase to the upper limit.

## Performance

The tuning factor `alpha` trades fewer segmented sieve blocks to iterate over for more computations per block. From the paper, it must be <= x^(1/6) and should grow like O(log^3 x).  However, empirical tuning and measurements show, at least for the limited counting ranges used for JavaScript, this relationship does not hold and this implementation uses a polynomial curve fit based on the log of the counting range to find something close to the optimum value.  Because sieving is relatively faster than counting for this implementation due to the JavaScript problems described above, it was found that this `alpha` tuning factor could be lowered somewhat below otherwise ideal values to slightly increase the amount of sieving in order to decrease the amount of counting for a slight decrease in computation time.

The following benchmarks use the fixed page-segment size of 512 KiloBytes as described in the sieving section above and `alpha` values optimized as descibed above.  These benchmarks were run on my mini-PC (AMD 7840HS CPU with 32 GigaBytes of RAM, Fedora 44, and Google Chrome Version 150.0.7871.128 (Official Build) (64-bit)); times are mimimums as reported by the web page:

| x                          | Time     | π(x)            |
|----------------------------|----------|-----------------|
| 10**9                      | 2 ms     | 50847534        |
| 10**10                     | 5 ms     | 455052511       |
| 10**11                     | 21 ms    | 4118054813      |
| 10**12                     | 89 ms    | 37607912018     |
| 10**13                     | 304 ms   | 346065536839    |
| 10**14                     | 769 ms   | 3204941750802   |
| 10**15                     | 3045 ms  | 29844570422669  |
| 2**53-1 (9007199254740991) | 11824 ms | 252252704148404 |

It was found that Firefox version 152.0.6 (64-bit) is about 40 percent slower than these Google Chrome results and that for the largest counting ranges that the time is quite variable for both browsers, likely due to Garbage Collection, with the first run after a page reload usually the fastest other than for very small counting ranges where repeated runs can take advantage of hot code reloading better optimizations.

## Conclusion

This is probably about the ultimate prime counting implementation written in JavaScript and is quite practical up to the given `2**53-1` counting range.  One can envisage an application where a client web page could calculate the "Pi" count as well as the nth prime within this range as does [this web page by Andrew Carr and Andrew Booker](https://t5k.org/nthprime/) without requiring the large data base of intermediate count values on a swerver as required by that web page.  There may be other applications...
