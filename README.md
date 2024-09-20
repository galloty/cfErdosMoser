# cfErdosMoser
The Erdős-Moser equation revisited using continued fractions

See Yves Gallot, Pieter Moree, Wadim Zudilin,  
The Erdős-Moser equation 1<sup>*k*</sup>&nbsp;+&nbsp;2<sup>*k*</sup>&nbsp;+&nbsp;...&nbsp;+&nbsp;(*m*−1)<sup>*k*</sup> = *m*<sup>*k*</sup> revisited using continued fractions.  
[Math. Comp. **80** (2011), 1221–1237](https://www.ams.org/journals/mcom/2011-80-274/S0025-5718-2010-02439-1/).

On September 2024, Robert Gahan established the new benchmark *m* > 10<sup>10<sup>10</sup></sup> with *cfErdosMoser*. The processor was an Intel&reg; Core&trade; i9-13980HX with 64GB RAM, running Windows 11 Pro. It took about 61 hours with 8 threads and the peak memory usage was 33GB to prove, using *N* = 3&nbsp;&middot;&nbsp;19#, that conditions (a), (b), (c) and partly (d) are not satisfied if *q*<sub>*j*</sub> < 7.76&nbsp;&middot;&nbsp;10<sup>&nbsp;10&nbsp;381&nbsp;635&nbsp;903</sup>. With *N* = 2<sup>8</sup>&nbsp;&middot;&nbsp;3<sup>5</sup>&nbsp;&middot;&nbsp;5<sup>4</sup>&nbsp;&middot;&nbsp;7, it took about 91 hours and the peak memory usage was 56GB to find the solution *q*<sub>*j*</sub> = 7.69&nbsp;&middot;&nbsp;10<sup>&nbsp;16&nbsp;771&nbsp;044&nbsp;760</sup> that fits conditions (a), (b), (c) and (d).

**Theorem&nbsp;3** (2024). If an integer pair (*m*,&nbsp;*k*) with *k*&nbsp;&ge;&nbsp;2 satisfies the Erdős-Moser equation, then *m*&nbsp;>&nbsp;*k*&nbsp;>&nbsp;10<sup>10<sup>10</sup></sup>.

## About
This new version of the program is not based on the binary expansion of log&nbsp;2. Instead, the main diagonal Padé approximant to log&nbsp;2 at the point one is calculated. It can be computed using a generalized continued fraction but then convergents are not irreducible fractions. However, the canonical form can be obtained by dividing the numerator and denominator by a known function.  
If the normalized convergents of the generalized continued fraction and the coefficients of the regular continued fraction are evaluated together then the size of the remaining *j*<sup>&nbsp;th</sup> and (*j*&nbsp;+&nbsp;1)<sup>th</sup> fractions is smaller than the size of the binary expansion needed for evaluating *j* coefficients. Because memory size is the main limitation of this computation, this method is more efficient.  
If *q*<sub>*j*</sub> is a *n*-digit number, memory usage is about 4*n*.

## Results

The following table provides the smallest integers *j* satisfying conditions (a), (b) and (c) of Theorem&nbsp;2 and (d) is checked for *p*&nbsp;&le;&nbsp;43 such that 3 is a primitive root modulo&nbsp;*p*.  
*n*# is the [primorial](https://en.wikipedia.org/wiki/Primorial) function.

The solution to 2<sup>8</sup>&nbsp;&middot;&nbsp;3<sup>5</sup>&nbsp;&middot;&nbsp;5<sup>4</sup>&nbsp;&middot;&nbsp;7 = 272160000, found by Robert Gahan, is the best known result.  
23# = 223092870 is the most wanted candidate to improve the bound of Theorem&nbsp;3.

| *N* | *j*<sub>*N*</sub> | *a*<sub>*j*+1</sub> | *q*<sub>*j*</sub> <span style="font-weight: normal">(rounded down)</span> | *q*<sub>*j*</sub> <span style="font-weight: normal">mod 6</span> | *p*<span style="font-weight: normal">(*q*<sub>*j*</sub>)</span> |
|:---:| ---:| ---:|:--- |:---:|:---:|
| 23# | ? | ? | | | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>4</sup> &middot; 7 | 13 988 361 422 | 107 695 241 629 | ~1.444542 &middot; 10<sup> 7 208 444 672</sup>~ | -1 | 7 |
|                                                                        | 32 544 869 308 |  65 416 187 470 | 7.699688 &middot; 10<sup> 16 771 044 760</sup> | -1 | |
| 3 &middot;19# | 20 145 968 414 | 8 176 658 150 | 7.766616 &middot; 10<sup> 10 381 635 903</sup> | +1 | |
| 2 &middot;19# | 13 078 292 184 | 5 890 259 269 | 6.066802 &middot; 10<sup>  6 739 544 561</sup> | -1 | |
| 19# | 13 050 765 960 |  4 245 878 986 | 1.202797 &middot; 10<sup> 6 725 335 083</sup> | -1 | |
| 11! |  5 151 307 662 | 14 371 269 642 | 1.086746 &middot; 10<sup> 2 654 644 376</sup> | -1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>4</sup> | 2 015 385 392 | 21 941 638 090 | ~5.565196 &middot; 10<sup> 1 038 523 018</sup>~ | -1 | 19 |
|                                                             | 3 236 257 942 | 11 535 576 954 |  5.427816 &middot; 10<sup> 1 667 658 416</sup>  | +1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>3</sup> | 2 015 279 170 |  4 388 327 617 | ~5.565196 &middot; 10<sup> 1 038 523 018</sup>~ | -1 | 19 |
|                                                             | 3 236 170 820 |  2 307 115 390 |  5.427816 &middot; 10<sup> 1 667 658 416</sup>  | +1 | | 
| 10! | 1 738 197 230 | 3 345 301 699 | ~1.594940 &middot; 10<sup>   895 721 905</sup>~ | +1 | 5 |
|     | 2 579 384 226 | 1 116 890 723 |  5.498581 &middot; 10<sup> 1 329 206 901</sup>  | +1 |   |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>2</sup> | 1 738 154 976 | 1 433 700 727 | ~1.594940 &middot; 10<sup>   895 721 905</sup>~ | +1 | 5 |
|                                                             | 1 977 626 256 |   853 324 651 |  1.196828 &middot; 10<sup> 1 019 133 881</sup>  | -1 | |
| 17# | 367 396 948 | 127 090 830 | ~4.724143 &middot; 10<sup> 189 330 916</sup>~ | -1 |  5 |
|     | 586 036 322 | 331 038 577 | ~9.914162 &middot; 10<sup> 302 014 571</sup>~ | +1 | 19 |
|     | 639 808 912 | 675 894 301 |  2.276299 &middot; 10<sup> 329 726 891</sup>  | +1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5 | 351 661 538 | 85 898 211 | 9.729739 &middot; 10<sup> 181 214 202</sup> | -1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> | 155 830 946 | 31 664 035 | 2.257099 &middot; 10<sup> 80 303 211</sup> | -1 | |
