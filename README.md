# cfErdosMoser
The Erdős-Moser equation revisited using continued fractions

See Yves Gallot, Pieter Moree, Wadim Zudilin,  
The Erdős-Moser equation 1<sup>*k*</sup> + 2<sup>*k*</sup> + ... + (*m*−1)<sup>*k*</sup> = *m*<sup>*k*</sup> revisited using continued fractions.  
[Math. Comp. **80** (2011), 1221–1237](https://www.ams.org/journals/mcom/2011-80-274/S0025-5718-2010-02439-1/).

## About
This new version of program is not based on the binary expansion of log 2. Instead, the main diagonal Padé approximant to log 2 at the point one is calculated. It can be computed using a generalized continued fraction but then convergents are not irreducible fractions. However, the canonical form can be obtained by dividing the numerator and denominator by a known function.  
If the normalized convergents of the generalized continued fraction and the coefficients of the regular continued fraction are evaluated together then the size of remaining *j*<sup> th</sup> fractions is smaller than the size of the binary expansion needed for *j* coefficients. Because memory size is the main limitation of this computation, this method is more efficient.  
If *q*<sub>*j*</sub> is a *n*-digit number, memory usage is about 3*n*.

## Results

The following table provides the smallest integers *j* satisfying conditions (a), (b) and (c) of Theorem 2 and (d) is checked for *p* &le; 43 such that 3 is a primitive root modulo *p*.  
2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>2</sup> &middot; 7,
2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>2</sup> &middot; 7<sup>2</sup> and
2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>2</sup> &middot; 7<sup>3</sup> can be tested to improve the bound of Theorem 3.


| *N* | *j*<sub>*N*</sub> | *a*<sub>*j*+1</sub> | *q*<sub>*j*</sub> <span style="font-weight: normal">(rounded down)</span> | *q*<sub>*j*</sub> <span style="font-weight: normal">mod 6</span> | *p*<span style="font-weight: normal">(*q*<sub>*j*</sub>)</span> |
|:---:| ---:| ---:|:--- |:---:|:--- |
| 1             |     642 |     764 | 2.383153 &middot; 10<sup>   330</sup> | -1 | |
| 2             |     664 |   1 529 | 2.383153 &middot; 10<sup>   330</sup> | -1 | |
| 2<sup>2</sup> |   1 254 |  21 966 | 1.132014 &middot; 10<sup>   638</sup> | +1 | 5 |
|               |   3 468 |   1 373 | 8.283906 &middot; 10<sup> 1 791</sup> | +1 | |
| 2<sup>3</sup> |   1 264 |  43 933 | 1.132014 &middot; 10<sup>   638</sup> | +1 | 5 |
|               |   4 876 |  25 832 | 6.968069 &middot; 10<sup> 2 489</sup> | -1 | |
| 2<sup>4</sup> |   1 280 |  87 866 | 1.132014 &middot; 10<sup>   638</sup> | +1 | 5 |
|               |   5 952 |   4 490 | 1.222504 &middot; 10<sup> 3 079</sup> | +1 | |
| 2<sup>5</sup> |   1 294 | 175 733 | 1.132014 &middot; 10<sup>   638</sup> | +1 | 5 |
|               |   6 062 |   8 981 | 1.222504 &middot; 10<sup> 3 079</sup> | +1 | |
| 2<sup>6</sup> |   8 950 |  26 416 | 3.458446 &middot; 10<sup> 4 589</sup> | -1 | |
| 2<sup>7</sup> |   8 926 |  52 834 | 3.458446 &middot; 10<sup> 4 589</sup> | -1 | |
| 2<sup>8</sup> | 119 476 | 122 799 | 1.374540 &middot; 10<sup>61 317</sup> | +1 | |
| 2<sup>8</sup> &middot; 3             |     119 008 |    368 398 | 1.374540 &middot; 10<sup>    61 317</sup> | +1 | |
| 2<sup>8</sup> &middot; 3<sup>2</sup> |     139 532 |    782 152 | 9.351282 &middot; 10<sup>    71 882</sup> | +1 | |
| 2<sup>8</sup> &middot; 3<sup>3</sup> |   6 168 634 |  1 540 283 | 8.220719 &middot; 10<sup> 3 177 670</sup> | +1 | |
| 2<sup>8</sup> &middot; 3<sup>4</sup> |  22 383 618 |  5 167 079 | 5.128265 &middot; 10<sup>11 538 265</sup> | +1 | 17 |
| 2<sup>8</sup> &middot; 3<sup>4</sup> |  54 192 270 |  4 646 574 | 1.060113 &middot; 10<sup>27 929 457</sup> | -1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> | 155 830 946 | 31 664 035 | 2.257099 &middot; 10<sup>80 303 211</sup> | -1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5             |   351 661 538 |    85 898 211 | 9.729739 &middot; 10<sup>  181 214 202</sup> | -1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>2</sup> | 1 738 154 976 | 1 433 700 727 | 1.594940 &middot; 10<sup>  895 721 905</sup> | +1 | 5 |
|                                                             | 1 977 626 256 |   853 324 651 | 1.196828 &middot; 10<sup>1 019 133 881</sup> | -1 | |
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>3</sup> | 2 015 279 170 |  4 388 327 617 | 5.565196 &middot; 10<sup>1 038 523 018</sup> | -1 | 19 |
|                                                             | 3 236 170 820 |  2 307 115 390 | 5.427816 &middot; 10<sup>1 667 658 416</sup> | +1 | | 
| 2<sup>8</sup> &middot; 3<sup>5</sup> &middot; 5<sup>4</sup> | 2 015 385 392 | 21 941 638 090 | 5.565196 &middot; 10<sup>1 038 523 018</sup> | -1 | 19 |
|                                                             | 3 236 257 942 | 11 535 576 954 | 5.427816 &middot; 10<sup>1 667 658 416</sup> | +1 | |
