# GhostBuster 0.1.1 - benchmark on real stacks

Automated mode, default selection (white stars).
Figures in the ghost zone (22-40 px at the reference scale). Real stacks have no ground truth: these
figures describe what the correction did and are meant for comparison between versions; the contact
sheets are the reference.

- **removed**: light removed, % of the light above the local background
- **uneven**: spread of the removed light over 8 sectors, (max - min) / mean in %; 0 = round ghost, large = lopsided
- **colour**: largest change of a channel share between raw and corrected (points)

## DRACO RAW Estern Veil Nebula 36x300s 

4096 x 3072, 3 channel(s); scale 1.00 (astrometric solution (1.45"/px)); 9 candidate stars, 1 corrected in 2.2 s.

![DRACO RAW Estern Veil Nebula 36x300s ](DRACO_RAW_Estern_Veil_Nebula_36x300s_.png)

Contact sheet: the 1 largest corrected stars, raw | corrected, two per row, in the order of the table.

| # | star (x, y) | confidence | removed % | uneven % | colour pts |
|---|---|---|---|---|---|
| 1 | (1298, 50) | 2.1 | 9.7 | 4.3 | 0.3 |

## DRACO RAW IC 1805 36x300s

4096 x 3072, 3 channel(s); scale 1.00 (image size (4096 x 3072)); 16 candidate stars, 11 corrected in 5.5 s.

![DRACO RAW IC 1805 36x300s](DRACO_RAW_IC_1805_36x300s.png)

Contact sheet: the 8 largest corrected stars, raw | corrected, two per row, in the order of the table.

| # | star (x, y) | confidence | removed % | uneven % | colour pts |
|---|---|---|---|---|---|
| 1 | (4027, 2269) | 4.6 | 20.5 | 65.9 | 0.4 |
| 2 | (2375, 954) | 4.8 | 31.9 | 35.7 | 2.5 |
| 3 | (32, 424) | 2.4 | 19.6 | 100.1 | 1.9 |
| 4 | (131, 1590) | 5.0 | 15.5 | 95.0 | 1.1 |
| 5 | (2526, 823) | 2.0 | 19.7 | 91.8 | 7.0 |
| 6 | (1959, 1332) | 3.6 | 9.3 | 5.8 | 2.4 |
| 7 | (2236, 2425) | 2.0 | 20.8 | 96.2 | 1.3 |
| 8 | (891, 1941) | 3.7 | 32.8 | 98.4 | 0.4 |
| 9 | (2271, 3042) | 1.6 | 1.6 | 60.1 | 0.5 |
| 10 | (172, 414) | 1.8 | 5.9 | 58.5 | 1.2 |
| 11 | (2426, 1318) | 1.7 | 26.6 | 30.9 | 18.9 |

## DRACO RAW M 42 36x300s_Duo-Band

4096 x 3072, 3 channel(s); scale 1.00 (astrometric solution (1.45"/px)); 20 candidate stars, 11 corrected in 7.5 s.

![DRACO RAW M 42 36x300s_Duo-Band](DRACO_RAW_M_42_36x300s_Duo-Band.png)

Contact sheet: the 8 largest corrected stars, raw | corrected, two per row, in the order of the table.

| # | star (x, y) | confidence | removed % | uneven % | colour pts |
|---|---|---|---|---|---|
| 1 | (1035, 2364) | 13.2 | 22.5 | 60.9 | 1.6 |
| 2 | (1031, 2695) | 8.2 | 1.9 | 153.4 | 0.2 |
| 3 | (3003, 587) | 3.5 | 15.2 | 148.2 | 1.3 |
| 4 | (2856, 493) | 5.9 | 22.8 | 157.3 | 0.7 |
| 5 | (269, 2101) | 6.6 | 20.5 | 71.1 | 0.6 |
| 6 | (138, 1499) | 5.3 | 19.6 | 70.5 | 1.5 |
| 7 | (3639, 16) | 3.7 | 20.3 | 190.0 | 1.3 |
| 8 | (884, 2648) | 3.5 | 29.5 | 57.9 | 0.8 |
| 9 | (3477, 974) | 2.9 | 20.2 | 92.1 | 1.4 |
| 10 | (1109, 2505) | 3.2 | 58.0 | 181.9 | 22.0 |
| 11 | (3997, 378) | 2.0 | 19.9 | 97.4 | 1.9 |

## DRACO RAW Veil Nebula_31x300s120_Duo-Band

4096 x 3072, 3 channel(s); scale 1.00 (astrometric solution (1.45"/px)); 4 candidate stars, 2 corrected in 1.3 s.

![DRACO RAW Veil Nebula_31x300s120_Duo-Band](DRACO_RAW_Veil_Nebula_31x300s120_Duo-Band.png)

Contact sheet: the 2 largest corrected stars, raw | corrected, two per row, in the order of the table.

| # | star (x, y) | confidence | removed % | uneven % | colour pts |
|---|---|---|---|---|---|
| 1 | (1672, 2324) | 7.0 | 9.2 | 109.1 | 2.5 |
| 2 | (1596, 739) | 2.6 | 17.8 | 22.2 | 0.4 |

