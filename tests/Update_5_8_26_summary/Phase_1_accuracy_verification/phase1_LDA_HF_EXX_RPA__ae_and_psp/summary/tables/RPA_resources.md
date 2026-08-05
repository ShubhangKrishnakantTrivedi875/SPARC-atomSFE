## RPA -- peak memory and wall time

### all-electron, loose

| element | new (MiB) | original (MiB) | oneatom (MiB) | new (s) | original (s) | oneatom (s) | mem orig/new | time orig/new |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| He (Z=2) | 397 | 880 | 1672 | 15 | 14 | 16 | 2.2x | 0.9x |
| N (Z=7) | 409 | 1577 | 1740 | 20 | 20 | 19 | 3.9x | 1.0x |
| Mn (Z=25) | 424 | 3247 | 1842 | 176 | 219 | 459 | 7.7x | 1.2x |
| Mo (Z=42) | 428 | 4285 | 1861 | 33 | 42 | 23 | 10.0x | 1.3x |
| Au (Z=79) | 433 | 5689 | 2069 | 40 | 54 | 23 | 13.1x | 1.3x |

### all-electron, tight

| element | new (MiB) | original (MiB) | oneatom (MiB) | new (s) | original (s) | oneatom (s) | mem orig/new | time orig/new |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| He (Z=2) | 1194 | 7080 | 4878 | 291 | 265 | 210 | 5.9x | 0.9x |
| N (Z=7) | 1188 | 19239 | 5173 | 383 | 437 | 233 | 16.2x | 1.1x |
| Mn (Z=25) | 1265 | - | - | 3610 | - | - | - | - |
| Mo (Z=42) | 1217 | - | 5588 | 560 | - | 284 | - | - |
| Au (Z=79) | 1228 | - | 5995 | 794 | - | 369 | - | - |

### pseudopotential, loose

| element | new (MiB) | original (MiB) | oneatom (MiB) | new (s) | original (s) | oneatom (s) | mem orig/new | time orig/new |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| He (Z=2) | 399 | 882 | 1704 | 14 | 14 | 17 | 2.2x | 1.0x |
| N (Z=7) | 418 | 1297 | 1796 | 21 | 24 | 18 | 3.1x | 1.2x |
| Mn (Z=25) | 433 | 2416 | 1853 | 149 | 206 | 424 | 5.6x | 1.4x |
| Mo (Z=42) | 432 | - | 1840 | 34 | - | 31 | - | - |
| Au (Z=79) | 431 | 2431 | 1843 | 35 | 221 | 29 | 5.6x | 6.3x |

### pseudopotential, tight

| element | new (MiB) | original (MiB) | oneatom (MiB) | new (s) | original (s) | oneatom (s) | mem orig/new | time orig/new |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| He (Z=2) | 1266 | 7035 | 4948 | 283 | 263 | 205 | 5.6x | 0.9x |
| N (Z=7) | 1280 | 12953 | 5257 | 382 | 538 | 248 | 10.1x | 1.4x |
| Mn (Z=25) | 1258 | - | - | 2989 | - | - | - | - |
| Mo (Z=42) | 1276 | - | 5556 | 762 | - | 481 | - | - |
| Au (Z=79) | 1218 | - | 5524 | 743 | - | 466 | - | - |

