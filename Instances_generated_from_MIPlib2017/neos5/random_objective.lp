\ENCODING=ISO-8859-1
\Problem name: neos5.mps

Maximize
 r_0: 0 d1               -1 d2               + 1 d3               + 0 d4               + 0 d5               + 1 d6               + 0 d7               + 0 d8               + 1 d9               -1 d10              -1 d11              + 0 d12              + 0 d13              -1 d14              + 1 d15              -1 d16              -1 d17              + 1 d18              -1 d19              + 0 d20              -1 d21              + 1 d22              + 0 d23              + 1 d24              -1 d25              -1 d26              -1 d27              -1 d28              + 1 d29              + 0 d30              + 0 d31              -1 d32              + 1 d33              -1 d34              + 1 d35              + 0 d36              + 1 d37              + 0 d38              -1 d39              + 1 d40              + 0 d41              -1 d42              + 0 d43              -1 d44              -1 d45              -1 d46              + 0 d47              + 0 d48              -1 d49              + 1 d50              + 1 d51              -1 d52              + 1 d53              + 1 d54              -1 d55              -1 d56              + 1 d57              + 0 d58              -1 d59              + 1 d60              + 0 d61              + 0 d62              + 1 d63
Subject To
 R1:  d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 + d13 + d14
      + d15 + d16 + d17 + d18 + d19 + d20 + d21 + d22 + d23 + d24 + d25 + d26
      + d27 + d28 + d29 + d30 + d31 + d32 >= 8
 R2:  d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 + d13 + d14
      + d15 + d16 + d33 + d34 + d35 + d36 + d37 + d38 + d39 + d40 + d41 + d42
      + d43 + d44 + d45 + d46 + d47 + d48 >= 8
 R3:  d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d17 + d18 + d19 + d20 + d21 + d22
      + d23 + d24 + d33 + d34 + d35 + d36 + d37 + d38 + d39 + d40 + d49 + d50
      + d51 + d52 + d53 + d54 + d55 + d56 >= 8
 R4:  d1 + d2 + d3 + d4 + d9 + d10 + d11 + d12 + d17 + d18 + d19 + d20 + d25
      + d26 + d27 + d28 + d33 + d34 + d35 + d36 + d41 + d42 + d43 + d44 + d49
      + d50 + d51 + d52 + d57 + d58 + d59 + d60 >= 8
 R5:  d1 + d2 + d5 + d6 + d9 + d10 + d13 + d14 + d17 + d18 + d21 + d22 + d25
      + d26 + d29 + d30 + d33 + d34 + d37 + d38 + d41 + d42 + d45 + d46 + d49
      + d50 + d53 + d54 + d57 + d58 + d61 + d62 >= 8
 R6:  d1 + d3 + d5 + d7 + d9 + d11 + d13 + d15 + d17 + d19 + d21 + d23 + d25
      + d27 + d29 + d31 + d33 + d35 + d37 + d39 + d41 + d43 + d45 + d47 + d49
      + d51 + d53 + d55 + d57 + d59 + d61 + d63 >= 8
 R7:  d17 + d18 + d19 + d20 + d21 + d22 + d23 + d24 + d25 + d26 + d27 + d28
      + d29 + d30 + d31 + d32 + d33 + d34 + d35 + d36 + d37 + d38 + d39 + d40
      + d41 + d42 + d43 + d44 + d45 + d46 + d47 + d48 >= 7
 R8:  d9 + d10 + d11 + d12 + d13 + d14 + d15 + d16 + d25 + d26 + d27 + d28
      + d29 + d30 + d31 + d32 + d33 + d34 + d35 + d36 + d37 + d38 + d39 + d40
      + d49 + d50 + d51 + d52 + d53 + d54 + d55 + d56 >= 7
 R9:  d9 + d10 + d11 + d12 + d13 + d14 + d15 + d16 + d17 + d18 + d19 + d20
      + d21 + d22 + d23 + d24 + d41 + d42 + d43 + d44 + d45 + d46 + d47 + d48
      + d49 + d50 + d51 + d52 + d53 + d54 + d55 + d56 >= 7
 R10: d5 + d6 + d7 + d8 + d13 + d14 + d15 + d16 + d21 + d22 + d23 + d24 + d29
      + d30 + d31 + d32 + d33 + d34 + d35 + d36 + d41 + d42 + d43 + d44 + d49
      + d50 + d51 + d52 + d57 + d58 + d59 + d60 >= 7
 R11: d5 + d6 + d7 + d8 + d13 + d14 + d15 + d16 + d17 + d18 + d19 + d20 + d25
      + d26 + d27 + d28 + d37 + d38 + d39 + d40 + d45 + d46 + d47 + d48 + d49
      + d50 + d51 + d52 + d57 + d58 + d59 + d60 >= 7
 R12: d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 + d21 + d22 + d23 + d24 + d25
      + d26 + d27 + d28 + d37 + d38 + d39 + d40 + d41 + d42 + d43 + d44 + d53
      + d54 + d55 + d56 + d57 + d58 + d59 + d60 >= 7
 R13: d3 + d4 + d7 + d8 + d11 + d12 + d15 + d16 + d19 + d20 + d23 + d24 + d27
      + d28 + d31 + d32 + d33 + d34 + d37 + d38 + d41 + d42 + d45 + d46 + d49
      + d50 + d53 + d54 + d57 + d58 + d61 + d62 >= 7
 R14: d3 + d4 + d7 + d8 + d11 + d12 + d15 + d16 + d17 + d18 + d21 + d22 + d25
      + d26 + d29 + d30 + d35 + d36 + d39 + d40 + d43 + d44 + d47 + d48 + d49
      + d50 + d53 + d54 + d57 + d58 + d61 + d62 >= 7
 R15: d3 + d4 + d7 + d8 + d9 + d10 + d13 + d14 + d19 + d20 + d23 + d24 + d25
      + d26 + d29 + d30 + d35 + d36 + d39 + d40 + d41 + d42 + d45 + d46 + d51
      + d52 + d55 + d56 + d57 + d58 + d61 + d62 >= 7
 R16: d3 + d4 + d5 + d6 + d11 + d12 + d13 + d14 + d19 + d20 + d21 + d22 + d27
      + d28 + d29 + d30 + d35 + d36 + d37 + d38 + d43 + d44 + d45 + d46 + d51
      + d52 + d53 + d54 + d59 + d60 + d61 + d62 >= 7
 R17: d2 + d4 + d6 + d8 + d10 + d12 + d14 + d16 + d18 + d20 + d22 + d24 + d26
      + d28 + d30 + d32 + d33 + d35 + d37 + d39 + d41 + d43 + d45 + d47 + d49
      + d51 + d53 + d55 + d57 + d59 + d61 + d63 >= 7
 R18: d2 + d4 + d6 + d8 + d10 + d12 + d14 + d16 + d17 + d19 + d21 + d23 + d25
      + d27 + d29 + d31 + d34 + d36 + d38 + d40 + d42 + d44 + d46 + d48 + d49
      + d51 + d53 + d55 + d57 + d59 + d61 + d63 >= 7
 R19: d2 + d4 + d6 + d8 + d9 + d11 + d13 + d15 + d18 + d20 + d22 + d24 + d25
      + d27 + d29 + d31 + d34 + d36 + d38 + d40 + d41 + d43 + d45 + d47 + d50
      + d52 + d54 + d56 + d57 + d59 + d61 + d63 >= 7
 R20: d2 + d4 + d5 + d7 + d10 + d12 + d13 + d15 + d18 + d20 + d21 + d23 + d26
      + d28 + d29 + d31 + d34 + d36 + d37 + d39 + d42 + d44 + d45 + d47 + d50
      + d52 + d53 + d55 + d58 + d60 + d61 + d63 >= 7
 R21: d2 + d3 + d6 + d7 + d10 + d11 + d14 + d15 + d18 + d19 + d22 + d23 + d26
      + d27 + d30 + d31 + d34 + d35 + d38 + d39 + d42 + d43 + d46 + d47 + d50
      + d51 + d54 + d55 + d58 + d59 + d62 + d63 >= 7
 R22: d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d25 + d26 + d27 + d28 + d29 + d30
      + d31 + d32 + d41 + d42 + d43 + d44 + d45 + d46 + d47 + d48 + d49 + d50
      + d51 + d52 + d53 + d54 + d55 + d56 >= 6
 R23: d1 + d2 + d3 + d4 + d9 + d10 + d11 + d12 + d21 + d22 + d23 + d24 + d29
      + d30 + d31 + d32 + d37 + d38 + d39 + d40 + d45 + d46 + d47 + d48 + d49
      + d50 + d51 + d52 + d57 + d58 + d59 + d60 >= 6
 R24: d1 + d2 + d3 + d4 + d13 + d14 + d15 + d16 + d17 + d18 + d19 + d20 + d29
      + d30 + d31 + d32 + d37 + d38 + d39 + d40 + d41 + d42 + d43 + d44 + d53
      + d54 + d55 + d56 + d57 + d58 + d59 + d60 >= 6
 R25: d1 + d2 + d3 + d4 + d13 + d14 + d15 + d16 + d21 + d22 + d23 + d24 + d25
      + d26 + d27 + d28 + d33 + d34 + d35 + d36 + d45 + d46 + d47 + d48 + d53
      + d54 + d55 + d56 + d57 + d58 + d59 + d60 >= 6
 R26: d1 + d2 + d5 + d6 + d9 + d10 + d13 + d14 + d19 + d20 + d23 + d24 + d27
      + d28 + d31 + d32 + d35 + d36 + d39 + d40 + d43 + d44 + d47 + d48 + d49
      + d50 + d53 + d54 + d57 + d58 + d61 + d62 >= 6
 R27: d1 + d2 + d5 + d6 + d11 + d12 + d15 + d16 + d17 + d18 + d21 + d22 + d27
      + d28 + d31 + d32 + d35 + d36 + d39 + d40 + d41 + d42 + d45 + d46 + d51
      + d52 + d55 + d56 + d57 + d58 + d61 + d62 >= 6
 R28: d1 + d2 + d5 + d6 + d11 + d12 + d15 + d16 + d19 + d20 + d23 + d24 + d25
      + d26 + d29 + d30 + d33 + d34 + d37 + d38 + d43 + d44 + d47 + d48 + d51
      + d52 + d55 + d56 + d57 + d58 + d61 + d62 >= 6
 R29: d1 + d2 + d7 + d8 + d9 + d10 + d15 + d16 + d17 + d18 + d23 + d24 + d25
      + d26 + d31 + d32 + d35 + d36 + d37 + d38 + d43 + d44 + d45 + d46 + d51
      + d52 + d53 + d54 + d59 + d60 + d61 + d62 >= 6
 R30: d1 + d2 + d7 + d8 + d9 + d10 + d15 + d16 + d19 + d20 + d21 + d22 + d27
      + d28 + d29 + d30 + d33 + d34 + d39 + d40 + d41 + d42 + d47 + d48 + d51
      + d52 + d53 + d54 + d59 + d60 + d61 + d62 >= 6
 R31: d1 + d2 + d7 + d8 + d11 + d12 + d13 + d14 + d17 + d18 + d23 + d24 + d27
      + d28 + d29 + d30 + d33 + d34 + d39 + d40 + d43 + d44 + d45 + d46 + d49
      + d50 + d55 + d56 + d59 + d60 + d61 + d62 >= 6
 R32: d1 + d3 + d5 + d7 + d9 + d11 + d13 + d15 + d18 + d20 + d22 + d24 + d26
      + d28 + d30 + d32 + d34 + d36 + d38 + d40 + d42 + d44 + d46 + d48 + d49
      + d51 + d53 + d55 + d57 + d59 + d61 + d63 >= 6
 R33: d1 + d3 + d5 + d7 + d10 + d12 + d14 + d16 + d17 + d19 + d21 + d23 + d26
      + d28 + d30 + d32 + d34 + d36 + d38 + d40 + d41 + d43 + d45 + d47 + d50
      + d52 + d54 + d56 + d57 + d59 + d61 + d63 >= 6
 R34: d1 + d3 + d5 + d7 + d10 + d12 + d14 + d16 + d18 + d20 + d22 + d24 + d25
      + d27 + d29 + d31 + d33 + d35 + d37 + d39 + d42 + d44 + d46 + d48 + d50
      + d52 + d54 + d56 + d57 + d59 + d61 + d63 >= 6
 R35: d1 + d3 + d6 + d8 + d9 + d11 + d14 + d16 + d17 + d19 + d22 + d24 + d25
      + d27 + d30 + d32 + d34 + d36 + d37 + d39 + d42 + d44 + d45 + d47 + d50
      + d52 + d53 + d55 + d58 + d60 + d61 + d63 >= 6
 R36: d1 + d3 + d6 + d8 + d9 + d11 + d14 + d16 + d18 + d20 + d21 + d23 + d26
      + d28 + d29 + d31 + d33 + d35 + d38 + d40 + d41 + d43 + d46 + d48 + d50
      + d52 + d53 + d55 + d58 + d60 + d61 + d63 >= 6
 R37: d1 + d3 + d6 + d8 + d10 + d12 + d13 + d15 + d17 + d19 + d22 + d24 + d26
      + d28 + d29 + d31 + d33 + d35 + d38 + d40 + d42 + d44 + d45 + d47 + d49
      + d51 + d54 + d56 + d58 + d60 + d61 + d63 >= 6
 R38: d1 + d4 + d5 + d8 + d9 + d12 + d13 + d16 + d17 + d20 + d21 + d24 + d25
      + d28 + d29 + d32 + d34 + d35 + d38 + d39 + d42 + d43 + d46 + d47 + d50
      + d51 + d54 + d55 + d58 + d59 + d62 + d63 >= 6
 R39: d1 + d4 + d5 + d8 + d9 + d12 + d13 + d16 + d18 + d19 + d22 + d23 + d26
      + d27 + d30 + d31 + d33 + d36 + d37 + d40 + d41 + d44 + d45 + d48 + d50
      + d51 + d54 + d55 + d58 + d59 + d62 + d63 >= 6
 R40: d1 + d4 + d5 + d8 + d10 + d11 + d14 + d15 + d17 + d20 + d21 + d24 + d26
      + d27 + d30 + d31 + d33 + d36 + d37 + d40 + d42 + d43 + d46 + d47 + d49
      + d52 + d53 + d56 + d58 + d59 + d62 + d63 >= 6
 R41: d1 + d4 + d6 + d7 + d9 + d12 + d14 + d15 + d17 + d20 + d22 + d23 + d25
      + d28 + d30 + d31 + d33 + d36 + d38 + d39 + d41 + d44 + d46 + d47 + d49
      + d52 + d54 + d55 + d57 + d60 + d62 + d63 >= 6
 R42: d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 + d17 + d18 + d19 + d20 + d29
      + d30 + d31 + d32 + d33 + d34 + d35 + d36 + d45 + d46 + d47 + d48 + d53
      + d54 + d55 + d56 + d57 + d58 + d59 + d60 >= 5
 R43: d3 + d4 + d7 + d8 + d9 + d10 + d13 + d14 + d17 + d18 + d21 + d22 + d27
      + d28 + d31 + d32 + d33 + d34 + d37 + d38 + d43 + d44 + d47 + d48 + d51
      + d52 + d55 + d56 + d57 + d58 + d61 + d62 >= 5
 R44: d3 + d4 + d5 + d6 + d11 + d12 + d13 + d14 + d17 + d18 + d23 + d24 + d25
      + d26 + d31 + d32 + d33 + d34 + d39 + d40 + d41 + d42 + d47 + d48 + d51
      + d52 + d53 + d54 + d59 + d60 + d61 + d62 >= 5
 R45: d3 + d4 + d5 + d6 + d9 + d10 + d15 + d16 + d19 + d20 + d21 + d22 + d25
      + d26 + d31 + d32 + d33 + d34 + d39 + d40 + d43 + d44 + d45 + d46 + d49
      + d50 + d55 + d56 + d59 + d60 + d61 + d62 >= 5
 R46: d3 + d4 + d5 + d6 + d9 + d10 + d15 + d16 + d17 + d18 + d23 + d24 + d27
      + d28 + d29 + d30 + d35 + d36 + d37 + d38 + d41 + d42 + d47 + d48 + d49
      + d50 + d55 + d56 + d59 + d60 + d61 + d62 >= 5
 R47: d2 + d4 + d6 + d8 + d9 + d11 + d13 + d15 + d17 + d19 + d21 + d23 + d26
      + d28 + d30 + d32 + d33 + d35 + d37 + d39 + d42 + d44 + d46 + d48 + d50
      + d52 + d54 + d56 + d57 + d59 + d61 + d63 >= 5
 R48: d2 + d4 + d5 + d7 + d10 + d12 + d13 + d15 + d17 + d19 + d22 + d24 + d25
      + d27 + d30 + d32 + d33 + d35 + d38 + d40 + d41 + d43 + d46 + d48 + d50
      + d52 + d53 + d55 + d58 + d60 + d61 + d63 >= 5
 R49: d2 + d4 + d5 + d7 + d9 + d11 + d14 + d16 + d18 + d20 + d21 + d23 + d25
      + d27 + d30 + d32 + d33 + d35 + d38 + d40 + d42 + d44 + d45 + d47 + d49
      + d51 + d54 + d56 + d58 + d60 + d61 + d63 >= 5
 R50: d2 + d4 + d5 + d7 + d9 + d11 + d14 + d16 + d17 + d19 + d22 + d24 + d26
      + d28 + d29 + d31 + d34 + d36 + d37 + d39 + d41 + d43 + d46 + d48 + d49
      + d51 + d54 + d56 + d58 + d60 + d61 + d63 >= 5
 R51: d2 + d3 + d6 + d7 + d10 + d11 + d14 + d15 + d17 + d20 + d21 + d24 + d25
      + d28 + d29 + d32 + d33 + d36 + d37 + d40 + d41 + d44 + d45 + d48 + d50
      + d51 + d54 + d55 + d58 + d59 + d62 + d63 >= 5
 R52: d2 + d3 + d6 + d7 + d9 + d12 + d13 + d16 + d18 + d19 + d22 + d23 + d25
      + d28 + d29 + d32 + d33 + d36 + d37 + d40 + d42 + d43 + d46 + d47 + d49
      + d52 + d53 + d56 + d58 + d59 + d62 + d63 >= 5
 R53: d2 + d3 + d6 + d7 + d9 + d12 + d13 + d16 + d17 + d20 + d21 + d24 + d26
      + d27 + d30 + d31 + d34 + d35 + d38 + d39 + d41 + d44 + d45 + d48 + d49
      + d52 + d53 + d56 + d58 + d59 + d62 + d63 >= 5
 R54: d2 + d3 + d5 + d8 + d10 + d11 + d13 + d16 + d18 + d19 + d21 + d24 + d26
      + d27 + d29 + d32 + d33 + d36 + d38 + d39 + d41 + d44 + d46 + d47 + d49
      + d52 + d54 + d55 + d57 + d60 + d62 + d63 >= 5
 R55: d2 + d3 + d5 + d8 + d10 + d11 + d13 + d16 + d17 + d20 + d22 + d23 + d25
      + d28 + d30 + d31 + d34 + d35 + d37 + d40 + d42 + d43 + d45 + d48 + d49
      + d52 + d54 + d55 + d57 + d60 + d62 + d63 >= 5
 R56: d2 + d3 + d5 + d8 + d9 + d12 + d14 + d15 + d18 + d19 + d21 + d24 + d25
      + d28 + d30 + d31 + d34 + d35 + d37 + d40 + d41 + d44 + d46 + d47 + d50
      + d51 + d53 + d56 + d57 + d60 + d62 + d63 >= 5
 R57: d1 + d2 + d7 + d8 + d11 + d12 + d13 + d14 + d19 + d20 + d21 + d22 + d25
      + d26 + d31 + d32 + d35 + d36 + d37 + d38 + d41 + d42 + d47 + d48 + d49
      + d50 + d55 + d56 + d59 + d60 + d61 + d62 >= 4
 R58: d1 + d3 + d6 + d8 + d10 + d12 + d13 + d15 + d18 + d20 + d21 + d23 + d25
      + d27 + d30 + d32 + d34 + d36 + d37 + d39 + d41 + d43 + d46 + d48 + d49
      + d51 + d54 + d56 + d58 + d60 + d61 + d63 >= 4
 R59: d1 + d4 + d5 + d8 + d10 + d11 + d14 + d15 + d18 + d19 + d22 + d23 + d25
      + d28 + d29 + d32 + d34 + d35 + d38 + d39 + d41 + d44 + d45 + d48 + d49
      + d52 + d53 + d56 + d58 + d59 + d62 + d63 >= 4
 R60: d1 + d4 + d6 + d7 + d9 + d12 + d14 + d15 + d18 + d19 + d21 + d24 + d26
      + d27 + d29 + d32 + d34 + d35 + d37 + d40 + d42 + d43 + d45 + d48 + d49
      + d52 + d54 + d55 + d57 + d60 + d62 + d63 >= 4
 R61: d1 + d4 + d6 + d7 + d10 + d11 + d13 + d16 + d17 + d20 + d22 + d23 + d26
      + d27 + d29 + d32 + d34 + d35 + d37 + d40 + d41 + d44 + d46 + d47 + d50
      + d51 + d53 + d56 + d57 + d60 + d62 + d63 >= 4
 R62: d1 + d4 + d6 + d7 + d10 + d11 + d13 + d16 + d18 + d19 + d21 + d24 + d25
      + d28 + d30 + d31 + d33 + d36 + d38 + d39 + d42 + d43 + d45 + d48 + d50
      + d51 + d53 + d56 + d57 + d60 + d62 + d63 >= 4
 R63: d2 + d3 + d5 + d8 + d9 + d12 + d14 + d15 + d17 + d20 + d22 + d23 + d26
      + d27 + d29 + d32 + d33 + d36 + d38 + d39 + d42 + d43 + d45 + d48 + d50
      + d51 + d53 + d56 + d57 + d60 + d62 + d63 >= 3
Bounds
 0 <= d1 <= 1
 0 <= d2 <= 1
 0 <= d3 <= 1
 0 <= d4 <= 1
 0 <= d5 <= 1
 0 <= d6 <= 1
 0 <= d7 <= 1
 0 <= d8 <= 1
 0 <= d9 <= 1
 0 <= d10 <= 1
 0 <= d11 <= 1
 0 <= d12 <= 1
 0 <= d13 <= 1
 0 <= d14 <= 1
 0 <= d15 <= 1
 0 <= d16 <= 1
 0 <= d17 <= 1
 0 <= d18 <= 1
 0 <= d19 <= 1
 0 <= d20 <= 1
 0 <= d21 <= 1
 0 <= d22 <= 1
 0 <= d23 <= 1
 0 <= d24 <= 1
 0 <= d25 <= 1
 0 <= d26 <= 1
 0 <= d27 <= 1
 0 <= d28 <= 1
 0 <= d29 <= 1
 0 <= d30 <= 1
 0 <= d31 <= 1
 0 <= d32 <= 1
 0 <= d33 <= 1
 0 <= d34 <= 1
 0 <= d35 <= 1
 0 <= d36 <= 1
 0 <= d37 <= 1
 0 <= d38 <= 1
 0 <= d39 <= 1
 0 <= d40 <= 1
 0 <= d41 <= 1
 0 <= d42 <= 1
 0 <= d43 <= 1
 0 <= d44 <= 1
 0 <= d45 <= 1
 0 <= d46 <= 1
 0 <= d47 <= 1
 0 <= d48 <= 1
 0 <= d49 <= 1
 0 <= d50 <= 1
 0 <= d51 <= 1
 0 <= d52 <= 1
 0 <= d53 <= 1
 0 <= d54 <= 1
 0 <= d55 <= 1
 0 <= d56 <= 1
 0 <= d57 <= 1
 0 <= d58 <= 1
 0 <= d59 <= 1
 0 <= d60 <= 1
 0 <= d61 <= 1
 0 <= d62 <= 1
 0 <= d63 <= 1
Binaries
 d1  d2  d3  d4  d5  d6  d7  d8  d9  d10  d11  d12  d13  d14  d15  d16  d17 
 d18  d19  d20  d21  d22  d23  d24  d25  d26  d27  d28  d29  d30  d31  d32 
 d33  d34  d35  d36  d37  d38  d39  d40  d41  d42  d43  d44  d45  d46  d47 
 d48  d49  d50  d51  d52  d53 
End
