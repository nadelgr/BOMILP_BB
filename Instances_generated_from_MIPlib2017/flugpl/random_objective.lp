\ENCODING=ISO-8859-1
\Problem name: flugpl.mps

Maximize
 KOSTEN: -686 STM1             + 2186 ANM1             -2251 UE1              + 888 STM2             -2280 ANM2             -2207 UE2              + 456 STM3             + 49 ANM3             + 1223 UE3              + 1113 STM4             + 2030 ANM4             -898 UE4              -353 STM5             + 1871 ANM5             -1287 UE5              + 1673 STM6             + 2596 ANM6             + 2328 UE6       
Subject To
 ANZ1: STM1  = 60
 STD1: 150 STM1 - 100 ANM1 + UE1 >= 8000
 UEB1: - 20 STM1 + UE1 <= 0
 ANZ2: 0.9 STM1 + ANM1 - STM2  = 0
 STD2: 150 STM2 - 100 ANM2 + UE2 >= 9000
 UEB2: - 20 STM2 + UE2 <= 0
 ANZ3: 0.9 STM2 + ANM2 - STM3  = 0
 STD3: 150 STM3 - 100 ANM3 + UE3 >= 8000
 UEB3: - 20 STM3 + UE3 <= 0
 ANZ4: 0.9 STM3 + ANM3 - STM4  = 0
 STD4: 150 STM4 - 100 ANM4 + UE4 >= 10000
 UEB4: - 20 STM4 + UE4 <= 0
 ANZ5: 0.9 STM4 + ANM4 - STM5  = 0
 STD5: 150 STM5 - 100 ANM5 + UE5 >= 9000
 UEB5: - 20 STM5 + UE5 <= 0
 ANZ6: 0.9 STM5 + ANM5 - STM6  = 0
 STD6: 150 STM6 - 100 ANM6 + UE6 >= 12000
 UEB6: - 20 STM6 + UE6 <= 0
Bounds
 0 <= ANM1 <= 18
 57 <= STM2 <= 75
 0 <= ANM2 <= 18
 57 <= STM3 <= 75
 0 <= ANM3 <= 18
 57 <= STM4 <= 75
 0 <= ANM4 <= 18
 57 <= STM5 <= 75
 0 <= ANM5 <= 18
 57 <= STM6 <= 75
 0 <= ANM6 <= 18
Generals
 ANM1  STM2  ANM2  STM3  ANM3  STM4  ANM4  STM5  ANM5  STM6  ANM6 
End
