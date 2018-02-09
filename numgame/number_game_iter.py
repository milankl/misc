## Number game

import numpy as np

N = 10**19
i = 0
sol = False
lim = 10**20
ii = 0

while i < lim-N and not(sol):
    s = str(N+i)
    S = int(s[-1]+s[:-1])
    sol = (S == 2*(N+i))
    #if np.mod(i,(lim-N)/100) == 0:
    #    ii = ii+1
    #    print(str(ii)+'% done.')
    i = i+1

if sol:
    print(N+i)
else:
    print(sol)