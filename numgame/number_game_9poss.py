import numpy as np

p = 2

for N in range(19):
    n = np.arange(1,10)
    
    for i in n:
        #produce number
        s = ''
        s = s + str(i)
        r1 = 0
        r2 = 0
        
        for j in range(N-1):
            r1 = (p*int(s[-1]))//10
            s = s + str(np.mod(p*int(s[-1])+r2,10))
            r2 = r1
    
        s = str(int(s[-1::-1]))
        S = int(s[-1]+s[:-1])
        s = int(s)
        #print(s)
        sol = (S == p*s)
        if sol:
            print(s)
         
    
    
    