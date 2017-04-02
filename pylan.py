## pylan.py - own functions

def rmean(x, N):
    """running mean / moving average
    
    y = rmean(x,N) returns a vector y that represents the running mean of vector x. Parameter N should be odd and defines the window length. The first as well as the last (N-1)/2 are replaced by NaNs."""
    
    import numpy as np
    
    a = (N-1)/2
    y = np.convolve(x, np.ones((N,))/N)[a:len(x)+a]
    y[:a] = None
    y[-a:] = None
    return y


def cca(x,y):
    
    """ canonical correlation analysis cca
    
    wx, wy, r = cca(x,y) returns wx, wy two matrices which columns [:,i] correspond to the canonical weights (normalized eigenvectors) and a vector r containing the canonical correlations, all sorted in decreasing order. cca assumes as input matrices x,y of size l*m (time*nvar), and l*n, that are centered (no mean along 1st axis) within the function. cca returns an error if either x,y are not full rank."""
    
    import numpy as np
    
    mx = x.shape[1]
    my = y.shape[1]
    l = x.shape[0] #needs to be the same for y
    if l != y.shape[0]:
        raise ValueError('Time dimension is not same length for x,y')
    xrank = np.linalg.matrix_rank(x)
    yrank = np.linalg.matrix_rank(y)
    
    if mx > xrank:
        raise ValueError('Matrix x is not full rank.')
    if my > yrank:
        raise ValueError("Matrix y is not full rank.")
    
    #no mean
    x = x - np.outer(x.mean(axis=0),np.ones(l)).transpose()
    y = y - np.outer(y.mean(axis=0),np.ones(l)).transpose()
    
    #covariance estimators
    Sxy = np.dot(x.transpose(),y) / l
    Sxx = np.dot(x.transpose(),x) / l
    Syy = np.dot(y.transpose(),y) / l

    B1 = np.dot(np.linalg.inv(Sxx),Sxy)
    B2 = np.dot(np.linalg.inv(Syy),Sxy.transpose())

    evalx, eigvx = np.linalg.eig(np.dot(B1,B2))
    evaly, eigvy = np.linalg.eig(np.dot(B2,B1))
    
    #normalize eigenvectors
    eigvx = eigvx / np.outer(np.ones((mx,1)),np.sqrt((eigvx**2).sum(axis=0)))
    eigvy = eigvy / np.outer(np.ones((my,1)),np.sqrt((eigvy**2).sum(axis=0)))

    # eigenvalues should be the same in evalx and evaly
    rx = np.sqrt(abs(evalx)) #correlation
    ry = np.sqrt(abs(evaly))
    
    #sort
    ordargx = np.argsort(rx)[-1:-mx-1:-1] #decreasing order
    ordargy = np.argsort(ry)[-1:-mx-1:-1]
    rx = rx[ordargx]
    ry = ry[ordargy]
    eigvx = eigvx[:,ordargx]
    eigvy = eigvy[:,ordargy]
    
    if mx >= my:
        r = rx
    else:
        r = ry

    return eigvx, eigvy, r
    
    
def rmse(x,y):
    
    """ root mean square error (RMSE)
    
    e = rmse(x,y) calculates the root mean square error of two vectors x,y of equal length"""
    
    import numpy as np
    
    e = np.sqrt(((x-y)**2).mean())
    return e
    
def nans(s, dtype=float):
    
    """ initialize a numpy array of shape s (=tuple) with all entries being nans. """
    
    import numpy as np
    
    a = np.empty(s, dtype)
    a.fill(np.nan)
    return a
    
def acf(x,l):
    
    """ Auto correlation function of a given vector x with a maximum lag of length l.
    """
    
    import numpy as np
    
    return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1] \
        for i in range(1, l)])
        
def detrend(t,x):

    """ Detrends the timeseries x with polyfit/polyval. The time step are either mono spaced or given in vector t. """
    
    import numpy as np
    return x - np.polyval(np.polyfit(t,x,1),t)
    
def whiten(t,x,n):

    """ Whitens the timeseries x with polyfit/polyval. The time step are either mono spaced or given in vector t. """
    
    import numpy as np
    return x - np.polyval(np.polyfit(t,x,n),t)

def closest(x,y):
    """ Finds the location of a value in vector x that is from all entries of x closest to a given value y. Returns the position as a boolean vector of same size as x. """
    
    import numpy as np
    
    return np.argmin(abs(x-y))
    
def round1(x):
    """ rounds floats to the most significant digit."""
    import numpy as np
    return np.round(x, -np.int(np.floor(np.log10(abs(x)))))
    
def findzc(x,a):
    """ finds crossing of vector x with value a, several values might be returned."""
    import numpy as np
    return np.where(abs(np.diff(np.sign(x-a))) > 0)[0]
    
    
## colormaps

exec(open('python/modules/colormaps.py').read())

def cmaps():
    import numpy as np
 
    cmat1 = [[0.,0.,0.7], [0.,0.,0.9], [0., .4, 1.], [0., .6, 1.], [0., .7, 1.], [0., .8, 1.], [0.5, .9, 1], [1., 1., 1.], [1., .9, 0.], [1., .8, 0.], [1., .7, 0.], [1., .6, 0.], [1., .4, 0.], [0.9, 0., 0.], [0.7, 0., 0.]]
 
    cmat2 = np.array([230,0,0, \
    230,3,1, \
    231,6,3, \
    231,10,4, \
    231,13,5, \
    231,16,7, \
    231,19,8, \
    231,22,10, \
    232,26,11, \
    232,29,12, \
    232,32,14, \
    232,35,15, \
    232,38,16, \
    232,42,18, \
    233,45,19, \
    233,48,21, \
    233,51,22, \
    233,54,23, \
    233,58,25, \
    233,61,26, \
    234,64,27, \
    234,67,29, \
    234,70,30, \
    234,74,32, \
    234,77,33, \
    235,80,34, \
    235,83,36, \
    235,86,37, \
    235,90,38, \
    235,93,40, \
    235,96,41, \
    236,99,42, \
    236,102,44, \
    236,106,45, \
    236,109,47, \
    236,112,48, \
    236,115,49, \
    237,118,51, \
    237,122,52, \
    237,125,53, \
    237,128,55, \
    237,131,56, \
    237,134,58, \
    238,138,59, \
    238,141,60, \
    238,144,62, \
    238,147,63, \
    238,150,64, \
    238,154,66, \
    239,157,67, \
    239,160,68, \
    239,163,70, \
    239,166,71, \
    239,170,73, \
    239,173,74, \
    240,176,75, \
    240,179,77, \
    240,182,78, \
    240,186,79, \
    240,189,81, \
    240,192,82, \
    241,195,84, \
    241,198,85, \
    241,202,86, \
    241,204,88, \
    241,205,91, \
    242,206,93, \
    242,207,96, \
    242,207,99, \
    242,208,101, \
    243,209,104, \
    243,210,106, \
    243,211,109, \
    243,212,112, \
    243,212,114, \
    244,213,117, \
    244,214,120, \
    244,215,122, \
    244,216,125, \
    245,216,128, \
    245,217,130, \
    245,218,133, \
    245,219,136, \
    246,220,138, \
    246,220,141, \
    246,221,143, \
    246,222,146, \
    246,223,149, \
    247,224,151, \
    247,225,154, \
    247,225,157, \
    247,226,159, \
    248,227,162, \
    248,228,165, \
    248,229,167, \
    248,229,170, \
    249,230,173, \
    249,231,175, \
    249,232,178, \
    249,233,180, \
    249,233,183, \
    250,234,186, \
    250,235,188, \
    250,236,191, \
    250,237,194, \
    251,237,196, \
    251,238,199, \
    251,239,202, \
    251,240,204, \
    252,241,207, \
    252,242,210, \
    252,242,212, \
    252,243,215, \
    252,244,217, \
    253,245,220, \
    253,246,223, \
    253,246,225, \
    253,247,228, \
    254,248,231, \
    254,249,233, \
    254,250,236, \
    254,250,239, \
    255,251,241, \
    255,252,244, \
    255,253,247, \
    255,254,249, \
    255,255,252, \
    255,255,254, \
    254,255,255, \
    252,255,255, \
    249,254,255, \
    246,253,254, \
    243,252,254, \
    241,251,253, \
    238,251,253, \
    235,250,253, \
    232,249,252, \
    230,248,252, \
    227,247,251, \
    224,247,251, \
    221,246,250, \
    219,245,250, \
    216,244,250, \
    213,243,249, \
    210,243,249, \
    208,242,248, \
    205,241,248, \
    202,240,248, \
    199,239,247, \
    197,239,247, \
    194,238,246, \
    191,237,246, \
    188,236,245, \
    186,235,245, \
    183,235,245, \
    180,234,244, \
    177,233,244, \
    175,232,243, \
    172,231,243, \
    169,231,242, \
    166,230,242, \
    164,229,242, \
    161,228,241, \
    158,227,241, \
    155,227,240, \
    153,226,240, \
    150,225,239, \
    147,224,239, \
    144,223,239, \
    142,223,238, \
    139,222,238, \
    136,221,237, \
    133,220,237, \
    131,219,236, \
    128,219,236, \
    125,218,236, \
    122,217,235, \
    120,216,235, \
    117,215,234, \
    114,215,234, \
    111,214,234, \
    109,213,233, \
    106,212,233, \
    103,211,232, \
    100,211,232, \
    98,210,231, \
    95,209,231, \
    92,208,231, \
    89,207,230, \
    87,207,230, \
    84,206,229, \
    81,205,229, \
    79,202,228, \
    78,199,228, \
    77,196,228, \
    76,193,227, \
    74,190,227, \
    73,186,226, \
    72,183,226, \
    71,180,226, \
    69,177,225, \
    68,173,225, \
    67,170,224, \
    66,167,224, \
    64,164,224, \
    63,161,223, \
    62,157,223, \
    61,154,222, \
    59,151,222, \
    58,148,222, \
    57,145,221, \
    55,141,221, \
    54,138,220, \
    53,135,220, \
    52,132,220, \
    50,128,219, \
    49,125,219, \
    48,122,218, \
    47,119,218, \
    45,116,218, \
    44,112,217, \
    43,109,217, \
    42,106,216, \
    40,103,216, \
    39,100,216, \
    38,96,215, \
    37,93,215, \
    35,90,214, \
    34,87,214, \
    33,84,214, \
    32,80,213, \
    30,77,213, \
    29,74,212, \
    28,71,212, \
    26,67,212, \
    25,64,211, \
    24,61,211, \
    23,58,210, \
    21,55,210, \
    20,51,210, \
    19,48,209, \
    18,45,209, \
    16,42,208, \
    15,39,208, \
    14,35,208, \
    13,32,207, \
    11,29,207, \
    10,26,206, \
    9,22,206, \
    8,19,206, \
    6,16,205, \
    5,13,205, \
    4,10,204, \
    3,6,204, \
    1,3,204, \
    0,0,203]).reshape((3,-1))

    return cmat1,cmat2
