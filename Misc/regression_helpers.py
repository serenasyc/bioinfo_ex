import numpy as np

class rotate_data:
    '''
    circular rotation of 10 rows of data
    '''
    def __init__(self, x):
        self.x = x
        
    def __iter__(self):
        return self
    
    def next(self):
        part = self.x[:10,:]
        rest = self.x[10:,:]
        self.x = np.concatenate((rest, part))
        return rest,part

def polybasis(x_arr, deg):
    '''
    x_arr   ARRAY  1xD
    deg     INT    degree to raise each element of x to 
    '''
    return np.power(x_arr, deg)
    
def sigbasis(x_arr, mu, s):
    '''
    x_arr  ARRAY  1xD
    mu,s   INT    params of sigmoidal function
    '''
    return np.array([1/(1+np.exp((mu-x)/s)) for x in x_arr])
    
def get_designmatrix_poly(xmtx, degree):
    '''
    xmtx  MATRIX  N by D matrix (N samples of D features)
    M     INT     degree of polynomial
    return the design matrix
    '''
    (N,D) = xmtx.shape
    # initialize matrix 
    dmtx = np.ones((1, degree*D+1))
    # for each row in the design matrix
    for row in range(N):
        # for each degree from 1 to M
        arr = np.array([1])
        for deg in range(1,degree+1):  
            # calculate the terms of x to the power of i
            arr = np.append(arr, polybasis(xmtx[row,:], deg))
        dmtx = np.vstack((dmtx, arr))
    #remove the extra row
    return dmtx[1:,]
    
def get_designmatrix_sig(xmtx, mu_values, s):
    '''
    xmtx  MATRIX  N by D matrix (N samples of D features)
    M     INT     degree of polynomial
    return the design matrix
    '''
    (N,D) = xmtx.shape
    # initialize matrix 
    dmtx = np.ones((1, D*2+1))
    # for each row in the design matrix
    for row in range(N):
        arr = np.array([1])
        for mu in mu_values:
            arr = np.append(arr, sigbasis(xmtx[row,:], mu, s))
        dmtx = np.vstack((dmtx,arr))
    #remove the extra row
    return dmtx[1:,]

def polyfunc(xmtx, w, degree):
    '''
    xmtx    MATRIX  input N rows of D features
    w       MATRIX  Mx1 weights
    degree  INT     degree of polynomial basis
    returns the estimates of t = y(x,w) = w^T*dmtx
    '''
    dmtx = get_designmatrix_poly(xmtx, degree)
    return np.dot(w.T, dmtx.T)
    
def sigfunc(xmtx, w, mu_values, s):
    dmtx = get_designmatrix_sig(xmtx, mu_values, s)
    return np.dot(w.T, dmtx.T)

def rms(w, x, t, func, params):
    '''
    w       MATRIX      1xM weights 
    x       MATRIX      NxD inputs
    t       MATRIX      Nx1 target values
    func    FUNCTION    used to calculate estimates of y
    params  LIST        parameters for func
    return root mean square error
    '''
    y = func(x, w, *params) 
    # sum of squared error
    sse = np.sum(np.power(y.T-t, 2))
    return np.sqrt(sse/len(t))    

def get_MLweights(xmtx, tvec, basisfn, params):
    '''
    xmtx    MATRIX  N rows of input arrays, each with D features
    tvec    ARRAY   M rows of target values
    basisfn STR     type of basis function used
    param   LIST    contains params for basis function
    return maximum likelihood weights
    calculated via Moore-Penrose pseudo-inverse & sum of squared error
    '''
    if basisfn == 'poly_reg':
        degree,lam = params
        dmtx = get_designmatrix_poly(xmtx, degree) 
        # if lam == 0, identity matrix not added, continues at bottom
        if lam != 0:
            wML  = np.dot(dmtx.T, dmtx) + lam*np.identity(dmtx.shape[1])
            # only invertible if regularized
            wML = np.linalg.inv(wML) 
            wML = np.dot(wML,dmtx.T)*tvec
            return wML
    elif basisfn == 'polynomial':
        degree = params[0]
        dmtx = get_designmatrix_poly(xmtx, degree)
    elif basisfn == 'sigmoidal':
        mu_values,s = params
        dmtx = get_designmatrix_sig(xmtx, mu_values, s)
    wML  = np.dot(np.linalg.pinv(dmtx), tvec)
    return wML 
