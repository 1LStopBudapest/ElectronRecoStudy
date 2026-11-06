import numpy as np





class Infix:
    def __init__(self, function):
        self.function = function
    def __ror__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __or__(self, other):
        return self.function(other)
    def __rlshift__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __rshift__(self, other):
        return self.function(other)
    def __call__(self, value1, value2):
        return self.function(value1, value2)

    
def qsum(*args):
    sm = 0
    for arg in args:
        sm += arg**2
    return np.sqrt(sm)
qs = Infix(qsum)





def overlapIndicesOfLists(a, b):
    from bisect import bisect_left, bisect_right

    overlap_start = max(a[0], b[0])
    overlap_end   = min(a[-1], b[-1])

    a_i_start = bisect_left(a, overlap_start)
    a_i_end   = bisect_right(a, overlap_end) - 1

    b_i_start = bisect_left(b, overlap_start)
    b_i_end   = bisect_right(b, overlap_end) - 1

    return a_i_start, a_i_end, b_i_start, b_i_end




def trapcalc(x,y):
    sum = 0
    for i in range(len(x)-1):
        if(y[i] < y[i+1]):
            sum += y[i] * (x[i+1]-x[i])  +  (y[i+1]-y[i]) * (x[i+1]-x[i])/2
        else:
            sum += y[i+1] * (x[i+1]-x[i])  +  (y[i]-y[i+1]) * (x[i+1]-x[i])/2
    return sum

def useBigger(a,b):
    return a if a >= b else b

def useSmaller(a,b):
    return a if a <= b else b