import numpy as np

def foo(func, arr):
    return func(arr)

bar = np.eye(2)

print(foo(np.any, bar))
