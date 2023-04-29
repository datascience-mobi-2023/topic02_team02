import numpy as np


def fib(x: int) -> int:
    # primitive implementation der Fibonacci Folge -> einfache Rekursion
    if x == 1:
        return 1
    elif x <= 1:
        return 1
    else:
        return x + fib(x-1)

def nums(n):
    yield n
    yield from nums(n+1)

s = nums(2)
def sieb(s):
    n = next(s)
    yield n
    yield from sieb(i for i in s if i%n!=0)

# generiert die ersten x Primzahlen durch das Sieb des Eratostenes.
def primzahlen(x: int) -> list:
    res = []
    prime = sieb(nums(2))
    for i in range(x):
        res.append(next(prime))
    return res