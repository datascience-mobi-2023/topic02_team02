# hier werden alle inputs zusammgengefügt.
import packages

# bitte alle funktionen, die verwendet werden sollen einzeln so importieren, nicht ganze dateien importieren
from test_func import fib, primzahlen

if __name__ == '__main__':
    print(fib(12))

    # die ersten x primzahlen
    print(primzahlen(100))