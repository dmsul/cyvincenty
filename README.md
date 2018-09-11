# cy-vincenty

`cyvincenty` is a simple [Cython](https://github.com/cython/cython)
implementation of Vincenty's inverse formula to find the distance between two
latitude-longitude points.

This package contains two functions, `vincenty` and `vincenty_cross`.

`vincenty` finds the distance between two points:

```python
from cyvincenty import vincenty

x1, y1 = -118, 32
x2, y2 = -117, 31

dist_km = vincenty(x1, y1, x2, y2)
```


`vincenty_cross` finds the all the distances between two sets of points (like
`scipy.spatial.distance.cdist`).

```python
import numpy as np
from cyvincenty import vincenty_cross

x1 = np.linspace(-118, -117, 100, dtype=np.float32)
y1 = np.linspace(32, 34, 100, dtype=np.float32)

x2 = np.linspace(-112, -110, 35, dtype=np.float32)
y2 = np.linspace(35, 37, 35, dtype=np.float32)

# A 100-by-35 numpy array
dist_km = vincenty_cross(x1, y1, x2, y2)
```

# Requirements

- Python 3 (tested with 3.6)
- Cython (tested with 0.27)
- [NumPy](http://www.numpy.org/) (tested with 1.13)
- C compiler (tested with [Visual Studio Community 2017](https://www.visualstudio.com/downloads/))

# Installation

Just clone the repository and run

```console
$ python setup.py install
```

If you're using Windows, you will need to jump through the necessary hoops so
Python can find your C compiler. Installing Visual Studio is the most
straightforward route, but make sure you install the C/C++ command line
interface (CLI) tools! If you're using Visual Studio, the easiest route is to
run `setup.py` in a Visual Studio Developer console.
