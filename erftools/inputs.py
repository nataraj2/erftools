from collections.abc import MutableMapping
import numpy as np


class ERFInputFile(MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.verbose = kwargs.pop('verbose',True)
        self.store = dict()
        # SET DEFAULTS HERE
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __str__(self):
        #return '\n'.join([f'{key} = {str(val)}' for key,val in self.items()])
        s = ''
        for key,val in self.items():
            if isinstance(val, (list,tuple,np.ndarray)):
                val = ' '.join([str(v) for v in val])
            s += f'{key} = {str(val)}\n'
        return s.rstrip()

    def __getitem__(self, key):
        return self.store[self._keytransform(key)]

    def __setitem__(self, key, value):
        try:
            oldval = self[key]
        except KeyError:
            pass
        else:
            if self.verbose:
                print(f'Overriding existing `{key}` with {val}')
        finally:
            self.store[self._keytransform(key)] = value

    def __delitem__(self, key):
        del self.store[self._keytransform(key)]

    def __iter__(self):
        return iter(self.store)
    
    def __len__(self):
        return len(self.store)

    def _keytransform(self, key):
        return key
    
