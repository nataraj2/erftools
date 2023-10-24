import os
import glob
import numpy as np
import pandas as pd
import yt
yt.set_log_level('error')

def load_pltfile_column(dpath,
                        zlevels=None,
                        zlevels_stag=None,
                        return_amrex_dataset=False):
    """Extract column of data from (ilo,jlo)"""
    amrds = yt.load(dpath)
    grid_spacing = (amrds.domain_right_edge.value
                    - amrds.domain_left_edge.value) / amrds.domain_dimensions
    # get vertical levels
    if zlevels is not None:
        assert len(zlevels) == amrds.domain_dimensions[2]
    elif zlevels_stag is not None:
        assert len(zlevels_stag) == amrds.domain_dimensions[2]+1
        zlevels = 0.5 * (zlevels_stag[1:] + zlevels_stag[:-1])
    else:
        zlevels = 0.5*grid_spacing[2] \
                + np.arange(amrds.domain_dimensions[2]) * grid_spacing[2]
    # extract column
    col_left_edge = amrds.domain_left_edge.value
    col_right_edge = [grid_spacing[0], grid_spacing[1], amrds.domain_right_edge.value[2]]
    col = amrds.box(col_left_edge, col_right_edge)
    # get column values
    fields = [name for (dtype,name) in amrds.field_list if dtype=='boxlib']
    coldata = {field: col[field].value for field in fields}
    df = pd.DataFrame(coldata, index=pd.Index(zlevels, name='height'))
    if return_amrex_dataset:
        return df, amrds
    else:
        return df

class Column(object):
    """Helper for loading in pltfile data from a column simulation, generally
    set up as a domain that has
        amr.n_cell = amr.blocking_factor amr.blocking_factor n_levels
        geometry.is_periodic = 1 1 0

    The data are stored in a pandas dataframe with a time-height multiindex,
    accessible by:
        Column.df  # return Column._df (with multiindex) if multiple times were
                   # loaded, otherwise a dataframe with height index
        Column.t  # equivalent to Column._df.index.levels[0]
        Column.z  # equivalent to Column._df.index.levels[1]
        Column['x_velocity']  # shorthand for Column._df['x_velocity']
    """
    def __init__(self,pltfiles,**kwargs):
        """Read one or more pltfiles and convert to dataframe format

        Parameters
        ----------
        pltfiles : str, list
            May be specified as a globstring, a single pltfile path,
            or a list of pltfile paths
        """
        if isinstance(pltfiles, str):
            # first see if input is a globstring
            globbed = glob.glob(pltfiles)
            if len(globbed) > 0:
                pltfiles = globbed
            elif os.path.isdir(pltfiles):
                pltfiles = [pltfiles]
            else:
                raise FileNotFoundError(
                        f'{pltfiles} is neither a valid directory path'
                         ' or glob string')
        else:
            # input checks
            assert hasattr(pltfiles,'__iter__')
            assert all([os.path.isdir(pltfile) for pltfile in pltfiles])
        self.pltfiles = pltfiles
        self._sort_pltfiles()
        self._load_pltfile_wrapper(**kwargs)

    def _sort_pltfiles(self):
        if any(['.old.' in dpath for dpath in self.pltfiles]):
            print('Warning: One or more input files is old -- ignoring')
        pltfiles = [dpath for dpath in self.pltfiles
                    if os.path.isdir(dpath) and (not '.old.' in dpath)]
        pltnames = [os.path.split(dpath)[1] for dpath in pltfiles]
        assert all([name.startswith('plt') for name in pltnames])
        outputsteps = [int(name[3:]) for name in pltnames]
        sortorder = np.argsort(outputsteps)
        self.pltfiles = [pltfiles[i] for i in sortorder]

    def _load_pltfile_wrapper(self,**kwargs):
        dflist = []
        times = []
        for pltfile in self.pltfiles:
            print('\rLoading',pltfile,end='')
            df,amrds = load_pltfile_column(pltfile, return_amrex_dataset=True,
                                           **kwargs)
            dflist.append(df)
            times.append(amrds.current_time.value.item())
        print('')
        self.t = np.array(times)
        tindex = pd.Index(times,name='time')
        self._df = pd.concat(dflist, axis=0, keys=tindex)
        self.z = self._df.index.levels[1].values

    @property
    def df(self):
        if len(self.t) == 1:
            return self._df.xs(self.t[0],level='time')
        else:
            return self._df

    def __getitem__(self,key):
        return self.df[key]
