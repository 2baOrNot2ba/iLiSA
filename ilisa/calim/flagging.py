"""Provides flagging functionality for visibilities"""
import sys
import numpy as np
from numpy import ma as ma
import ilisa.operations.data_io as data_io
from ilisa.calim.visibilities import baseline_distances, rot2uv

class Flags:
    def __init__(self, vis=None, nrelems=None):
        if vis is not None:
            self.vis = vis
            self.shape = self.vis.shape
        if nrelems is not None:
            self.shape = (nrelems, nrelems)
        self.bl_mask = None
        self.pol_mask = None

    def select_cov_mask(self, selections):
        """
        From a covariance selection specification, make a mask matrix

        This function can be used for covariances such as baselines or
        polarizations.

        Parameters
        ----------
        selections: list
            Covariance selection specification.
            Integeters in list select corresponding elements.
            Tuples of length 1 select auto-correlation elements.
            Tuples of length 2, where 2nd item is None, select all covariance
            between all elements in 1st item slot.
            Tuples of length 2, select all covariances with one element from 1st
            list and the other element from 2nd list.
            If first element is None this means the final result is inverted,
            so rather than being a selection list this argument is interpreted as a
            deselection list (starting with all selected).
        nr_cov_el: int
            Number of covariance elements

        Examples
        --------
        >>> from ilisa.calim.flagging import Flags
        Select antenna 2 and 3:
        >>> Flags(nrelems=4).select_cov_mask([2,3]).bl_mask
        array([[False, False,  True,  True],
           [False, False,  True,  True],
           [ True,  True,  True,  True],
           [ True,  True,  True,  True]])
        Select auto-correlation of antenna 1:
        >>> Flags(nrelems=4).select_cov_mask([(1,)]).bl_mask
        array([[False, False, False, False],
           [False,  True, False, False],
           [False, False, False, False],
           [False, False, False, False]])
        Select baseline 0-3:
        >>> Flags(nrelems=4).select_cov_mask([(0,3)]).bl_mask
        array([[False, False, False,  True],
           [False, False, False, False],
           [False, False, False, False],
           [ True, False, False, False]])
        Select all auto-correlations:
        >>> Flags(nrelems=4).select_cov_mask([(None,)]).bl_mask
        array([[ True, False, False, False],
               [False,  True, False, False],
               [False, False,  True, False],
               [False, False, False,  True]])
        >>> Flags(nrelems=4).select_cov_mask([None, 0]).bl_mask
        array([[ True,  True,  True,  True],
               [ True, False, False, False],
               [ True, False, False, False],
               [ True, False, False, False]])
        """
        _invert = False
        maskmat = np.zeros(self.shape[-2:], dtype=bool)
        if len(selections)>0 and selections[0] is None:
            _invert = True
            selections.pop(0)
        for sel in selections:
            if type(sel) is int:
                # Antenna select
                maskmat[sel, :] = True
                maskmat[:, sel] = True
            elif type(sel) is tuple:
                if len(sel) == 1:
                    if sel[0] is None:
                        sel =  (range(self.shape[-1]),)
                    # Autocorrelations
                    maskmat[sel[0], sel[0]] = True
                elif len(sel) == 2:
                    if sel[1] is None:
                        if sel[0] is not None:
                            _sel = tuple(np.meshgrid(sel[0], sel[0]))
                            maskmat[_sel] = True
                        else:
                            maskmat[:, :] = True
                    else:
                        maskmat[sel[0], sel[1]] = True
                        maskmat[sel[1], sel[0]] = True
        if _invert:
            maskmat = np.logical_not(maskmat)
        self.bl_mask = maskmat
        return self

    def apply_vispol_flags(self):
        """
        Apply flags to visibilities

        Returns
        -------
        vispol_flagged: array
            Flagged polarized visibilities.
        """
        flagbls = self.bl_mask
        flagpols = self.pol_mask
        if flagbls is None and flagpols is None:
            vispol_flagged = ma.array(self.vis)
            return vispol_flagged
        if flagpols is None:
            polshape = self.vis.shape[:-2]
            flagpols = np.zeros(polshape, dtype=bool)
        else:
            # flagbls is None, so set to unselected
            blshape = self.vis.shape[-2:]
            flagbls = np.zeros(blshape, dtype=bool)
        flagpolmat = np.expand_dims(flagpols, axis=(-1, -2))
        flagblspolmat = flagbls + flagpolmat
        vispol_flagged = ma.array(self.vis, mask=flagblspolmat)
        return vispol_flagged

    def set_blflagargs(self, blflagargs):
        if type(blflagargs) == str:
            blflagargs = eval(blflagargs)
        self.select_cov_mask(blflagargs)
        return self

    def apply_blflagargs(self, blflagargs):
        self.set_blflagargs(blflagargs)
        vis_pol = self.apply_vispol_flags()
        return vis_pol

    def save(self, filename='flags_bl.txt'):
        np.savetxt(filename, self.bl_mask, fmt='%i')

    def load(self, filename='flags_bl.txt'):
        self.bl_mask = np.loadtxt(filename, dtype=bool)
        self.shape = self.bl_mask.shape
        return self


def main_cli():
    flagtype = sys.argv[1]
    dataff = sys.argv[3]
    cvcobj = data_io.CVCfiles(dataff)
    fv = Flags(nrelems=cvcobj.cvcdim1//2)
    if flagtype == 'manual':
        sel = eval(sys.argv[2])
        fv.set_blflagargs(sel)
    elif flagtype == 'shortbl':
        distance = float(sys.argv[2])
        uv_m = rot2uv(cvcobj.stn_antpos, cvcobj.stn_rot)
        bl_dist = baseline_distances(uv_m)
        short_bls = bl_dist < distance
        fv.bl_mask = short_bls
    else:
        print("flagtype {} not recognized".format(flagtype))
        sys.exit()
    fv.save(dataff+'/flags_{}.txt'.format(flagtype))


if __name__ == "__main__":
    main_cli()
