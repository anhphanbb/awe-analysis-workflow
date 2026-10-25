"""Copied from local juwavelet/utils.py; only change: stop after requested cluster count.
Original algorithm and peak-amplitude ranking retained.
"""
import numpy as np
import itertools
def log(*args, **kwargs):
    pass

def identify_cluster2d(decomposition, min_amp=2, thr=0.5, max_clusters=4):
    """
    Identify clusters of high magnitude

    Parameters
    ----------
    decomposition : dictionary
    min_amp : float
        smallest amplitude to consider a new cluster
    thr : float
        smallest amplitude to gather up into cluster

    Returns
    -------
    list
        maximal amplitude of cluster
    ndarray
        array with number of cluster at its position in the
        wavelet decomposition
    """

    adec = abs(decomposition["decomposition"])
    shape = adec.shape
    iwave = np.full_like(adec, -1, dtype=np.int16)
    amps = []
    idxs = []

    sel = np.zeros(shape, dtype=bool)

    sel1 = np.zeros((shape[0], shape[1]), dtype=bool)
    new1 = np.zeros((shape[0], shape[1]), dtype=bool)

    sel2 = np.zeros((shape[2], shape[3]), dtype=bool)
    new2 = np.zeros((shape[2], shape[3]), dtype=bool)

    for j in range(max_clusters):
        wmax = adec.max()
        if wmax < min_amp:
            break
        amps.append(wmax)
        udx = np.unravel_index(adec.argmax(), adec.shape)
        idxs.append(udx)
        log("------------")
        log("cluster=", j)
        log("amplitude=", wmax, "idx=", udx)
        sel[:] = False
        sel[udx] = True
        last_max = wmax

        awave2 = adec[:, :, udx[2], udx[3]]
        min_sel = awave2 > thr / 2
        sel1[:] = False
        sel1[udx[0], udx[1]] = True
        for i in itertools.count():
            new1[:] = False
            new1[1:, :] |= sel1[:-1, :]
            new1[:-1, :] |= sel1[1:, :]
            new1[:, 1:] |= sel1[:, :-1]
            new1[:, :-1] |= sel1[:, 1:]
            new1 &= min_sel
            new1 &= ~sel1
            new1 &= awave2 <= last_max
            if new1.sum() == 0:
                break
            new_last_max = awave2[new1].max()
            assert new_last_max <= last_max
            last_max = new_last_max
            if new1.sum() == 0:
                break
            sel1 |= new1

        for ix, iy in zip(*np.where(sel1)):
            sel[ix, iy, udx[2], udx[3]] = True
        for ix, iy in zip(*np.where(sel1)):
            sel2[:] = sel[ix, iy]
            awave2 = adec[ix, iy]
            min_sel = awave2 > thr
            last_max = awave2.max()
            for i in itertools.count():
                new2[:] = sel2
                new2[1:, :] |= sel2[:-1, :]
                new2[:-1, :] |= sel2[1:, :]
                new2[:, 1:] |= sel2[:, :-1]
                new2[:, :-1] |= sel2[:, 1:]
                new2 &= min_sel
                last_max = wmax
                new2 &= awave2 <= last_max
                new2 &= ~sel2
                if new2.sum() == 0:
                    break
                new_last_max = awave2[new2].max()
                assert new_last_max <= last_max
                last_max = new_last_max
                sel2 |= new2
            sel[ix, iy, :, :] = sel2

        log("selected elements:", sel.sum())
        adec[sel] = 0
        iwave[sel] = j

    return amps, idxs, iwave

