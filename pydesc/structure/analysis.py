# Copyright 2017 Giovanni Nico, Tymoteusz Oleniecki
#
# This file is part of PyDesc.
#
# PyDesc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyDesc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.

"""
Tools for structure and trajectory analysis.

created: 10.07.2013 - , Giovanni Nico, Tymoteusz 'hert' Oleniecki
"""

import numpy
from pydesc.mers import Residue

eps = 0.00000001
rtod = 180. / numpy.pi


def calc_R(vec):
    """Returns Re and Im parts for sum of given angles vector."""
    return numpy.sum(numpy.sin(vec)), numpy.sum(numpy.cos(vec))


def adi(vec):  # ADI
    """Calculates and returns ADI value.

    Argument:
    vec -- numpy.array of omega angles (sum of psi and phi) over trajectory.
    """
    sums, sumc = calc_R(vec)
    sumss, sumcc = calc_R(2 * vec)

    mod = numpy.sqrt((sums ** 2) + (sumc ** 2))
    modd = numpy.sqrt((sumss ** 2) + (sumcc ** 2))
    nval = len(vec)
    res = nval * (nval - modd) / (2 * (mod ** 2))

    return res


def adis(vec):
    """Calculates and returns ADI with sign.

    Argument:
    vec -- numpy.array of angles.
    """
    n = int(.05 * len(vec))
    r1 = calc_R(vec[:n])
    r2 = calc_R(vec[len(vec) - n:])
    return numpy.sign(numpy.cross(r1, r2)) * adi(vec)


def adicum(vec):  # ADI
    """Calculates cumulative ADI for vector of phi and psi angles sums."""

    sums = 0
    sumc = 0
    sumss = 0
    sumcc = 0
    sums1 = 0
    sumc1 = 0
    nlines = len(vec)

    for angle in vec:
        sums1 += numpy.sin(angle)
        sumc1 += numpy.cos(angle)
        angle = numpy.arctan2(sums1, sumc1)
        sums += numpy.sin(angle)
        sumc += numpy.cos(angle)
        sumss += numpy.sin(2 * angle)
        sumcc += numpy.cos(2 * angle)

    mod = numpy.sqrt(sums ** 2 + sumc ** 2)
    modd = numpy.sqrt(sumss ** 2 + sumcc ** 2)
    nval = nlines
    res = nval * (nval - modd) / (2 * mod ** 2)

    return res


def atip(vec):  # ATIP
    """Calculates ATIP index for vector of phi and psi angles sums."""

    nlines = len(vec)
    minframes = min(20, nlines * 0.25)
    numbseries = int(nlines / float(minframes))
    rangeold = 0

    # ~ import pdb; pdb.set_trace()

    for p in range(2, numbseries):
        numbstep = int(nlines / float(p))
        adimax = -1
        adimin = 9999
        for imin in (0, int(0.5 * numbstep)):
            imax = imin + numbstep
            for k in range(0, p):
                iadi = adi(vec[imin: imax])
                adimax = max(iadi, adimax)
                adimin = min(iadi, adimin)
                imin = imax
                imax = imin + numbstep
                if imax > nlines:
                    break

        range_adi = adimax - adimin

        if range_adi < rangeold:
            break
        rangeold = range_adi

    return rangeold


def atin(vec):  # ATIN
    """Calculates ATIN index for vector of phi and psi angles sums."""

    nlines = len(vec)
    adimin = 9999
    numbstep = 20
    numbstep1 = nlines * 0.25
    numbstep = min(numbstep, numbstep1)
    numbseries = int(nlines / float(numbstep))

    for itry in (0, 1):
        if itry == 0:
            imin = 0
        else:
            imin = int(0.5 * numbstep)

        imax = int(imin + numbstep)
        for k in range(0, numbseries):
            iadi = adi(vec[imin:imax])
            if iadi < adimin:
                adimin = iadi
            imin = imax
            imax = imin + numbstep
            if imax > nlines:
                break

    numbstep = 2
    numbseries = int(nlines / float(numbstep))

    for itry in (0, 1):
        if itry == 0:
            imin = 2
            sdadi1 = 0
            aveadi1 = 0
        else:
            imin = int(0.5 * numbstep)
            sdadi2 = 0
            aveadi2 = 0

        imax = int(imin + numbstep)
        for k in range(0, numbseries):
            iadi = adi(vec[imin:imax])
            if itry == 0:
                aveadi1 = aveadi1 + iadi
                sdadi1 = sdadi1 + iadi ** 2
            else:
                aveadi2 = aveadi2 + iadi
                sdadi2 = sdadi2 + iadi ** 2
            imin = imax
            imax = imin + numbstep
            if imax > nlines:
                break

        if itry == 0:
            aveadi1 = aveadi1 / float(numbseries)
            sdadi1 = sdadi1 / float(numbseries)
            sdadi1 = numpy.sqrt(sdadi1 - aveadi1 ** 2)
        else:
            aveadi2 = aveadi2 / float(numbseries)
            sdadi2 = sdadi2 / float(numbseries)
            sdadi2 = numpy.sqrt(sdadi2 - aveadi2 ** 2)

    sdadi = max(sdadi1, sdadi2)
    q = adimin

    m = 0.0014
    atin = (sdadi - q) / m

    return atin


def tai(adi1, adi2, ati1, ati2, mai, pai):  # TAI
    ndef = 3
    pai1 = numpy.arctan2(ati1, adi1)
    pai2 = numpy.arctan2(ati2, adi2)

    if ndef == 1:
        mai = numpy.sqrt(mai1 * mai2)
        if mai1 < eps:
            mai = mai2
        if mai2 < eps:
            mai = mai1
            pai = numpy.sqrt(pai1 * pai2)
    elif ndef == 2:
        sumwt = adi1 + adi2
        mai = (mai1 * adi1 + mai2 * adi2) / sumwt
        pai = numpy.sqrt(pai1 * pai2)
    elif ndef == 3:
        sums = mai1 * numpy.sin(pai1) + mai2 * numpy.sin(pai2)
        sumc = mai1 * numpy.cos(pai1) + mai2 * numpy.cos(pai2)
        mai = numpy.sqrt(numpy.pow(sums, 2) + numpy.pow(sumc, 2))
        pai = numpy.arctan2(sums, sumc)

    mai = rtod * numpy.arccos((1 - mai) / (1 + mai))
    pai = rtod * pai

    return atin


def pad(adi_i):
    """Returns PAD value.

    Argument:
    adi_i -- ADI calue.
    """
    return rtod * numpy.arccos((1 - adi_i) / (1 + adi_i))


def pai(atip_i, adi_i):
    """Returns PAI.

    Arguments:
    atip_i -- ATIP value.
    adi_i -- ADI value.
    """
    return rtod * numpy.arctan2(atip_i, adi_i)


def padcum(adicum_i):
    """Return PADcum value.

    Argument:
    adicum_i -- ADIcum value.
    """
    return rtod * numpy.arccos((1 - adicum_i) / (1 + adicum_i))


def tag(padcum, pai):
    """Returns tag for given values of PADcum and pai.

    Returns 'T', 't' or 'F' for transition, small transition and fluctuation respectively.
    """
    if padcum < 10:
        return "F"
    if pai < 30:
        return "F"
    elif pai > 60:
        return "t"
    return "T"


def corre(thetas):
    """

    Argument:
    thetas -- output from calculate_theta function.
    """
    lst = [thetas[k] for k in thetas]
    res = numpy.zeros(len(lst), len(lst))
    for i, mer in enumerate(lst):
        for k, mer2 in enumerate(lst[i + 1:]):
            res[i, k] = numpy.mean(numpy.exp((0 + 1j) * mer) * numpy.exp((0 - 1j) * mer2))

    return res


def corre1(vec, vec1):
    """

    Argument:
    vec and vec1 contains time series of angular values, with the same number of elements
    """
    return numpy.mean(numpy.exp((0 + 1j) * numpy.array(vec[1:])) * numpy.exp((0 - 1j) * numpy.array(vec1)))


def calculate_omega(stc):
    """Returns dictionary of all dihedral angles for each mer in subsequent trajectory frames.

    Argument:
    stc-- pydesc.structure.Structure linked to a trajectory.
    """

    mangs = {m: [] for m in stc}

    stc.frame = 0
    # ~ stc.frame = 1
    while True:
        Residue.calculate_angles_static(stc)

        for mer in stc:
            mangs[mer].append(sum(mer.angles))

        try:
            stc.next_frame()
        except IndexError:
            break

    for k, v in mangs.items():
        mangs[k] = numpy.array(v)

    return mangs


def calculate_theta(stc, mer_ind):
    """Returns theta angles for given mer over all frames of its structure trajectory."""
    mer = stc[mer_ind]

    res = {m: [] for m in stc}

    stc.frame = 1
    while True:
        rvec = (mer.ca - mer.cbx).get_unit_vector()
        for m in stc:
            if m == mer: continue
            vec = (m.ca - mer.ca).get_unit_vector()
            res[m].append(numpy.arccos(numpy.dot(vec, rvec)))
        try:
            stc.next_frame()
        except IndexError:
            break

    for k, v in res.items():
        res[k] = numpy.array(v)

    return res


def calculate_indices(angle_vec):
    """Returns array of indices for each mer of given structure.

    Argument:
    angle_vec -- list of angles for one mer over its trajectory.
    """
    adi_i = adi(angle_vec)
    pad_i = pad(adi_i)
    adicum_i = adicum(angle_vec)
    padcum_i = padcum(adicum_i)
    atip_i = atip(angle_vec)
    erre = numpy.sqrt(adi_i ** 2 + atip_i ** 2)
    tai_i = rtod * numpy.arccos((1 - erre) / (1 + erre))
    pai_i = rtod * numpy.arctan2(atip_i, adi_i)
    atin_i = atin(angle_vec)
    tag_v = tag(padcum_i, pai_i)

    return (pad_i, padcum_i, tai_i, pai_i, atin_i, tag_v)
