# -*- coding: utf-8 -*-
# Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option, any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from .Utilities import vect2str, BezierToTrajectoryString
import string
import numpy as np
from pylab import double, array, random
from .TOPPbindings import TOPPInstance

################# Reading from string #####################

def ProfileFromLines(lines):
    l = lines[0]
    [duration, dt] = [double(x) for x in l.split(' ')]
    if duration<= 0 :
        return None
    l = lines[1]
    sarray = array([double(x) for x in l.split(' ')])
    l = lines[2]
    sdarray = array([double(x) for x in l.split(' ')])
    return [duration, dt, sarray, sdarray]


def ProfilesFromString(s):
    s = s.strip(" \n")
    profileslist = []
    lines = [l.strip(" \n") for l in s.split('\n')]
    n = int(len(lines) / 3)
    for i in range(n):
        l = ProfileFromLines(lines[3 * i:3 * i + 3])
        if l is None:
            continue
        profileslist.append(ProfileFromLines(lines[3 * i:3 * i + 3]))
    return profileslist


def ExtraFromString(s):
    s = s.strip(" \n")
    lines = [l.strip(" \n") for l in s.split('\n')]
    lines.pop(0)
    tvect = []
    torques = []
    for i in range(len(lines) / 2):
        tvect.append(double(lines[2 * i]))
        torques.append(array([double(x) for x in lines[2 * i + 1].split(' ')]))
    return array(tvect), array(torques)


def SwitchPointsFromString(s):
    if len(s) == 0:
        return []
    s = s.strip(" \n")
    switchpointslist = []
    lines = [l.strip(" \n") for l in s.split('\n')]
    for l in lines:
        switchpointslist.append(VectorFromString(l))
    return switchpointslist


def VectorFromString(s):
    # left for compatibility TODO: remove?
    s = s.strip(" \n")
    return array([double(x) for x in s.split(' ')])


def GenerateRandomTrajectory(ncurve, ndof, bound):
    def vector2string(v):
        ndof = len(v)
        s = str(ndof)
        for a in v:
            s += ' %f' % a
        return s

    p0a = vector2string(random(ndof) * 2 * bound - bound)
    p0b = vector2string(random(ndof) * 2 * bound - bound)
    p1a = vector2string(random(ndof) * 2 * bound - bound)
    p1b = vector2string(random(ndof) * 2 * bound - bound)
    s = '%d' % ncurve
    s += '\n1.0 ' + p0a + ' ' + p0b
    for k in range(ncurve - 1):
        a = random(ndof) * 2 * bound - bound
        b = random(ndof) * 2 * bound - bound
        c = 2 * b - a
        pa = vector2string(a)
        pb = vector2string(b)
        pc = vector2string(c)
        s += ' ' + pa + ' ' + pb + '\n1.0 ' + pb + ' ' + pc
    s += ' ' + p1a + ' ' + p1b
    Tv, p0v, p1v, p2v, p3v = string2p(s)
    return BezierToTrajectoryString(Tv, p0v, p1v, p2v, p3v)


################# Compute Kinematic Constraints #####################

def ComputeKinematicConstraints(traj, amax, discrtimestep):
    # Sample the dynamics constraints
    ndiscrsteps = int((traj.duration + 1e-10) / discrtimestep) + 1
    constraintstring = ""
    for i in range(ndiscrsteps):
        t = i * discrtimestep
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        constraintstring += "\n" + vect2str(+qd) + " " + vect2str(-qd)
        constraintstring += "\n" + vect2str(+qdd) + " " + vect2str(-qdd)
        constraintstring += "\n" + vect2str(-amax) + " " + vect2str(-amax)
    return constraintstring


################# Compute MMR Constraints #####################

def ComputeMaterialRemovalConstraints(traj, amax, mrr_desired, volume_rate,
                                      discrtimestep):
    # Sample the MMR constraints. Some linear interpolation for volume.
    ndiscrsteps = int((traj.duration + 1e-10) / discrtimestep) + 1
    constraintstring = ""
    for i in range(ndiscrsteps):
        t = i * discrtimestep
        qd = traj.Evald(t)
        qdd = traj.Evaldd(t)
        volume = volume_rate.Eval(t)
        speed2 = np.linalg.norm(qd)**2
        constraintstring += "\n" + vect2str(+qd) + " " + vect2str(-qd) + " " + str(0)
        # MRR = Vol / Time = vol / (dist/speed) = Vol * speed / dist
        constraintstring += "\n" + vect2str(+qdd) + " " + vect2str(-qdd) + " " + str(volume**2)
        constraintstring += "\n" + vect2str(-amax) + " " + vect2str(-amax) + " " + str(-(mrr_desired**2))
    return constraintstring


######################## Plots ############################

def PlotProfiles(profileslist0, switchpointslist=[], max_sdot=20, figstart=None, colorscheme = 1):
    from pylab import figure, clf, hold, plot, gca, axis, title, xlabel, ylabel, cycler
    profileslist = list(profileslist0)
    if figstart is not None:
        figure(figstart)
        clf()
    mvcbobrow = profileslist.pop(0)
    mvcdirect = profileslist.pop(0)
    if colorscheme == 1:
        plot(mvcbobrow[2], mvcbobrow[3], 'm', linewidth=1)
        plot(mvcdirect[2], mvcdirect[3], 'g--', linewidth=1)
    else:
        plot(mvcbobrow[2], mvcbobrow[3], 'm', linewidth=1)
        plot(mvcdirect[2], mvcdirect[3], 'g--', linewidth=1)
    colorcycle = cycler('color', ['r', 'g', 'b', 'c', 'm'])
    ax = gca()
    ax.set_prop_cycle(colorcycle)
    for p in profileslist:
        plot(p[2], p[3], 'k',linewidth=2)
    if len(profileslist) > 0:
        M = 2 * max([max(p[3]) for p in profileslist])
    else:
        M = max_sdot
        bobrow = list(filter((lambda x: x < M), mvcbobrow[3]))
        direct = list(filter((lambda x: x < M), mvcdirect[3]))
        if len(bobrow) > 0:
            M = max(M, max(bobrow))
        if len(direct) > 0:
            M = max(M, max(direct))
    for sw in switchpointslist:
        if sw[2] == 0:
            plot(sw[0], sw[1], 'ro', markersize=8)
        if sw[2] == 1:
            plot(sw[0], sw[1], 'go', markersize=8)
        if sw[2] == 2:
            plot(sw[0], sw[1], 'bo', markersize=8)
        if sw[2] == 3:
            plot(sw[0], sw[1], 'yo', markersize=8)
    smax, sdmax = mvcbobrow[0], M
    axis([0, smax, 0, sdmax])
    if colorscheme == 1:
        title('Maximum Velocity Curves and profiles', fontsize=20)
        xlabel('$s$', fontsize=22)
        ylabel('$\dot s$', fontsize=22)
    return smax, sdmax  # return this for PlotPhase (yurk!)


def PlotTSMap(traj, svalues, figstart=1):
    from pylab import figure, clf, hold, plot, gca, axis, title, xlabel, ylabel, cycler
    figure(figstart)
    clf()
    times = []
    current_time = 0
    for chunk in traj.chunkslist:
        times.append(current_time)
        current_time = current_time + chunk.duration
    plot(svalues, times)
    #title("Variation of t with s", fontsize=20)
    hfont = {'fontname':'Arial'}
    xlabel("$s$", fontsize=18, **hfont)
    ylabel("$t$", fontsize=18, **hfont)

def PlotComputedProfiles(topp_bind, figstart=1):
    topp_bind.WriteProfilesList()
    topp_bind.WriteSwitchPointsList()
    profileslist = ProfilesFromString(topp_bind.resprofilesliststring)
    switchpointslist = SwitchPointsFromString(topp_bind.switchpointsliststring)
    PlotProfiles(profileslist, switchpointslist, figstart)


def PlotAlphaBeta(topp_inst, prec=30):
    from pylab import axis, linspace, sqrt, plot
    smin, smax, sdmin, sdmax = axis()
    if sdmin <= 0.:
        sdmin = 1e-2
    s_coord = linspace(smin, smax, prec)
    sd_coord = linspace(sdmin, sdmax, prec)
    ds0 = s_coord[1] - s_coord[0]
    dsd0 = sd_coord[1] - sd_coord[0]
    alpha = lambda s, sd: topp_inst.GetAlpha(s, sd) / sd
    beta = lambda s, sd: topp_inst.GetBeta(s, sd) / sd
    yscl = dsd0 / ds0
    for s in s_coord:
        for sd in sd_coord:
            ds = ds0 / 2
            a, b = alpha(s, sd), beta(s, sd)
            na, nb = 1. / sqrt(1. + a ** 2), 1. / sqrt(1. + b ** 2)
            na = 1 / sqrt(1 + (a / yscl) ** 2)
            nb = 1 / sqrt(1 + (b / yscl) ** 2)
            plot([s, s + na * ds], [sd, sd + na * a * ds], 'b', alpha=.3)
            plot([s, s + nb * ds], [sd, sd + nb * b * ds], 'r', alpha=.3)
            if a > b:
                plot([s, s], [sd, sd], 'ko', alpha=.3, markersize=3)
    axis([smin, smax, sdmin, sdmax])


def PlotKinematics(traj0, traj1, dt=0.01, vmax=[], amax=[], figstart=0):
    from pylab import figure, clf, hold, gca, title, xlabel, ylabel, plot, axis, cycler, show
    x = ['r', 'g', 'b', 'c', 'm']
    colorcycle = cycler('color', x[0:traj0.dimension])
    Tmax = max(traj0.duration, traj1.duration)

    # Joint angles
    figure(figstart)
    clf()
    ax = gca()
    ax.set_prop_cycle(colorcycle)
    traj0.Plot(dt, '--')
    ax.set_prop_cycle(colorcycle)
    traj1.Plot(dt)
    title('Joint values', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('Joint values', fontsize=18)

    # Velocity
    figure(figstart + 1)
    clf()
    ax = gca()
    ax.set_prop_cycle(colorcycle)
    traj0.Plotd(dt, '--')
    ax.set_prop_cycle(colorcycle)
    traj1.Plotd(dt)
    for v in vmax:
        plot([0, Tmax], [v, v], '-.')
    for v in vmax:
        plot([0, Tmax], [-v, -v], '-.')
    if len(vmax) > 0:
        Vmax = 1.2 * max(vmax)
        if Vmax < 0.1:
            Vmax = 10
        axis([0, Tmax, -Vmax, Vmax])
    title('Joint velocities', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('Joint velocities', fontsize=18)

    # Acceleration
    figure(figstart + 2)
    clf()
    ax = gca()
    ax.set_prop_cycle(colorcycle)
    traj0.Plotdd(dt, '--')
    ax.set_prop_cycle(colorcycle)
    traj1.Plotdd(dt)
    for a in amax:
        plot([0, Tmax], [a, a], '-.')
    for a in amax:
        plot([0, Tmax], [-a, -a], '-.')
    if len(amax) > 0:
        Amax = 1.2 * max(amax)
        axis([0, Tmax, -Amax, Amax])
    title('Joint accelerations', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('Joint accelerations', fontsize=18)

def PlotMRR(traj, volume_rates, svalues, dt=0.01, mrr_desired=[], figstart=0):
    from pylab import figure, clf, hold, gca, title, xlabel, ylabel, plot, axis, cycler
    x = ['r', 'g', 'b', 'c', 'm']
    colorcycle = cycler('color', x[0:traj.dimension])
    Tmax = traj.duration

    # Material removal rate
    figure(figstart)
    clf()
    ax = gca()
    ax.set_prop_cycle(colorcycle)
    tvect = np.arange(0, traj.duration + dt, dt)

    # Compute times for tsmap
    times = []
    current_time = 0
    for chunk in traj.chunkslist:
        times.append(current_time)
        current_time = current_time + chunk.duration

    # Compute for each time, if it is contained in interval [tcur, tcur+dt] for
    # tvect[0:-1]
    time_contained_in_dt = np.asarray([np.logical_and(times >= tvect[index],
                                                      times < tvect[index+1])
                                       for index in np.arange(tvect.size-1)])
    # For each interval of tvect[0:-1], compute s values that are traversed in
    # interval. Correspond these s-values to the beginning of the interval
    s_for_dts = [svalues[time_contained_in_dt[index, :]] for index in np.arange(time_contained_in_dt.shape[0])]
    # Use the s values to compute the ds values for each interval
    ds_for_dts_firsts = [0] + [s_for_dt[0]-s_for_dt_last[-1] for s_for_dt, s_for_dt_last in
                         zip(s_for_dts[1:], s_for_dts[:-1])]
    ds_for_dts_others = [tuple(s_for_dt[1:]-s_for_dt[:-1]) for s_for_dt in s_for_dts]
    ds_for_dts = [np.asarray((ds_for_dt_first,) + ds_for_dt_others) for
                  ds_for_dt_first, ds_for_dt_others in zip(ds_for_dts_firsts,
                                                           ds_for_dts_others)]

    volumes_dt = np.zeros_like(tvect)
    volumes_dt[1:] = [np.sum(np.asarray([volume_rates.Eval(s_for_dt[index]) *
                           ds_for_dt[index] for index in
                           np.arange(s_for_dt.size)])) for s_for_dt,
                           ds_for_dt in zip(s_for_dts, ds_for_dts)]

    qdvect = array([traj.Evald(t) for t in tvect])
    #plot(tvect, np.linalg.norm(qdvect, axis=1), 'r--', linewidth=2)
    #plot(tvect, volumes_dt, 'b--', linewidth=2)
    plot(tvect, (volumes_dt/dt), 'g', linewidth=2)
    for mrr in mrr_desired:
        plot([0, Tmax], [mrr, mrr], 'b-.')
    for mrr in mrr_desired:
        plot([0, Tmax], [0, 0], 'b-.')
    #if len(mrr_desired) > 0:
    #    Vmax = 1.2 * max(vmax)
    #    if Vmax < 0.1:
    #        Vmax = 10
    #    axis([0, Tmax, -Vmax, Vmax])
    title('Material removal rate', fontsize=20)
    xlabel('Time (s)', fontsize=18)
    ylabel('MRR', fontsize=18)

def string2p(s):
    lines = [l.strip(" \n") for l in s.split('\n')]
    Tv = []
    p0v = []
    p1v = []
    p2v = []
    p3v = []
    for i in range(1, len(lines)):
        l = [float(x) for x in lines[i].split(' ')]
        Tv.append(l.pop(0))
        ndof = int(l[0])
        p0v.append(l[1:ndof + 1])
        p1v.append(l[ndof + 2:2 * (ndof + 1)])
        p2v.append(l[2 * (ndof + 1) + 1:3 * (ndof + 1)])
        p3v.append(l[3 * (ndof + 1) + 1:4 * (ndof + 1)])
    return Tv, p0v, p1v, p2v, p3v

