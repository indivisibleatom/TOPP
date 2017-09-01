
# Copyright (C) 2016 Quang-Cuong Pham <cuong.pham@normalesup.org>
#
# This file is part of the Time-Optimal Path Parameterization (TOPP) library.
# TOPP is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import string, time
from pylab import *
from numpy import *
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import Utilities

# A two-dof path going through 5 viapoints (0,1) - (1,1) - (5,1) - (3,2) - (5,4)
data = np.loadtxt(sys.argv[1]).transpose()
path = data[0:3,:]
volumes = data[-1,:]
traj0 = Utilities.InterpolateViapoints(path) # Interpolate using splines
volume_rates = Utilities.InterpolateVolumeRates(volumes, path) # Interpolate using
                                                          # linear splines

# Display results
plot_dt = 0.001
ion()
figure_index = 0
from pylab import figure, clf, hold, plot, gca, axis, title, xlabel, ylabel, cycler
figure(figure_index)
volume_rates.Plot(plot_dt)
xlabel('$s$')
ylabel('Volume($V$) and DV/Ds')
figure_index = figure_index + 1

# Constraints
vmax = float(sys.argv[2])*ones(traj0.dimension)
amax = float(sys.argv[3])*ones(traj0.dimension)
mrr_desired = float(sys.argv[4])

# Set up the TOPP instance
trajectorystring = str(traj0)
discrtimestep = float(sys.argv[5])

ndiscrsteps = int((traj0.duration + 1e-10) / discrtimestep) + 1

uselegacy = False
TOPPbindings.passswitchpointnsteps = 100
if uselegacy: #Using the legacy KinematicLimits (a bit faster but not fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    constraintstring += "\n" + string.join([str(a) for a in amax])
    x = TOPPbindings.TOPPInstance(None,"KinematicLimits",constraintstring,trajectorystring);
else: #Using the general QuadraticConstraints (fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    #constraintstring += TOPPpy.ComputeKinematicConstraints(traj0, amax, discrtimestep)
    constraintstring += TOPPpy.ComputeMaterialRemovalConstraints(traj0, amax,
                                                                 mrr_desired,
                                                                 volume_rates,
                                                                 discrtimestep)
    x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring);

# Run TOPP
ret = x.RunComputeProfiles(0,0)
x.ReparameterizeTrajectory(0.0005)

x.WriteProfilesList()
x.WriteSwitchPointsList()

x.WriteTSMap()
svalues = np.fromstring(x.tsmapstring, sep='\n')

profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
TOPPpy.PlotProfiles(profileslist,switchpointslist,figure_index)

x.WriteResultTrajectory()
traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
TOPPpy.PlotTSMap(traj1,svalues,figure_index+1)
TOPPpy.PlotKinematics(traj1,traj1,plot_dt,vmax,amax,figure_index+2)
TOPPpy.PlotMRR(traj1,volume_rates,svalues,plot_dt,[mrr_desired],5)
np.savetxt(sys.argv[6], np.asarray([traj0.duration, traj1.duration]))

raw_input()
