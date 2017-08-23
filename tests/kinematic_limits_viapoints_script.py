# -*- coding: utf-8 -*-
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
path = data[1:4,:]
volumes = data[4,:]
traj0 = Utilities.InterpolateViapoints(path) # Interpolate using splines

# Constraints
vmax = float(sys.argv[2])*ones(traj0.dimension)
amax = float(sys.argv[3])*ones(traj0.dimension)
mrr_desired = float(sys.argv[4])

# Set up the TOPP instance
trajectorystring = str(traj0)
discrtimestep = float(sys.argv[5])
uselegacy = False
if uselegacy: #Using the legacy KinematicLimits (a bit faster but not fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    constraintstring += "\n" + string.join([str(a) for a in amax])
    x = TOPPbindings.TOPPInstance(None,"KinematicLimits",constraintstring,trajectorystring);
else: #Using the general QuadraticConstraints (fully supported)
    constraintstring = str(discrtimestep)
    constraintstring += "\n" + string.join([str(v) for v in vmax])
    #constraintstring += TOPPpy.ComputeKinematicConstraints(traj0, amax, discrtimestep)
    constraintstring += TOPPpy.ComputeMaterialRemovalConstraints(traj0,
                                                        amax, mrr_desired, volumes, discrtimestep)
    x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring);

# Run TOPP
ret = x.RunComputeProfiles(0,0)
x.ReparameterizeTrajectory()

# Display results
ion()
x.WriteProfilesList()
x.WriteSwitchPointsList()
profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
#TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
x.WriteResultTrajectory()
traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
np.savetxt(sys.argv[6], np.asarray([traj0.duration, traj1.duration]))
dtplot = 0.01
TOPPpy.PlotKinematics(traj1,traj1,dtplot,vmax,amax)
TOPPpy.PlotMRR(traj1,volumes,dtplot,[mrr_desired])

raw_input()
