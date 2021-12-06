# coding: utf8
import sys
import time

import crocoddyl
import numpy as np

from quadruped_walkgen import ActionModelQuadruped

N = 16  # number of nodes
T = int(sys.argv[1]) if (len(sys.argv) > 1) else int(5000)  # number of trials
MAXITER = (int(sys.argv[2]) if (len(sys.argv) > 1) else int(1))  # number of maximum iteration

# Gait Matrix
gait = np.zeros((6, 5))
fsteps = np.full((6, 13), np.nan)

# MpcInterface object that contains information about the current state of the robot
# Change the initial conditions here
lC = np.array([[0.0, 0.0, 0.2]]).transpose()  # CoM centered and at 20 cm above the ground
abg = np.array([[0.0, 0.0, 0.0]]).transpose()  # horizontal base (roll, pitch, 0.0)
lV = np.array([[0.2, 0.0, 0.0]]).transpose()  # motionless base (linear velocity)
lW = np.array([[0.0, 0.0, 0.0]]).transpose()  # motionless base (angular velocity)
l_feet = np.array([
    [0.19, 0.19, -0.19, -0.19],
    [0.15005, -0.15005, 0.15005, -0.15005],
    [0.0, 0.0, 0.0, 0.0],
])  # position of feet in local frame

x0 = np.vstack((lC, abg, lV, lW))  # Current state vector
# The reference state, copy of the initial position,
xref = np.repeat(x0, N + 1, axis=1)  # Desired future state vectors
xref[6, :] = 0.0  # Target linear velocity to Vx = 0

fsteps = np.array([
    [
        1.00000000e00,
        1.90776486e-01,
        1.48962816e-01,
        4.22498932e-03,
        1.90060159e-01,
        -1.50265109e-01,
        0.00000000e00,
        -1.89740429e-01,
        1.50467686e-01,
        1.78713224e-06,
        -1.90692335e-01,
        -1.48946056e-01,
        4.22561856e-03,
    ],
    [
        7.00000000e00,
        1.90776486e-01,
        1.48962816e-01,
        4.22498932e-03,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        -1.90692335e-01,
        -1.48946056e-01,
        4.22561856e-03,
    ],
    [
        1.00000000e00,
        1.90776486e-01,
        1.48962816e-01,
        4.22498932e-03,
        1.90000000e-01,
        -1.50050000e-01,
        0.00000000e00,
        -1.90000000e-01,
        1.50050000e-01,
        0.00000000e00,
        -1.90692335e-01,
        -1.48946056e-01,
        4.22561856e-03,
    ],
    [
        7.00000000e00,
        np.nan,
        np.nan,
        np.nan,
        1.90000000e-01,
        -1.50050000e-01,
        0.00000000e00,
        -1.90000000e-01,
        1.50050000e-01,
        0.00000000e00,
        np.nan,
        np.nan,
        np.nan,
    ],
    [
        0.00000000e00,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    [
        0.00000000e00,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
])


def createProblem():
    # Cannot use 1 model for the whole control cycle, because each model depends on the position of the feet
    # And the inertia matrix depends on the reference state (approximation )
    runningModels = []
    x0 = np.zeros(12)
    x0[2] = 0.2
    for i in range(N):
        model = ActionModelQuadruped()

        # Add model to the list of model
        runningModels.append(model)

    # Terminal Model
    terminalModel = ActionModelQuadruped()

    terminalModel.forceWeights = np.zeros((12, 1))
    terminalModel.frictionWeights = 0

    problem = crocoddyl.ShootingProblem(x0, runningModels, terminalModel)

    xs = [x0] * (N + 1)
    us = [np.zeros(12)] * N

    return xs, us, problem


def updateProblem(fsteps, xref, x0, problem):
    # Update the model according to the position of the feet in local frame (fstep)
    # Construction of the gait matrix representing the feet in contact with the ground
    index = next((idx for idx, val in np.ndenumerate(fsteps[:, 0]) if val == 0.0), 0.0)[0]
    gait[:, 0] = fsteps[:, 0]
    gait[:index, 1:] = 1.0 - (np.isnan(fsteps[:index, 1::3]) | (fsteps[:index, 1::3] == 0.0))

    # Replace NaN values by zeroes
    fsteps[np.isnan(fsteps)] = 0.0
    j = 0
    k_cum = 0
    # L = []

    # Iterate over all phases of the gait
    # The first column of xref correspond to the current state
    while gait[j, 0] != 0:
        for i in range(k_cum, k_cum + np.int(gait[j, 0])):

            # Update model
            problem.runningModels[i].updateModel(
                np.reshape(fsteps[j, 1:], (3, 4), order="F"),
                xref[:, i + 1],
                gait[j, 1:],
            )

        k_cum += np.int(gait[j, 0])
        j += 1

    # Update model of the terminal model
    problem.terminalModel.updateModel(np.reshape(fsteps[j - 1, 1:], (3, 4), order="F"), xref[:, -1], gait[j - 1, 1:])

    # update initial state of the problem
    problem.x0 = x0
    xs = [x0] * (N + 1)
    us = [np.zeros(12)] * N

    return xs, us, problem


def runUpdateProblemBenchmark(fsteps, xref, x0, problem):
    xs = [x0] * (N + 1)
    us = [np.zeros(12)] * N
    duration = []
    for i in range(T):
        c_start = time.time()
        xs, us, problem = updateProblem(fsteps, xref, x0, problem)
        c_end = time.time()
        duration.append(1e3 * (c_end - c_start))

    avrg_duration = sum(duration) / len(duration)
    min_duration = min(duration)
    max_duration = max(duration)
    return avrg_duration, min_duration, max_duration


def runDDPSolveBenchmark(xs, us, problem):
    ddp = crocoddyl.SolverDDP(problem)
    duration = []
    for i in range(T):
        c_start = time.time()
        ddp.solve(xs, us, MAXITER)
        c_end = time.time()
        duration.append(1e3 * (c_end - c_start))

    avrg_duration = sum(duration) / len(duration)
    min_duration = min(duration)
    max_duration = max(duration)
    return avrg_duration, min_duration, max_duration


def runShootingProblemCalcBenchmark(xs, us, problem):
    duration = []
    for i in range(T):
        c_start = time.time()
        problem.calc(xs, us)
        c_end = time.time()
        duration.append(1e3 * (c_end - c_start))

    avrg_duration = sum(duration) / len(duration)
    min_duration = min(duration)
    max_duration = max(duration)
    return avrg_duration, min_duration, max_duration


def runShootingProblemCalcDiffBenchmark(xs, us, problem):
    duration = []
    for i in range(T):
        c_start = time.time()
        problem.calcDiff(xs, us)
        c_end = time.time()
        duration.append(1e3 * (c_end - c_start))

    avrg_duration = sum(duration) / len(duration)
    min_duration = min(duration)
    max_duration = max(duration)
    return avrg_duration, min_duration, max_duration


print("Python bindings:")
xs, us, problem = createProblem()
xs, us, problem = updateProblem(fsteps, xref, x0, problem)
avrg_duration, min_duration, max_duration = runUpdateProblemBenchmark(fsteps, xref, x0, problem)
print("  UpdateModel [ms]: {0} ({1}, {2})".format(avrg_duration, min_duration, max_duration))
avrg_duration, min_duration, max_duration = runDDPSolveBenchmark(xs, us, problem)
print("  DDP.solve [ms]: {0} ({1}, {2})".format(avrg_duration, min_duration, max_duration))
avrg_duration, min_duration, max_duration = runShootingProblemCalcBenchmark(xs, us, problem)
print("  ShootingProblem.calc [ms]: {0} ({1}, {2})".format(avrg_duration, min_duration, max_duration))
avrg_duration, min_duration, max_duration = runShootingProblemCalcDiffBenchmark(xs, us, problem)
print("  ShootingProblem.calcDiff [ms]: {0} ({1}, {2})".format(avrg_duration, min_duration, max_duration))
