import time
import numpy as np
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from functions_influence import quilt, wake_rollup
from terminal_output import print_output as po
import functions_graphics as graph
from functions_general import archive

# The BIG TODO list
# TODO: Compact variable names to use arrays ex. Body.x, Body.y, Body.z -> Body[0], Body[1], Body[2]
# TODO: swimmer_class.py, line 37: Create new 3D Geometries and add them here
# TODO: swimmer_subclasses.py, line 94: Add 3D Geometry creation here for different profile shapes
# TODO: swimmer_subclasses.py, line 277: Change neutral_axis to neutral_plane
# TODO: swimmer_subclasses.py, line 312: Update panel_positions to account for 3D objects
# TODO: parameter_classes.py, line 20: Update GeoParameters for the new 3D shapes
# TODO: functions_general.py, line 4: panel_vectors now need a second tangential vector and area of panel
# TODO: functions_general.py, line 14: point_vectors now need a second tangential vector
# TODO: functions_general.py, line 45: Determine if there needs to be a major change in transformation module
# TODO: functions_influence.py, line5: FIX EVERYTHING!!!!!!11!!!ONE!!!
# TODO: functions_influence.py, line 147: Break up wake_rollup into smaller modules and make 3D

def main():

    po().prog_title('1.0.0')
    start_time = time.time()

    COUNTER = P['COUNTER']
    DEL_T = P['DEL_T']
    DSTEP = P['DSTEP']
    TSTEP = P['TSTEP']
    T = [DEL_T*i for i in xrange(COUNTER)]
    RHO = P['RHO']
    RE = P['RE']

    SwiP = PC.SwimmerParameters(P['CE'], P['DELTA_CORE'], P['SW_KUTTA'])
    GeoP = PC.GeoVDVParameters(P['N_BODY'], P['S'], P['C'], P['K'], P['EPSILON'])
    MotP1 = PC.MotionParameters(0., 0., P['V0'], P['THETA_MAX'], P['F'], P['PHI'])

    S1 = Swimmer(SwiP, GeoP, MotP1, COUNTER)
    Swimmers = [S1]

    po().calc_input(MotP1.THETA_MAX/np.pi*180.,RE,MotP1.THETA_MAX/np.pi*180.,DEL_T)

    # Data points per cycle == 1/(F*DEL_T)
    for i in xrange(COUNTER):
        if i == 0:
            po().initialize_output(T[i])
        else:
            if np.fmod(i,P['VERBOSITY']) == 0:
                po().timestep_header(i,T[i])

        for Swim in Swimmers:
                Swim.Body.panel_positions(DSTEP, T[i])
                Swim.Body.surface_kinematics(DSTEP, TSTEP, DEL_T, T[i], i)
                Swim.edge_shed(DEL_T, i)
                Swim.wake_shed(DEL_T, i)
        quilt(Swimmers, RHO, DEL_T, i)
        wake_rollup(Swimmers, DEL_T, i)
        archive(S1.Body.AF.x_mid)
        archive(S1.Body.AF.z_mid)

        if np.fmod(i,P['VERBOSITY']) == 0:
            po().solution_output(0,0,0,0,0,0)
            po().solution_complete_output(i/float(COUNTER-1)*100.)

    total_time = time.time()-start_time
    print "Simulation time:", np.round(total_time, 3), "seconds"

    graph.body_wake_plot(Swimmers)
#    graph.cp_plot(S1.Body)

if __name__ == '__main__':
    main()