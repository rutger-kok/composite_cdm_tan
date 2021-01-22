'''
This script runs each single element verifcation simulation to check the
subroutine correctly predicts failure stresses and dissipates the correct
amount of energy.

To run this script:
- in Abaqus CAE go to File -> Run script
- in a terminal window navigate to the directory containing this scipt and the
  required input files, then enter: abaqus cae noGUI=verify_subroutine.py

Last updated: 22/01/2021

(c) Rutger Kok 2021
'''

from abaqus import *
from abaqusConstants import *
from caeModules import *
from visualization import *
import numpy as np
import subprocess


def run_job(input_file):
    '''Use subproces module to run Abaqus job given input file'''
    subprocess.call(["abaqus", "job={}".format(input_file),
                     "user=composite_cdm.for", "double=both"])


def post_process(input_file):
    '''
    Open output database for given input file, extract and plot results,
    then export data to .csv
    '''
    # Open the output database
    odb = openOdb(path='{}.odb'.format(input_file))

    # extract stress-strain data from single element
    stress_data = []  # initialize lists
    strain_data = []
    for frame in odb.steps['Step-1'].frames:
        stress = frame.fieldOutputs['S']
        s_11 = stress.values[0].data[0]
        strain = frame.fieldOutputs['LE']
        le_11 = strain.values[0].data[0]
        stress_data.append(s_11)
        strain_data.append(le_11)

    # call plot_data function to plot stress-strain curves
    data = zip(strain_data, stress_data)
    plot_data(input_file, data)

    # check the dissipated energy is the same as the fracture energy
    # first obtain fracture energy from material properties
    energy_table = {'tension_1': 15, 'tension_2': 17, 'compression_1': 16,
                    'compression_2': 18} # indexes of fracture energies
    mat_props = odb.materials['VTC401'].userMaterial.mechanicalConstants
    fracture_energy = mat_props[energy_table[input_file]]
    # integrate stress-strain data to determine dissipated energy
    energy_tol = 1E-3  # tolerance for comparison
    dissipated_energy = np.trapz(stress_data, strain_data)
    if (fracture_energy - dissipated_energy) < energy_tol:
        print '{} energy check = PASS'.format(input_file)
    else:
        print '{} energy check = FAIL'.format(input_file)
        print 'Fracture Energy = {}'.format(fracture_energy)
        print 'Dissipated Energy = {}'.format(dissipated_energy)

    # check failure stresses
    strength_table = {'tension_1': 9, 'tension_2': 11, 'compression_1': 10,
                      'compression_2': 12}
    strength = mat_props[strength_table[input_file]]
    max_stress = max(stress_data)
    strength_tol = 1E-2
    if (strength - max_stress) < strength_tol:
        print '{} strength check = PASS'.format(input_file)
    else:
        print '{} strength check = FAIL'.format(input_file)
        print 'Strength = {}'.format(strength)
        print 'Max Stress = {}'.format(max_stress)


def plot_data(input_file, data):
    plot_name = 'S11 vs LE11'
    x_axis_title = 'LE11'
    y_axis_title = 'S11'
    xy_data = session.XYData(plot_name, data)
    curve = session.Curve(xy_data)
    xy_plot = session.XYPlot(plot_name)
    chart = xy_plot.charts.values()[0]
    chart.setValues(curvesToPlot=(plot_name, ))
    chart.axes1[0].axisData.setValues(useSystemTitle=False,
                                      title=x_axis_title)
    chart.axes2[0].axisData.setValues(useSystemTitle=False,
                                      title=y_axis_title)

    # Display the XY Plot in the current viewport
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=xy_plot)


def main():
    for job in ('tension_1', 'compression_1', 'tension_2'):
        run_job(job)
        post_process(job)


if __name__ == '__main__':
    main()
