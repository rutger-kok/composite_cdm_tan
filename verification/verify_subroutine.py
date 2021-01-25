'''
This script runs each single element verifcation simulation to check the
subroutine correctly predicts failure stresses and dissipates the correct
amount of energy.

To run this script:
- in Abaqus CAE first set the working directory to the directory containing
  the Python script and the required input files, then go to File -> Run script
OR
- in a terminal window navigate to the directory containing this script and the
  required input files, then enter: abaqus cae noGUI=verify_subroutine.py

Last updated: 25/01/2021

(c) Rutger Kok 2021
'''

from abaqus import *
from abaqusConstants import *
from caeModules import *
from visualization import *
import numpy as np


def run_job(input_file):
    '''Run Abaqus job given input file'''
    mdb.ModelFromInputFile(name=input_file, inputFileName=input_file + '.inp')
    mdb.Job(
        name=input_file, model=input_file, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF,
        userSubroutine='..\\composite_cdm.for', numCpus=1,
        scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
        numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT)
    mdb.jobs[input_file].submit(consistencyChecking=OFF)
    mdb.jobs[input_file].waitForCompletion()


def post_process(input_file):
    '''
    Open output database for given input file, extract and plot results,
    then export data to .csv
    '''
    # Open the output database
    odb = openOdb(path='{}.odb'.format(input_file))

    # get direction of loading from input file name
    _, loading_direction = input_file.split('_')
    # table linking loading directions to their corresponding stress and strain
    # indices in an Abaqus ODB
    ss_indices = {'11': 0, '22': 1, '33': 2, '12': 3, '23': 4, '31': 5}
    ss_index = ss_indices[loading_direction]

    # extract stress-strain data (NOTE: currently single element only)
    stress_data = []  # initialize lists
    strain_data = []
    for frame in odb.steps['Step-1'].frames:
        try:
            stress = frame.fieldOutputs['S']
            s = stress.values[0].data[ss_index]
            strain = frame.fieldOutputs['LE']
            le = strain.values[0].data[ss_index]
            stress_data.append(s)
            strain_data.append(le)
        except KeyError:  # ignore frames with the stres-strain data
            continue

    # call plot_data function to plot stress-strain curves
    data = zip(strain_data, stress_data)
    plot_data(input_file, data)

    # check the dissipated energy is the same as the fracture energy
    # first obtain fracture energy from material properties
    energy_table = {'tension_11': 15, 'tension_22': 17, 'compression_11': 16,
                    'compression_22': 18}  # indices of fracture energies
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
    strength_table = {'tension_11': 9, 'tension_22': 11, 'compression_11': 10,
                      'compression_22': 12}
    strength = mat_props[strength_table[input_file]]
    max_stress = max([abs(x) for x in stress_data])
    strength_tol = 1E-2
    if (strength - max_stress) < strength_tol:
        print '{} strength check = PASS'.format(input_file)
    else:
        print '{} strength check = FAIL'.format(input_file)
        print 'Strength = {}'.format(strength)
        print 'Max Stress = {}'.format(max_stress)


def plot_data(input_file, data):
    plot_name = input_file
    x_axis_title = 'Logarithmic Strain (-) '
    y_axis_title = 'Stress (GPa)'
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
    # set working directory to the directory this script is run from
    for job in ('tension_11', 'tension_22', 'compression_11'):
        run_job(job)
        post_process(job)


if __name__ == '__main__':
    main()
