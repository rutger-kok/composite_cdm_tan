'''
This script runs single element verification simulations to check a
subroutine correctly predicts failure stresses and dissipates the correct
amount of energy.

To run this script:
- in Abaqus CAE first set the working directory to the directory containing
  the Python script and the required input files, then go to File -> Run script
OR
- in a terminal window navigate to the directory containing this script and the
  required input files, then enter: abaqus cae noGUI=verify_subroutine.py

Last updated: 19/04/2021

(c) Rutger Kok 2021
'''

from abaqus import *
from abaqusConstants import *
from caeModules import *
from visualization import *
import time
import csv
import numpy as np
import shutil
import os


def run_job(input_file):
    '''Run Abaqus job given input file'''

    mdb.ModelFromInputFile(name=input_file, inputFileName=input_file + '.inp')
    subroutine_path = "C:\\GitHub\\composite_cdm\\composite_cdm.for"
    mdb.Job(
        name=input_file, model=input_file, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF,
        userSubroutine=subroutine_path, numCpus=1,
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
    ss_indices = {'11': 0, '22': 1, '33': 2, '12': 3, '23': 4, '13': 5}
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
    # indices of fracture energies
    energy_table = {'tension_11': [18, 19], 'tension_22': 21,
                    'compression_11': 20, 'compression_22': None,
                    'shear_12': 22, 'shear_23': 22}
    mat_props = odb.materials['VTC401'].userMaterial.mechanicalConstants
    if isinstance(energy_table[input_file], list):
        fracture_energy = 0.0
        for energy_index in energy_table[input_file]:
            fracture_energy += mat_props[energy_index]
    elif energy_table[input_file] is None:
        fracture_energy = mat_props[22] / cos(radians(mat_props[16]))
    else:
        fracture_energy = mat_props[energy_table[input_file]]

    # integrate stress-strain data to determine dissipated energy
    dissipated_energy = np.trapz(stress_data, strain_data)

    # check failure stresses
    strength_table = {'tension_11': 9, 'tension_22': 12, 'compression_11': 11,
                      'compression_22': 13, 'shear_12': 14, 'shear_23': 15}
    strength = mat_props[strength_table[input_file]]
    max_stress = max([abs(x) for x in stress_data])
    data = [input_file, fracture_energy, dissipated_energy, strength,
            max_stress]
    return data


def plot_data(input_file, data):
    plot_name = input_file
    x_axis_title = 'Logarithmic Strain (-) '
    y_axis_title = 'Stress (GPa)'
    xy_data = session.XYData(plot_name, data)
    curve = session.Curve(xy_data)
    xy_plot = session.XYPlot(plot_name)
    chart = xy_plot.charts.values()[0]
    chart.setValues(curvesToPlot=(curve, ))
    chart.axes1[0].axisData.setValues(useSystemTitle=False,
                                      title=x_axis_title)
    chart.axes2[0].axisData.setValues(useSystemTitle=False,
                                      title=y_axis_title)

    # Display the XY Plot in the current viewport
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=xy_plot)


def main():
    # create verification directory
    cwd = os.getcwd()
    wd = 'C:\\Workspace\\verification_' + str(int(time.time()))
    if not os.path.exists(wd):
        os.makedirs(wd)

    # copy input files to verification directory
    for f in os.listdir(cwd):
        if f.endswith(".inp"):
            shutil.copy(os.path.join(cwd, f), wd)
    os.chdir(wd)

    all_data = []
    for job in ('tension_11', 'tension_22',
                'compression_11', 'compression_22',
                'shear_12', 'shear_23'):
        run_job(job)
        test_data = post_process(job)
        all_data.append(test_data)

    csv_headings = ['Name', 'Fracture Energy', 'Dissipated Energy',
                    'Strength', 'Max Stress', 'Test Status']

    with open('verification_log.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(csv_headings)
        writer.writerows(all_data)


if __name__ == '__main__':
    main()
