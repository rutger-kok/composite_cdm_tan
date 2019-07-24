from odbAccess import *
import os 
import csv 

# Search for .odb files in C:\Workspace
odb_file_path_list = []
root_path_list = []
for root, dirs, files in os.walk(r'path to directory containing odbs here'):
    for file in files:
        if file.endswith('.odb'):
            root_path_list.append(root)
            odb_file_path_list.append(os.path.join(root, file))

data = []
for path in odb_file_path_list:
    odb = openOdb(path)

    # whatever you want to do with the odb data

    # I left my code here as an example
    name = path[path.rfind('\\')+1:path.rfind('.')]

    stress = odb.steps['Loading Step'].frames[-1].fieldOutputs['S']
    undulationInstance = odb.rootAssembly.instances['Undulation Instance']
    undulationCells = undulationInstance.elementSets['All Undulation Cells']
    subset_stress = stress.getSubset(region=undulationCells)
    s11 = [subset_stress.values[x].data[0] for x in range(
            len(subset_stress.values))]
    s11_avg = (sum(s11)/len(s11))*10**3
    s11_max = max(s11)*10**3
    s11_min = min(s11)*10**3
    sif = s11_max/s11_avg
    relrange = (s11_max-s11_min)/s11_avg

    u =  odb.steps['Loading Step'].frames[-1].fieldOutputs['U']
    left_side = odb.rootAssembly.nodeSets['Left Side']
    left_side_deflection = u.getSubset(region=left_side)
    left_side_u1 = [left_side_deflection.values[z].data[0] for z 
                        in range(len(left_side_deflection.values))]
    left_side_avg_u1 = sum(left_side_u1)/len(left_side_u1)

    right_side = odb.rootAssembly.nodeSets['Right Side']
    right_side_deflection = u.getSubset(region=right_side)
    right_side_u1 = [right_side_deflection.values[zz].data[0] for zz 
                        in range(len(right_side_deflection.values))]
    right_side_avg_u1 = sum(right_side_u1)/len(right_side_u1)

    total_u1 = right_side_avg_u1 - left_side_avg_u1

    rf =  odb.steps['Loading Step'].frames[-1].fieldOutputs['RF']
    refpoint2 = odb.rootAssembly.nodeSets['Coupling Reference Point 2']
    subset_force = rf.getSubset(region=refpoint2)
    rf_1 = subset_force.values[0].data[0]

    L_u = (right_side.nodes[0][0].coordinates[0]
            -left_side.nodes[0][0].coordinates[0])
    data.append(
        [name, s11_avg, s11_max, s11_min, sif, rf_1, total_u1, L_u])

# create headings for csv file
headings = ['Name', 'Avg S11', 'Max S11', 'Min S11', 'SIF', 'RF1', 'Displacement', 'Lu']

# write data and heading to csv
with open('undulation_ps_results.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow(headings)
    writer.writerows(data)
