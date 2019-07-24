from odbAccess import *
from visualization import XYData, USER_DEFINED

# filePath = r"C:\Workspace\check\cube_compression1.odb"
filePath = r"C:\Workspace\check\cube_compression1.odb"
# filePath = r"C:\Workspace\catalanotti_failure_criteria\catalanotti_tension2\cube_tension2.odb"
# filePath = r"C:\Workspace\catalanotti_failure_criteria\catalanotti_compression2\cube_compression2.odb"
# filePath = r"C:\Workspace\check\cube_tension1.odb"
odb = openOdb(filePath)
name = filePath[filePath.rfind('\\')+1:filePath.rfind('.')]

data = []
for frame in odb.steps['Step-1'].frames:
    stress = frame.fieldOutputs['S']
    s11 = stress.values[0].data[0]
    strain = frame.fieldOutputs['LE']
    le11 = strain.values[0].data[0]
    data.append((le11,s11))

plotName = 'S11 vs LE11'
xAxisTitle = 'LE11'
yAxisTitle = 'S11'
xyData = session.XYData(plotName, data)
curve = session.Curve(xyData)
xyPlot = session.XYPlot(plotName)

data2 = []
for frame2 in odb.steps['Step-1'].frames:
    sdv13 = frame.fieldOutputs['SDV13']
    d1Plus = sdv13.values[0].data
    strain = frame.fieldOutputs['LE']
    le11 = strain.values[0].data[0]
    data2.append((le11,d1Plus))

plotName2 = 'd1Plus vs LE11'
xyData2 = session.XYData(plotName2, data2)
curve2 = session.Curve(xyData2)
xyPlot2 = session.XYPlot(plotName2)

chart = xyPlot.charts.values()[0]
chart.setValues(curvesToPlot=(plotName,plotName2, ))
chart.axes1[0].axisData.setValues(useSystemTitle=False,
                                    title=xAxisTitle)
chart.axes2[0].axisData.setValues(useSystemTitle=False,
                                    title=yAxisTitle)

# Display the XY Plot in the current viewport
vp = session.viewports[session.currentViewportName]
vp.setValues(displayedObject=xyPlot)
