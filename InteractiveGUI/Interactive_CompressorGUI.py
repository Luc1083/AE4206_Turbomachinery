import PySimpleGUI as sg
from Turbo_Emulators import CompressorEstimator
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import matplotlib
matplotlib.use('TkAgg')

# Interactive interface for students to visualize the effects of compressor geometry and operation on the flow regime.
# Authors: ir. E.C.Bunschoten
# Delft University of Technology - All rights reserved


sw_default = 0.5 

# Define compressor estimation object.
C = CompressorEstimator()
C.ComputePerformance()

# Prepare all the sliders
options_column = [
    [sg.Text("RPM"),sg.Slider(range=(0.01*C.GetRPM(), 2*C.GetRPM()), default_value=C.GetRPM(), orientation='horizontal',key='-CHANGE_RPM-',enable_events=True)],
    [sg.Text("Total temperature[K]"),sg.Slider(range=(0.01*C.GetTt(), 2*C.GetTt()), default_value=C.GetTt(), orientation='horizontal',key='-CHANGE_Tt-',enable_events=True)],
    [sg.Text("Total pressure[Pa]"),sg.Slider(range=(0.01*C.GetPt(), 2*C.GetPt()), default_value=C.GetPt(), orientation='horizontal',key='-CHANGE_Pt-',enable_events=True)],
    [sg.Text("Mass flow rate[kg/s]"), sg.Slider(range=(0.01*C.GetMdot(), 2*C.GetMdot()), default_value=C.GetMdot(), orientation='horizontal',resolution=1,enable_events=True,key='-CHANGE_mdot-')],
    [sg.Text("t/c ratio"),sg.Slider(range=(0,3*C.GetThicknessRatio()), resolution=0.001, default_value=C.GetThicknessRatio(), orientation='horizontal',key='-CHANGE_tc-',enable_events=True)],
    [sg.Text("Beta root[deg]"),sg.Slider(range=(-89, 0), default_value=C.GetBetaBlade(0), orientation='horizontal',key='-CHANGE_beta_0-',enable_events=True)],
    [sg.Text("Beta mid[deg]"),sg.Slider(range=(-89, 0), default_value=C.GetBetaBlade(1), orientation='horizontal',key='-CHANGE_beta_1-',enable_events=True)],
    [sg.Text("Beta tip[deg]"),sg.Slider(range=(-89, 0), default_value=C.GetBetaBlade(2), orientation='horizontal',key='-CHANGE_beta_2-',enable_events=True)],
    [sg.Text("Hub radius[m]"),sg.Slider(range=(0.05, 0.30), resolution=0.01,default_value=C.GetHubRadius(),orientation='horizontal',key='-CHANGE_R_HUB-',enable_events=True)],
    [sg.Text("Blade height[m]"),sg.Slider(range=(0.05, 0.3),resolution=0.01,default_value=C.GetShroudRadius()-C.GetHubRadius(),orientation='horizontal',key='-CHANGE_BH-',enable_events=True)],
    [sg.Text("Blade count"),sg.Slider(range=(5, 2*C.GetBladeCount()),resolution=1,default_value=C.GetBladeCount(),orientation='horizontal',key='-CHANGE_BC-',enable_events=True)],
    [sg.Text("Span-wise location"),sg.Slider(range=(0, 1), default_value=0.5, orientation='vertical', enable_events=True,resolution=0.01,key='-SWSL-')]

]

# Column with visualization
vis_column = [[sg.Text("Blade-to-blade view")],
    [sg.Graph((640, 480), (0, 0), (640, 480),key='Window')]
]

# Prepare window layout
layout = [
    [
        sg.Column(vis_column),
        sg.VSeperator(),
        sg.Column(options_column),
    ]
]

# Generate GUI window
window = sg.Window("Compressor Performance", layout,finalize=True,font=("monospace",12))

# Functions to link figure with window.
def draw_figure(canvas, figure):
   figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
   figure_canvas_agg.draw()
   figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
   return figure_canvas_agg

def pack_figure(graph, figure):
    canvas = FigureCanvasTkAgg(figure, graph.Widget)
    plot_widget = canvas.get_tk_widget()
    plot_widget.pack(side='top', fill='both', expand=1)
    return plot_widget

# Prepare matploblib subplot for visualizations.
fig, axs = plt.subplots(1,2)
graph1 = window['Window']
pack_figure(graph1, fig)

# Plot flow and blades mid-span.
current_sw = sw_default
C.PlotBladeAndFlow(current_sw, axs[0])
C.PlotSpanWiseResults(axs[1], current_sw)

# Function for updating the visualization window.
def UpdatePlots(spanwise_val):
    axs[0].cla()
    axs[1].cla()
    C.PlotBladeAndFlow(spanwise_val, axs[0])
    C.PlotSpanWiseResults(axs[1], spanwise_val)
    fig.canvas.draw()

while True:

    # Read inputs.
    event, values = window.read() 
    if event == '-SWSL-':
        spanwise_val = float(values['-SWSL-'])
        current_sw = spanwise_val
        UpdatePlots(current_sw)

    if event == '-CHANGE_RPM-':
        val_RPM = float(values['-CHANGE_RPM-'])
        C.SetRPM(val_RPM)
        C.ComputePerformance()
        UpdatePlots(current_sw)

    if event == '-CHANGE_Tt-':
        val_Tt = float(values['-CHANGE_Tt-'])
        C.SetTt(val_Tt)
        C.ComputePerformance()
        UpdatePlots(current_sw)  

    if event == '-CHANGE_Pt-':
        val_Pt = float(values['-CHANGE_Pt-'])
        C.SetPt(val_Pt)
        C.ComputePerformance()
        UpdatePlots(current_sw)  
    
    if event == '-CHANGE_mdot-':
        val_mdot = float(values['-CHANGE_mdot-'])
        C.SetMdot(val_mdot)
        C.ComputePerformance()
        UpdatePlots(current_sw)  

    if event == '-CHANGE_tc-':
        val_tc = float(values['-CHANGE_tc-'])
        C.SetThicknessRatio(val_tc)
        C.ComputePerformance()
        UpdatePlots(current_sw)  

    if event == '-CHANGE_beta_0-':
        val_beta_0 = float(values['-CHANGE_beta_0-'])
        C.SetBetaBlade(val_beta_0, 0)
        C.ComputePerformance()
        UpdatePlots(current_sw)  

    if event == '-CHANGE_beta_1-':
        val_beta_1 = float(values['-CHANGE_beta_1-'])
        C.SetBetaBlade(val_beta_1, 1)
        C.ComputePerformance()
        UpdatePlots(current_sw) 

    if event == '-CHANGE_beta_2-':
        val_beta_2 = float(values['-CHANGE_beta_2-'])
        C.SetBetaBlade(val_beta_2, 2)
        C.ComputePerformance()
        UpdatePlots(current_sw) 
    
    if event == '-CHANGE_R_HUB-':
        val_R_hub = float(values['-CHANGE_R_HUB-'])
        C.SetHubRadius(val_R_hub)
        C.ComputePerformance()
        UpdatePlots(current_sw)
    
    if event == '-CHANGE_BH-':
        val_BH = float(values['-CHANGE_BH-'])
        R_tip = C.GetHubRadius()+val_BH
        C.SetShroudRadius(R_tip)
        C.ComputePerformance()
        UpdatePlots(current_sw)

    if event == '-CHANGE_BC-':
        val_BC = int(values['-CHANGE_BC-'])
        C.SetBladeCount(val_BC)
        C.ComputePerformance()
        UpdatePlots(current_sw)
        
    if event == "Exit" or event == sg.WIN_CLOSED:

        break