# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_daq as daq
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import math
import re
import os
pd.options.plotting.backend = "plotly"

# Create the dashboard
app = dash.Dash(__name__)



# Define here some parameters
Rcst=8.314
Wmlr=44.01e-3
Cp=40.0/Wmlr
Cv=Cp-Rcst/Wmlr
Gamma=Cp/Cv
Tinlet=430
Tinit=300
MFR=0.2

# Temperature graph
def create_Tfig():
    
    # Read the data
    df=pd.read_csv('monitor_Twall300/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5],names=['Time','Temp','Mass','Pres'])
    df['Tadia']=(df['Mass'].iloc[0]*df['Temp'].iloc[0]+Gamma*Tinlet*(df['Mass']-df['Mass'].iloc[0]))/df['Mass']
    
    # Put in fig
    Tfig=go.Figure()
    Tfig.add_trace(go.Scatter(name='0D adiabatic model',x=df['Time']/60,y=df['Tadia'],mode='lines',showlegend=True,line=dict(color='firebrick',width=2,dash='dot')))
    Tfig.add_trace(go.Scatter(name='NGA2 with Twall=300K',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(color='navy',width=2)))
    Tfig.update_layout(width=800,height=600)
    Tfig.update_xaxes(title_text='Time (min)',title_font_size=24,tickfont_size=24)
    Tfig.update_yaxes(title_text='Temperature (K)',title_font_size=24,tickfont_size=24,range=[280,450])
    Tfig.add_shape(type='line',x0=0,y0=Tinit,x1=df['Time'].iloc[-1]/60,y1=Tinit,line_color='black')
    Tfig.add_annotation(x=7.5,y=Tinit-7,text='Tinit',showarrow=False,font_size=16,font_color='black')
    Tfig.add_shape(type='line',x0=0,y0=Tinlet,x1=df['Time'].iloc[-1]/60,y1=Tinlet,line_color='green')
    Tfig.add_annotation(x=7.5,y=Tinlet+7,text='Tinlet',showarrow=False,font_size=16,font_color='green')
    Tfig.update_layout(legend=dict(font=dict(size=14)))
    
    # Add temperature data
    dfnow=pd.read_csv('monitor/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5],names=['Time','Temp','Mass','Pres'])
    Tfig.add_trace(go.Scatter(name='NGA2 with adbiatatic walls',x=dfnow['Time']/60,y=dfnow['Temp'],mode='lines',showlegend=True,line=dict(color='firebrick',width=2)))
    
    
    # Various debugging tests
    #dfnow['rhoCvT']=Cv*dfnow['Mass']*dfnow['Temp']
    #dfnow['rhoCpT']=Cp*dfnow['Mass']*dfnow['Temp']
    #dfnow['rhoCpT_model']=Cp*dfnow['Mass']*Tinit+(dfnow['Mass']-dfnow['Mass'].iloc[0])*Cp*Tinlet
    #dfnow['rhoCvT_model']=Cv*dfnow['Mass'].iloc[0]*Tinit+(dfnow['Mass']-dfnow['Mass'].iloc[0])*Cp*Tinlet
    # Put in fig
    #Tfig.add_trace(go.Scatter(name='Energy abiabatic wall',x=dfnow['Time']/60,y=dfnow['rhoCvT'],mode='lines',showlegend=True,line=dict(width=4)))
    #Tfig.add_trace(go.Scatter(name='model abiabatic wall',x=dfnow['Time']/60,y=dfnow['rhoCvT_model'],mode='lines',showlegend=True,line=dict(width=4)))
    
    return Tfig




# This is where we define the dashboard layout
def serve_layout():
    return html.Div(style={"margin-left": "15px"},children=[
    # Title of doc
    dcc.Markdown('''# Farther Farms Project'''),
    dcc.Markdown('''*NGA2 Dashboard written by O. Desjardins, last updated 02/06/2021*'''),
    # Intro
    dcc.Markdown('''
    ## Overview
    In this dashboard, we post-process the raw data generated by NGA2's pvessel
    case. This simulation is based on an experiment done by Farther Farms where
    a pressure vessel is filled with heated CO2.
    '''),
    # Imbibed volume over time
    dcc.Markdown(f'''
    ---
    ## Average temperature in the vessel
    The graph below shows the evolution of the average temperature inside the pressurized vessel.
    '''),
    #html.Div(create_Tfig(),style={'display':'none'}),
    dcc.Graph(id='Tgraph',figure=create_Tfig()),
])


# This is where we set the layout and run the server
app.layout = serve_layout
if __name__ == '__main__':
    app.run_server(debug=True)
