# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from graph_virus import Virus
from pathlib import Path

from Bio import Phylo
import pandas as pd
from plotly.grid_objs import Column, Grid
from plotly.offline import download_plotlyjs, init_notebook_mode,  iplot, plot
init_notebook_mode(connected=True)

import plotly.figure_factory as ff
import numpy as np

app = dash.Dash()

virus_name = "zika"
#species = zip(list(([('avian', 'avian'), ('dengue', 'dengue'), ('ebola', 'ebola'), ('flu', 'flu'), ('lassa', 'lassa'), ('measles', 'measles'), ('mumps', 'mumps'), ('zika', 'zika')])))
species = ['avian', 'dengue', 'ebola', 'flu', 'lassa', 'measles', 'mumps', 'zika']


def compute_expensive_data(chemin):
    dir = dir + chemin
    return dir

def create_map():
    df_airports = pd.read_csv(
        'https://raw.githubusercontent.com/plotly/datasets/master/2011_february_us_airport_traffic.csv')
    df_airports.head()

    df_flight_paths = pd.read_csv(
        'https://raw.githubusercontent.com/plotly/datasets/master/2011_february_aa_flight_paths.csv')
    df_flight_paths.head()

    airports = [dict(
        type='scattergeo',
        locationmode='USA-states',
        lon=df_airports['long'],
        lat=df_airports['lat'],
        hoverinfo='text',
        text=df_airports['airport'],
        mode='markers',
        marker=dict(
            size=2,
            color='rgb(255, 0, 0)',
            line=dict(
                width=3,
                color='rgba(68, 68, 68, 0)'
            )
        ))]

    flight_paths = []
    for i in range(len(df_flight_paths)):
        flight_paths.append(
            dict(
                type='scattergeo',
                locationmode='USA-states',
                lon=[df_flight_paths['start_lon'][i], df_flight_paths['end_lon'][i]],
                lat=[df_flight_paths['start_lat'][i], df_flight_paths['end_lat'][i]],
                mode='lines',
                line=dict(
                    width=1,
                    color='red',
                ),
                opacity=float(df_flight_paths['cnt'][i]) / float(df_flight_paths['cnt'].max()),
            )
        )

    layout = dict(
        title='Feb. 2011 American Airline flight paths<br>(Hover for airport names)',
        showlegend=False,
        geo=dict(
            scope='north america',
            projection=dict(type='azimuthal equal area'),
            showland=True,
            landcolor='rgb(243, 243, 243)',
            countrycolor='rgb(204, 204, 204)',
        ),
    )

    fig_map = dict(data=flight_paths + airports, layout=layout)
    return fig_map


def create_fig(tree_file, metadata_file):
    tree = Virus.read_treefile(tree_file)
    x_coords = Virus.get_x_coordinates(tree)
    y_coords = Virus.get_y_coordinates(tree)
    line_shapes = []
    Virus.draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords, y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []

    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        text.append(cl.name)

    df = Virus.read_metadata(metadata_file)
    df.columns
    nb_genome = len(df)
    print(nb_genome)

    graph_title = Virus.create_title(virus_name, nb_genome)
    intermediate_node_color='rgb(100,100,100)'

    NA_color={'Cuba': 'rgb(252, 196, 174)',#from cm.Reds color 0.2, ... 0.8
                 'Dominican Republic': 'rgb(201, 32, 32)',
                 'El Salvador': 'rgb(253, 202, 181)',
                 'Guadeloupe': 'rgb(253, 202, 181)',
                 'Guatemala': 'rgb(252, 190, 167)',
                 'Haiti': 'rgb(252, 145, 114)',
                 'Honduras': 'rgb(239, 66, 49)',
                 'Jamaica': 'rgb(252, 185, 161)',
                 'Martinique': 'rgb(252, 190, 167)',
                 'Mexico': 'rgb(247, 109, 82)',
                 'Nicaragua': 'rgb(249, 121, 92)',
                 'Panama': 'rgb(252, 185, 161)',
                 'Puerto Rico': 'rgb(252, 174, 148)',
                 'Saint Barthelemy': 'rgb(253, 202, 181)',
                 'USA': 'rgb(188, 20, 26)',
                 'Canada': 'rgb(188, 20, 26)',
                 'USVI': 'rgb(206, 36, 34)'
              }


    SAmer_color={'Brazil': 'rgb(21, 127, 59)',# from cm.Greens colors 0.2, 0.4, 0.6, 0.8
                 'Colombia': 'rgb(153, 213, 149)',
                 'Ecuador': 'rgb(208, 237, 202)',
                 'French Guiana': 'rgb(211, 238, 205)',
                 'Peru': 'rgb(208, 237, 202)',
                 'Suriname': 'rgb(206, 236, 200)',
                 'Venezuela': 'rgb(202, 234, 196)',
                 'Puerto Rico': 'rgb(201, 235, 199)',
                 'Argentina': 'rgb(203, 225, 185)'
                 }


    SAsia_color={'Singapore': '#0000EE',
                 'Vietnam': '#1E90FF',
                 'Malaysia': '#1E90AF',
                 'Philippines': '#1E90AE',
                 'Thailand': '#1E90AB',
                 'Myanmar': '#1E90AC',
                 'Cambodia': '#1E90AA',
                 'Indonesia': '#1E90AA'
                 }

    pl_SAsia=[[0.0, '#1E90FF'], [0.5, '#1E90FF'], [0.5, '#0000EE'], [1.0,'#0000EE' ]]


    Oceania_color={'American Samoa': 'rgb(209,95,238)',
                     'Fiji': 'rgb(238,130, 238)',
                     'French Polynesia': 'rgb(148,0,211)',
                     'Tonga': 'rgb(238,130, 238)',
                     'Australia': 'rgb(233,125, 235)',
                     'Micronesia': 'rgb(231,123, 235)'
                   }


    China_color={'China': 'rgb(255,185,15'}

    JapanKorea_color={'Japan': '#fcdd04'}

    SubsaharanAfrica_color={'Guinea': 'rgb(209,95,238)',
                             'Liberia': 'rgb(238,130, 238)',
                             'Sierra Leone': 'rgb(148,0,211)',
                             'Cote D Ivoire': 'rgb(145,0,209)',
                             'Angola': 'rgb(143,0,207)',
                             'Seychelles': 'rgb(145,10,217)',
                             'Comoros': 'rgb(141,5,203)'
                            }


    Africa_color={'Sudan': 'rgb(209,95,238)',
                     'Gambia': 'rgb(238,130, 238)',
                     'Nigeria': 'rgb(235,135, 233)',
                     'Mali': 'rgb(235,131, 229)'
                  }


    Europe_color={'France': 'rgb(209,95,238)',
                 'Germany': 'rgb(238,130, 238)',
                 'Italy': 'rgb(238,130, 238)',
                 'United Kingdom': 'rgb(238,130, 238)',
                 'Netherlands': 'rgb(148,0,211)',
                 'Spain': 'rgb(141,7,221)'
                  }

    country = []
    region = []
    color = [intermediate_node_color] * len(X)
    print(set(list(df['Region'])))
    print(set(list(df['Country'])))

    for k, strain in enumerate(df['Strain']):

        i = text.index(strain)

        text[i] = text[i] + '<br>Country: ' + '{:s}'.format(df.loc[k, 'Country']) + '<br>Region: ' + '{:s}'.format(
            df.loc[k, 'Region']) + \
                  '<br>Collection date: ' + '{:s}'.format(df.loc[k, 'Date']) + \
                  '<br>Journal: ' + '{:s}'.format(df.loc[k, 'Journal']) + '<br>Authors: ' + '{:s}'.format(
            df.loc[k, 'Authors'])
        country.append(df.loc[k, 'Country'])
        region.append(df.loc[k, 'Region'])
        if df.loc[k, 'Region'] == 'North America':
            color[i] = NA_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'South America':
            color[i] = SAmer_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Southeast Asia':
            color[i] = SAsia_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Oceania':
            color[i] = Oceania_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'China':
            color[i] = '#fecc00'
        elif df.loc[k, 'Region'] == 'Japan Korea':
            color[i] = '#dc7928'
        if df.loc[k, 'Region'] == 'Subsaharan Africa':
            color[i] = SubsaharanAfrica_color[df.loc[k, 'Country']]
        if df.loc[k, 'Region'] == 'Africa':
            color[i] = Africa_color[df.loc[k, 'Country']]
        if df.loc[k, 'Region'] == 'Europe':
            color[i] = Europe_color[df.loc[k, 'Country']]
        else:
            pass

    print(graph_title)
    axis = dict(showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title='' #y title
              )

    nodes = dict(type='scatter',
               x=X,
               y=Y,
               mode='markers',
               marker=dict(color=color,
                           size=5),
               text=text, #vignet information of each node
               hoverinfo='')

    layout = dict(title=graph_title,
                font=dict(family='Balto',size=14),
                width=1000,
                height=3000,
                autosize=False,
                showlegend=False,
                xaxis=dict(showline=True,
                           zeroline=False,
                           showgrid=False,
                           ticklen=4,
                           showticklabels=True,
                           title='branch length'),
                yaxis=axis,
                hovermode='closest',
                shapes=line_shapes,
                plot_bgcolor='rgb(250,250,250)',
                margin=dict(l=10)
               )

    fig = dict(data=[nodes], layout=layout)
    return fig


#TO DO validation file and directory exist
def create_paths_file(virus_name, level1="", level2="", level3=""):
    dir = "data/" + virus_name + "/"
    if level1 == "" and level2 == "" and level3 == "":
        tree_file = dir + "nextstrain_" + virus_name + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_metadata.csv"
        return tree_file, metadata_file
    elif level2 == "" and level3 == "":
        dir = dir + "/"+level1+"/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_metadata.csv"
        return tree_file, metadata_file
    elif level3 == "":
        dir = dir + "/" + level1 + "/"+level2+"/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_metadata.csv"
        return tree_file, metadata_file
    else:
        dir = dir + "/" + level1 + "/"+level2+"/"+level3+"/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_metadata.csv"
        return tree_file, metadata_file


tree_file, metadata_file = create_paths_file(virus_name, level1="", level2="", level3="")
fig = create_fig(tree_file, metadata_file)
fig_map = create_map()

def serve_layout():
    return html.Div([
        html.Div(
            className="row",
            children=[
                html.Div(
                    className="one columns"
                ),
                html.Div(
                    className="two columns",
                    children=[
                        html.Div(
                            children=html.Div([
                                html.H1(children='Criterion'),
                                html.H1(children=''),
                                html.H6(children='Dataset'),
                                dcc.Dropdown(
                                    id='my-dropdown1',
                                    options=[{'label': species[i], 'value': species[i]} for i in range(len(species))],
                                    value='zika',
                                ),
                                html.Div(id='output-container'),

                                html.Div(id='controls-container_mumps', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown2',
                                        options=[{'label': i, 'value': i} for i in ['global', 'na']],
                                        value='global',
                                    ),
                                ]),

                                html.Div(id='controls-container_dengue', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown3',
                                        options=[{'label': i, 'value': i} for i in ['all', 'denv1', 'denv2', 'denv3', 'denv4']],
                                        value='all',
                                    ),
                                ]),

                                html.Div(id='controls-container_lassa', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown4',
                                        options=[{'label': i, 'value': i} for i in ['s', 'l']],
                                        value='s',
                                    ),
                                ]),

                                html.Div(id='controls-container_avian', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown5',
                                        options=[{'label': i, 'value': i} for i in ['h7n9']],
                                        value='h7n9',
                                    ),
                                    dcc.Dropdown(
                                        id='my-dropdown6',
                                        options=[{'label': i, 'value': i} for i in ['ha', 'mp', 'na', 'ns', 'np', 'pa', 'pb2', 'pb1']],
                                        value='ha',
                                    ),
                                ]),

                                html.Div(id='controls-container_flu', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown7',
                                        options=[{'label': i, 'value': i} for i in ['h3n2', 'h1n1pdm', 'vic', 'yam']],
                                        value='h3n2',
                                    ),
                                    dcc.Dropdown(
                                        id='my-dropdown8',
                                        options=[{'label': i, 'value': i} for i in
                                                 ['ha', 'na']],
                                        value='ha',
                                    ),
                                    dcc.Dropdown(
                                        id='my-dropdown9',
                                        options=[{'label': i, 'value': i} for i in
                                                 ['2y', '3y', '6y', '12y']],
                                        value='3y',
                                    ),
                                ]),

                                html.H1(children=''),
                                html.H1(children=''),
                                html.H6(children='Date Range'),
                                dcc.RangeSlider(
                                    count=1,
                                    min=0,
                                    max=10,
                                    step=0.5,
                                    marks={
                                        0: '0 °F',
                                        3: '3 °F',
                                        5: '5 °F',
                                        7.65: '7.65 °F',
                                        10: '10 °F'
                                    },
                                    value=[3, 7.65]
                                )
                            ])
                        )
                    ]
                ),
                html.Div(
                    className="five columns",
                    children=html.Div([
                        dcc.Graph(
                            id='right-top-graph',
                            figure=fig
                        ),
                        dcc.Graph(
                            id='right-bottom-graph',
                            figure={
                                'data': [{
                                    'x': [1, 2, 3],
                                    'y': [3, 1, 2],
                                    'type': 'bar'
                                }],
                                'layout': {
                                    'height': 400,
                                    'margin': {'l': 10, 'b': 20, 't': 0, 'r': 0}
                                }
                            }
                        )

                    ])
                ),
                html.Div(
                    className="four columns",
                    children=html.Div([
                        dcc.Graph(
                            id='right-mid-graph',
                            figure=fig_map
                        )
                    ])
                )
            ]
        )
    ])


app.layout = serve_layout()


app.css.append_css({
    'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
})


@app.callback(
    dash.dependencies.Output('output-container', 'children'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def update_output(value):
    global virus_name
    virus_name = value
    return 'You have selected "{}" virus'.format(value)


@app.callback(
    dash.dependencies.Output('controls-container_mumps', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def update_output(value):
    global virus_name
    virus_name = value
    if virus_name == "mumps":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_dengue', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def update_output(value):
    global virus_name
    virus_name = value
    if virus_name == "dengue":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_lassa', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def update_output(value):
    global virus_name
    virus_name = value
    if virus_name == "lassa":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_avian', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def update_output(value):
    global virus_name
    virus_name = value
    if virus_name == "avian":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_flu', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def update_output(value):
    global virus_name
    virus_name = value
    if virus_name == "flu":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('right-top-graph', 'figure'),
    [dash.dependencies.Input('my-dropdown1', 'value'),
     dash.dependencies.Input('my-dropdown2', 'value'),
     dash.dependencies.Input('my-dropdown3', 'value'),
     dash.dependencies.Input('my-dropdown4', 'value'),
     dash.dependencies.Input('my-dropdown5', 'value'), dash.dependencies.Input('my-dropdown6', 'value'),
     dash.dependencies.Input('my-dropdown7', 'value'), dash.dependencies.Input('my-dropdown8', 'value'), dash.dependencies.Input('my-dropdown9', 'value')])
def update_fig(value, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3):
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        tree_file, metadata_file = create_paths_file(virus_name, level1="", level2="", level3="")
        return create_fig(tree_file, metadata_file)
    elif virus_name == "mumps":
        tree_file, metadata_file = create_paths_file(virus_name, level1=mumps, level2="", level3="")
        return create_fig(tree_file, metadata_file)
    elif virus_name == "dengue":
        tree_file, metadata_file = create_paths_file(virus_name, level1=dengue, level2="", level3="")
        return create_fig(tree_file, metadata_file)
    elif virus_name == "lassa":
        tree_file, metadata_file = create_paths_file(virus_name, level1=lassa, level2="", level3="")
        return create_fig(tree_file, metadata_file)
    elif virus_name == "avian":
        tree_file, metadata_file = create_paths_file(virus_name, level1=avian_opt1, level2=avian_opt2, level3="")
        return create_fig(tree_file, metadata_file)
    elif virus_name == "flu":
        tree_file, metadata_file = create_paths_file(virus_name, level1=flu_opt2, level2=flu_opt2, level3=flu_opt3)
        return create_fig(tree_file, metadata_file)


if __name__ == '__main__':
    app.run_server(debug=True, port=5557)
