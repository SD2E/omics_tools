import plotly.graph_objs as go
from plotly.offline import iplot
import plotly.figure_factory as ff
import numpy as np

def plot_interactive(dataframe):
    # Create Top Dendrogram
    data_array = dataframe.fillna(value=0)
    figure = ff.create_dendrogram(data_array, orientation='bottom')
    dendro_leaves = figure['layout']['xaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    labels = data_array.index
    figure = ff.create_dendrogram(data_array, orientation='bottom', labels=labels)
    #labels = [go_kegg_to_func[j] for j in data_array.columns]
    figure['layout']['yaxis']['ticktext'] = np.asarray(labels)

    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    data_array_T = data_array.T
    labels = data_array_T.index
    #labels = [go_kegg_to_func[j] for j in data_array_T.index]
    dendro_side = ff.create_dendrogram(data_array_T, orientation='right', )
    dendro_leaves_T = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves_T = list(map(int, dendro_leaves_T))
    figure['layout']['yaxis']['tickvals'] = np.asarray(dendro_side['layout']['yaxis']['tickvals'])

    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        figure.add_trace(data)

    # Create Heatmap
    ## reindex according to heirarchical clustering
    heat_data = data_array.iloc[dendro_leaves, :]
    heat_data = heat_data.T
    heat_data = heat_data.iloc[dendro_leaves_T, :]
    heatmap = [
        go.Heatmap(
            x=dendro_leaves,
            y=dendro_leaves_T,
            z=heat_data,
            zmin=-5.0,
            zmax=5,
            colorscale=[
                [0, 'rgb(255, 0, 0)'],
                [0.25, 'rgb(255, 204, 204)'],
                [0.5, 'rgb(255, 255, 255)'],
                [0.75, 'rgb(204, 255, 204)'],
                [1, 'rgb(0, 255, 0)']
            ],
            colorbar=dict(
                tickmode='array',
                tickvals=[-4, -2, 0, 2, 4],
                ticktext=['4', '2', '<b>-log10</b>',
                          '-2', '-4'],
                x=-0.15
            )
        )
    ]
    heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        figure.add_trace(data)

    # Edit Layout
    figure['layout'].update({'width': 1000, 'height': 1000,
                             'showlegend': False, 'hovermode': 'closest',
                             })

    figure['layout']['xaxis'].update({'domain': [.1, 0.9],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks': ""})

    figure['layout'].update({'xaxis2': {'domain': [0, .1],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""}})

    # ‘showticklabels’: True in figure[‘layout’][‘yaxis’].update
    figure['layout']['yaxis'].update({'domain': [0, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      # 'showticklabels': False,
                                      'side': 'right',
                                      'ticks': ""})

    figure['layout'].update({'yaxis2': {'domain': [.825, 1],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""}})

    iplot(figure, show_link=True, link_text='Export to plot.ly')
    return
