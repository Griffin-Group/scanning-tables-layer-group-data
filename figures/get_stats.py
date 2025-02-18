#!/usr/bin/env python3
# Copyright 2025 Bernard Field

import argparse
import json

import plotly.express as px
import plotly.graph_objects as go

import numpy as np

def load_json(fname):
    """Loads a JSON from file"""
    with open(fname, 'r') as f:
        return json.load(f)

def scanning_tables_condensed(table:list, aux:list) -> tuple[list,list]:
    """
    Get data from LayerScanningTables.json and ObliqueScanning.json
    Return condensed form showing allowed scanning groups and rod groups.
    Two lists: one from the layer group perspective, one from just the
    scanning groups.
    
    N.B. The groups are 1-indexed.
    Each entry in the returned list is 0-index of the layer group
    [{"H":{1,2,3,...}, "rod":{1,2,3,...}}, ...], [{1,2,3,...},...]
    """
    # Get the data
    L2H = scanning_flow_layer_to_H((table, aux))
    L2R = scanning_flow_layer_to_rod((table, aux))
    H2R = scanning_flow_H_to_rod((table, aux))
    # Put into expected form
    data = [{"H":a, "rod":b} for a,b in zip(L2H,L2R)]
    data_H = H2R
    return data, data_H

# Define layer group colors
layer_colors = px.colors.sample_colorscale('magenta',2,0.2,0.5) # Oblique/Triclinic
layer_colors += px.colors.sample_colorscale('purples',5,0.5) # Oblique/monoclinic
layer_colors += px.colors.sample_colorscale('blues',11,0.5) # Rectangular/monoclinic
layer_colors += px.colors.sample_colorscale('greens',30,0.5) # Rectangular/orthorhombic
layer_colors += px.colors.sample_colorscale('solar_r',16,0,0.5) # Square/tetragonal
layer_colors += px.colors.sample_colorscale('oranges',8,0.5) # Hexagonal/trigonal
layer_colors += px.colors.sample_colorscale('reds',8,0.5) # Hexagonal/hexagonal

# Define rod group colors
rod_colors = layer_colors[0:2] # Triclinic
rod_colors += px.colors.sample_colorscale('purples',5,0.5) # Monoclinic (inclined)
rod_colors += px.colors.sample_colorscale('blues',5,0.5) # Monoclinic (orthogonal)
rod_colors += px.colors.sample_colorscale('greens',10,0.5) # Orthorhombic

# Define greyscale friendly molecular foundry themed colours for crystal classes.
class_colors = ['#041c2c', '#b7312c', '#0076a8', '#e07e3c', '#8db7ca', '#f0b323', '#fbd872']
# Dark navy, dark red, slightly-dark blue, orange, slightly-light blue, sun(yellow), light yellow
# These colours are in order of increasing brightness. That is, sequential if convert to greyscale.
# (It isn't quite sequential in saturation, but that's a tricky one to balance while also
# conforming to the brand guideline pallete. But it doesn't *have* to be sequential, because it is
# categorical data.)
# layer group crystal classes.
layer_type_colors = [class_colors[0]]*2 # Oblique/Triclinic
layer_type_colors += [class_colors[1]]*5 # Oblique/monoclinic
layer_type_colors += [class_colors[2]]*11 # Rectangular/monoclinic
layer_type_colors += [class_colors[3]]*30 # Rectangular/orthorhombic
layer_type_colors += [class_colors[4]]*16 # Square/tetragonal
layer_type_colors += [class_colors[5]]*8 # Hexagonal/trigonal
layer_type_colors += [class_colors[6]]*8 # Hexagonal/hexagonal

def scanning_flow_layer_to_H(table:tuple) -> list:
    """
    From (scanning table,oblique scanning), get links from each scanned group
    to scanning group.

    List of sets, [{1,2,3,...},...]
    where index is layer group of layer (0-based) and entries of sets are
    scanning group (1-based).
    """
    data = []
    for layer, aux_group in zip(*table):
        # Intialise empty data
        # Using sets automatically avoids duplication.
        layer_data = set()
        # Go over each orientation orbit
        for entry in layer:
            layer_data.add(entry["H"])
        # Go over auxiliary tables, if we are beyond inclined scanning
        if len(data) > aux_group:
            # Add the data from the auxiliary scanning
            layer_data.update(data[aux_group-1])
        data.append(layer_data)
    return data

def scanning_flow_layer_to_rod(table:tuple) -> list:
    """
    From (scanning table,oblique scanning), get links from each scanned group
    to sectional rod group.

    List of sets, [{1,2,3,...},...]
    where index is layer group of layer (0-based) and entries of sets are
    rod group (1-based).
    """
    # Set up data structure
    data = []
    # Go over each layer group
    for layer, aux_group in zip(*table):
        # Intialise empty data
        # Using sets automatically avoids duplication.
        layer_data = set()
        # Go over each orientation orbit
        for orient in layer:
            layer_data.add(orient["general"])
            # Go over each special position
            for spec in orient["special"]:
                layer_data.add(spec["number"])
        # Go over auxiliary tables, if we are beyond inclined scanning
        if len(data) > aux_group:
            # Add the data from the auxiliary scanning
            layer_data.update(data[aux_group-1])
        # Add the data to the record
        data.append(layer_data)
    return data

def scanning_flow_H_to_rod(table:tuple) -> list:
    """
    From (scanning table,oblique scanning), get links from each scanning group
    to sectional rod group.

    List of sets, [{1,2,3,...},...]
    where index is layer group of scanning group (0-based) and entries of sets
    are rod group (1-based).
    """
    # Set up data structure
    data = [set() for i in range(48)]
    # Go over each layer group
    for layer, aux_group in zip(*table):
        # Go over each orientation orbit
        for orient in layer:
            data[orient["H"]-1].add(orient["general"])
            # Go over each special position
            for spec in orient["special"]:
                data[orient["H"]-1].add(spec["number"])
        # Go over auxiliary tables, if we are beyond inclined scanning
        if len(data) > aux_group:
            # Add the data from the auxiliary scanning
            data[orient["H"]-1].update(data[aux_group-1])
    return data

def scanning_class_flow_layer_to_H(table:tuple) -> list[set[int]]:
    """
    From (scanning table,oblique scanning), get links from each class of
    scanned group to class of scanning group.

    1: Oblique/Triclinic (L1-L2)
    2: Oblique/monoclinic (L3-L7)
    3: Rectangular/monoclinic (L8-L18)
    4: Rectangular/orthorhombic (L19-L48)
    5: Square/tetragonal (L49-L64)
    6: Hexagonal/trigonal (L65-L72)
    7: Hexagonal/hexagonal (L73-L80)

    List of sets, [{1,2,3,...},...]
    where index is layer group of layer (0-based) and entries of sets are
    scanning group (1-based).
    """
    # Get the granular data at the layer group level.
    flow_layers = scanning_flow_layer_to_H(table)
    # Mapping from individual layer groups to crystal classes.
    layer_to_class = [0]*2 + [1]*5 + [2]*11 + [3]*30 + [4]*16 + [5]*8 + [6]*8
    # Initialise empty data.
    data = [set() for _ in range(7)]
    for i, layer in enumerate(flow_layers):
        for entry in layer:
            # We want things 1-indexed, so must do some conversions.
            data[layer_to_class[i]].add(layer_to_class[entry-1]+1)
    return data

def scanning_class_flow_layer_to_rod(table) -> list[set[int]]:
    # Get the granular data at the layer group level.
    flow_layers = scanning_flow_layer_to_rod(table)
    # Mapping from individual layer groups to crystal classes.
    layer_to_class = [0]*2 + [1]*5 + [2]*11 + [3]*30 + [4]*16 + [5]*8 + [6]*8
    # And rod groups to crystal classes
    # Triclinic, Monoclinic/Oblique, Monoclinic/Rectangular, Orthorhombic
    rod_to_class = [0]*2 + [1]*5 + [2]*5 + [3]*10
    # Initialise empty data.
    data = [set() for _ in range(7)]
    for i, layer in enumerate(flow_layers):
        for entry in layer:
            # We want things 1-indexed, so must do some conversions.
            data[layer_to_class[i]].add(rod_to_class[entry-1]+1)
    return data

def scanning_class_flow_H_to_rod(table) -> list[set[int]]:
    # Get the granular data at the layer group level.
    flow_layers = scanning_flow_H_to_rod(table)
    # Mapping from individual layer groups to crystal classes.
    layer_to_class = [0]*2 + [1]*5 + [2]*11 + [3]*30 + [4]*16 + [5]*8 + [6]*8
    # And rod groups to crystal classes
    # Triclinic, Monoclinic/Oblique, Monoclinic/Rectangular, Orthorhombic
    rod_to_class = [0]*2 + [1]*5 + [2]*5 + [3]*10
    # Initialise empty data.
    data = [set() for _ in range(7)]
    for i, layer in enumerate(flow_layers):
        for entry in layer:
            # We want things 1-indexed, so must do some conversions.
            data[layer_to_class[i]].add(rod_to_class[entry-1]+1)
    return data

def invert_flow(func):
    """
    Wrapper which inverts a list of sets indicating links
    """
    def inverted_flow(table) -> list:
        data = func(table)
        # Find maximum number
        max_num = max([max(x) for x in data])
        # Generate new data
        new_data = [set() for i in range(max_num)]
        # Iterate
        for i, entry in enumerate(data):
            for j in entry:
                new_data[j-1].add(i+1) # 0-based to 1-based
        return new_data
    return inverted_flow

def generic_sankey(table, columns:list, norm:str, column_flow_functions,
                   column_colors, show=True, labels=None) -> go.Figure:
    # Sort columns
    columns = columns.copy()
    columns.sort()
    # Configure
    ndims = len(column_flow_functions)
    column_lengths = [len(c) for c in column_colors]
    # Confirm sensible inputs
    assert len(column_colors) == ndims, f"Expected {ndims} column colors."
    for i in range(ndims):
        assert len(column_flow_functions[i]) == ndims, "column_flow_functions not square array."
    # Cache results of flow calculation
    column_flow = [[None]*ndims for i in range(ndims)]
    def get_flow(c1:int, c2:int):
        assert c1 != c2
        if column_flow[c1][c2] is None:
            column_flow[c1][c2] = column_flow_functions[c1][c2](table)
        return column_flow[c1][c2]
    def count_paths(index:int, ci:int) -> int:
        # Recursively count how many paths lead to the end.
        if ci == len(columns)-1:
            # You've hit the end, that is one path.
            return 1
        # Count how many paths lead from each of the nodes in the next column over to the end
        # Then sum up each of these for each node the current node is connected to.
        return sum([count_paths(i-1, ci+1) for i in get_flow(columns[ci],columns[ci+1])[index]])
    # Column heights of each node (initialised)
    column_heights = [[0]*L for L in column_lengths]
    # Normalise
    assert norm == "path" or norm == "column"
    # We have to figure out how high the initial columns are
    if norm == "path":
        column_heights[columns[0]] = [count_paths(i,0) for i in range(column_lengths[columns[0]])]
    elif norm == "column":
        column_heights[columns[0]] = [1] * column_lengths[columns[0]]
    assert len(column_heights[columns[0]]) == column_lengths[columns[0]]
    # Generate links
    # Initialise
    source = []
    target = []
    value = []
    # Plotly Sankey indexes all columns together, so we want to track the
    # number of nodes in prior columns
    prior_indices = 0
    for ci in range(len(columns)-1):
        # Flow from column 1 to column 2
        c1 = columns[ci]
        c2 = columns[ci+1]
        flow = get_flow(c1, c2)
        for a in range(len(flow)): # 0-index
            for b in flow[a]: # 1-index
                source.append(a + prior_indices)
                target.append(b-1 + prior_indices + column_lengths[c1])
                # Figure out the value of this link
                if norm == "path":
                    # Weighting: the number of links going *into* this node.
                    if ci == 0:
                        # Initial column
                        weight = 1
                    else:
                        weight = len(get_flow(c1, columns[ci-1])[a])
                    # Now we sum the paths from b to the end
                    val = weight * count_paths(b-1, ci+1)
                elif norm == "column":
                    # Weight by height divided by immediate paths out
                    val = column_heights[c1][a] / len(get_flow(c1,c2)[a])
                value.append(val)
                column_heights[c2][b-1] += val
        prior_indices += column_lengths[c1]
    # Build the node information
    # Colours
    colors = []
    for c in columns:
        colors += column_colors[c].copy()
    # Label
    if labels is None:
        labels = []
        for c in columns:
            labels += [str(i) for i in range(1,column_lengths[c]+1)]
    # Coordinates for plotting
    # They are figure coords, so between 0 and 1.
    # https://stackoverflow.com/questions/72749062/how-to-set-order-of-the-nodes-in-sankey-diagram-plotly says coords get buggy if exactly 0 or 1.
    x = []
    for i,c in enumerate(columns):
        x += [0.001 + 0.998 * i/(len(columns)-1)] * column_lengths[c]
    # y should scale roughly with size of the bars.
    # While plotly has some magic arrangements which can handle adding padding,
    # it tends to get ordering wrong if the bars overlap.
    y = np.empty((sum([column_lengths[c] for c in columns]),))
    prior_indices = 0
    for c in columns:
        heights = np.asarray(column_heights[c])
        # y is sum of prior node widths plus half current node width.
        y[prior_indices:prior_indices+column_lengths[c]] = (np.cumsum(heights) - heights/2) / heights.sum()
        prior_indices += column_lengths[c]
    # Because plotly is buggy, we need to delete x and y for empty nodes
    x = np.asarray(x)
    heights = []
    for c in columns:
        heights += column_heights[c]
    x = x[np.nonzero(heights)]
    y = y[np.nonzero(heights)]
    # Get link colours
    # Do it by source
    link_colors = [colors[i] for i in source]
    # Add transparency
    link_colors = [c.replace('rgb','rgba').replace(')',',0.5)') for c in link_colors]
    # Create plot
    fig = go.Figure(data=[go.Sankey(
        arrangement = "fixed",
        node = dict(
            label = labels,
            color = colors,
            y = y,
            x = x,
            pad = 1,
        ),
        link = dict(
            source = source,
            target = target,
            value = value,
            color = link_colors,
        ))])
    if show:
        fig.show()
    return fig

def layer_sankey(table:list, aux_table:list, columns:list=[0,1,2], norm:str="path"):
    """
    Sankey plot of scanning tables.

    Currently supported columns are: layer, H, rod.

    """
    column_flow_functions = [
            [None, scanning_flow_layer_to_H, scanning_flow_layer_to_rod],
            [invert_flow(scanning_flow_layer_to_H), None, scanning_flow_H_to_rod],
            [invert_flow(scanning_flow_layer_to_rod), invert_flow(scanning_flow_H_to_rod), None]
            ]
    column_colors = [layer_colors, layer_colors[0:48], rod_colors]
    generic_sankey((table, aux_table), columns, norm, column_flow_functions, column_colors)

def layer_class_sankey(table:list, aux_table:list, columns:list=[0,1,2],
                       norm:str="column", filename=None, width=330, height=None,
                       scale=1, font_family="Arial"):
    column_flow_functions = [
            [None, scanning_class_flow_layer_to_H, scanning_class_flow_layer_to_rod],
            [invert_flow(scanning_class_flow_layer_to_H), None, scanning_class_flow_H_to_rod],
            [invert_flow(scanning_class_flow_layer_to_rod), invert_flow(scanning_class_flow_H_to_rod), None]
            ]
    column_colors = [class_colors.copy(), class_colors[0:4], class_colors[0:4]]
    for c in column_colors:
        c[0] = "#d1dde6" # Replace dark navy with light grey (fog), so I can still read the labels.
    labels = ["Oblique/Triclinic","Oblique/Monoclinic","Rectangular/Monoclinic","Rectangular/Orthorhombic","Square/Tetragonal","Hexagonal/Trigonal","Hexagonal/Hexagonal"]
    labels += labels[0:4]
    labels[-1] = "Rectangular/<br>Orthorhombic"
    labels += ["Triclinic","Monoclinic/Oblique", "Monoclinic<br>/Rectangular","Orthorhombic"]
    if norm == "column":
        labels[-3] = "Monoclinic<br>/Oblique"
        labels[-6] = "Rectangular/<br>Monoclinic"
    fig = generic_sankey((table, aux_table), columns, norm, column_flow_functions, column_colors, show=False, labels=labels)
    # Decorate with column labels
    headers = ["Scanned group", "Scanning group", "Penetration<br>rod group"]
    annotations = []
    for i,c in enumerate(columns):
        annotations.append(dict(
            x=i/(len(columns)-1),
            y=1.12,
            xref="x domain",
            yref="paper",
            text=headers[c],
            showarrow=False,
            font=dict(size=10),
            ))
    fig.update_layout(annotations=annotations)
    # Hide the text shadow because it is ugly in PDF rendering.
    # It does no change when "none", but anything else makes the shadow go away.
    fig.update_traces(textfont_shadow="thisisinvalidsonoshadowappears")
    # Fix layout for saving.
    fig.update_layout(margin=dict(b=10,l=10,r=10,t=40), plot_bgcolor='#fff',
                      font=dict(size=10,color='#000',family=font_family))
    # Save
    # Plotly measures all distances in pixels.
    # According to my outputs and Preview, 100 px = 26.5 mm = 1.043 in.
    if filename:
        fig.write_image(filename, width=width, height=height, scale=scale)
    # Output
    fig.show()

def bar_chart(table, flow_function, name:str, column_colors:list,
              xlabel:str=None, width:float=330, height:float=200,
              font_family:str=None):
    """
    Create a bar chart
    
    table: data
    flow_function: gives list of x to y from table.
    name: save plot here
    column_colors: list of colours for each column
    xlabel: label for the x-axis
    """
    data = flow_function(table)
    counts = [len(x) for x in data]
    fig = go.Figure(data=[go.Bar(x=np.arange(1,len(counts)+1), y=counts,
                          marker_color=column_colors, marker_line_width=0)])
    fig.update_yaxes(title_text="Scanned groups", minor=dict(dtick=1,showgrid=True),
                     dtick=(int(max(counts)/40)+1)*5,
                     gridcolor='#999', title_standoff=5)
    if not xlabel is None:
        fig.update_xaxes(title_text=xlabel)
    fig.update_xaxes(dtick=int(len(counts)/20)+1, title_standoff=5)
    fig.update_layout(margin=dict(b=10,l=10,r=10,t=10), plot_bgcolor='#fff',
                      font=dict(size=10,color='#000',family=font_family))
    # Plotly measures all distances in pixels.
    # According to my outputs and Preview, 100 px = 26.5 mm = 1.043 in.
    fig.write_image(name, width=width, height=height, scale=1)
    fig.show()

def stacked_bar_chart(table, flow_function, name:str, group_colors:list,
                      xlabel:str=None, ylabel:str="Scanned groups",
                      width:float=330, height:float=200, scale:float=1,
                      font_family:str=None):
    """
    Create a bar chart of counts of y for each x, coloured by y.
    e.g. if flow_function is scanning groups to scanned groups,
    then the height of the bars is the number of scanned groups per scanning
    group, and the bars are coloured by the scanned groups which contribute to
    the counts.
    """
    fig = _stacked_bar_chart_core(table, flow_function, group_colors,
                                  xlabel=xlabel, ylabel=ylabel,
                                  font_family=font_family)
    # Plotly measures all distances in pixels.
    # According to my outputs and Preview, 100 px = 26.5 mm = 1.043 in.
    if name:
        fig.write_image(name, width=width, height=height, scale=scale)
    fig.show()

def stacked_bar_chart_with_legend(table, flow_function, filename:str, group_colors:list,
                                  xlabel:str=None, ylabel:str="Scanned groups",
                                  legend_names:list=None, legend_colors:list=None,
                                  width:float=330, height:float=250, scale:float=1,
                                  font_family:str=None):
    fig = _stacked_bar_chart_core(table, flow_function, group_colors,
                                  xlabel=xlabel, ylabel=ylabel,
                                  font_family=font_family)
    if legend_colors is not None:
        # Make legend
        if legend_names is None:
            # Have default names, just numbers
            legend_names = [str(i) for i in range(len(legend_colors))]
        for (name, color) in zip(legend_names, legend_colors):
            fig.add_trace(go.Bar(x=[None], y=[None], marker_color=color,
                                 marker_line_width=0, showlegend=True, name=name))
        # Show legend and adjust position to top-right of plot
        fig.update_layout(showlegend=True, legend=dict(
            yanchor="top", y=1, xanchor="right", x=1))
    # Save
    # Plotly measures all distances in pixels.
    # According to my outputs and Preview, 100 px = 26.5 mm = 1.043 in.
    if filename:
        fig.write_image(filename, width=width, height=height, scale=scale)
    fig.show()


def _stacked_bar_chart_core(table, flow_function, group_colors:list,
                            xlabel:str=None, ylabel:str=None,
                            font_family:str=None) -> go.Figure:
    """
    Helper function, returns Figure object for further manipulation

    Create a bar chart of counts of y for each x, coloured by y.
    e.g. if flow_function is scanning groups to scanned groups,
    then the height of the bars is the number of scanned groups per scanning
    group, and the bars are coloured by the scanned groups which contribute to
    the counts.
    """
    #data = flow_function(table)
    # Outer index: index (0-based) over x.
    # Inner list: list (1-based) of y's.
    # Invert it,
    x_data = invert_flow(flow_function)(table)
    # Find number of x values:
    nx = max([max(x) for x in x_data])
    # Find maximum y value
    ymax = max([len(x) for x in flow_function(table)])
    # Plot
    fig_data = [go.Bar(x=list(d), y=np.ones(len(d)), marker_color=group_colors[i],
                       marker_line_width=0.2, marker_line_color=group_colors[i],
                       showlegend=False)
                for i,d in enumerate(x_data)]
    # (A finite small marker_line_width ensures no hairline gaps between blocks)
    fig = go.Figure(data=fig_data)
    fig.update_yaxes(minor=dict(dtick=1,showgrid=True),
                     dtick=(int(ymax/40)+1)*5,
                     gridcolor='#999', title_standoff=5)
    if not xlabel is None:
        fig.update_xaxes(title_text=xlabel)
    if not ylabel is None:
        fig.update_yaxes(title_text=ylabel)
    fig.update_xaxes(dtick=int(nx/20)+1, title_standoff=5)
    fig.update_layout(margin=dict(b=10,l=10,r=10,t=10), plot_bgcolor='#fff',
                      font=dict(size=10,color='#000',family=font_family),
                      barmode='stack', showlegend=False)
    # Return to calling function
    return fig



def main(print_stats:bool=True, include_H:bool=True, norm:str="path",
         bars:bool=True, sankey:bool=True, stacked_bars:bool=True, font:str=None):
    scanning_table = load_json("../TableData/LayerScanningTables.json")
    oblique_scanning = load_json("../TableData/ObliqueScanning.json")
    # Display
    if print_stats:
        scanning_summary, H_summary = scanning_tables_condensed(scanning_table, oblique_scanning)
        for i,x in enumerate(scanning_summary):
            print(i+1, ": H: ", sorted(x["H"]), " rod: ", sorted(x["rod"]))
        for i,x in enumerate(H_summary):
            print("H=", i+1, " rod: ", sorted(x))
    if bars:
        bar_chart((scanning_table,oblique_scanning),
                  invert_flow(scanning_flow_layer_to_H),
                  "scanning_H_bars.pdf", layer_colors,
                  xlabel="Scanning group (layer IT number)")
        bar_chart((scanning_table,oblique_scanning),
                  invert_flow(scanning_flow_layer_to_rod),
                  "scanning_rod_bars.pdf", rod_colors,
                  xlabel="Penetration rod group (rod IT number)")
    if stacked_bars:
        legend_colors = class_colors
        legend_names = ['Oblique/Triclinic', 'Oblique/Monoclinic', 'Rectangular/Monoclinic',
                        'Rectangular/Orthorhombic', 'Square/Tetragonal', 'Hexagonal/Trigonal',
                        'Hexagonal/Hexagonal']
        
        stacked_bar_chart_with_legend((scanning_table, oblique_scanning),
                  invert_flow(scanning_flow_layer_to_H),
                  "scanning_H_stacked_bars.pdf", layer_type_colors,
                  xlabel="Scanning group (layer IT number)",
                  ylabel="Number of scanned groups",
                  legend_colors=legend_colors, legend_names=legend_names,
                  font_family=font)
        if font=="Times":
            height = 260
        else:
            height = 280
        stacked_bar_chart_with_legend((scanning_table,oblique_scanning),
                  invert_flow(scanning_flow_layer_to_rod),
                  "scanning_rod_stacked_bars.pdf", layer_type_colors,
                  xlabel="Penetration rod group (rod IT number)",
                  ylabel="Number of scanned groups",
                  legend_colors=legend_colors, legend_names=legend_names,
                  height=height, font_family=font)
    if include_H:
        columns = [0,1,2]
    else:
        columns = [0,2]
    if sankey:
        layer_sankey(scanning_table, oblique_scanning, columns=columns,
                     norm=norm)
        layer_class_sankey(scanning_table, oblique_scanning, columns=columns,
                           norm=norm, filename="layer_class_sankey.pdf",
                           height=280, font_family=font)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--print',action='store_true',
                        help="Print condensed table.")
    parser.add_argument('-H',action='store_true',
                        help="Include H in Sankey plot.")
    parser.add_argument('-n','--norm',type=str,
                        choices=["path","column"],
                        default="path",
                        help="Normalisation mode for Sankey.")
    parser.add_argument('-b','--bar', action='store_true',
                        help="Show bar plots.")
    parser.add_argument('--sankey', action='store_true',
                        help="Show Sankey plots.")
    parser.add_argument('--stacked', action='store_true',
                        help="Show stacked bar plots.")
    parser.add_argument('-f','--font',
                        help="Names of fonts to use for plots. Maybe be a comma-separated string for multiple options.")
    args = parser.parse_args()
    main(print_stats=args.print, include_H=args.H, norm=args.norm,
         bars=args.bar, sankey=args.sankey, stacked_bars=args.stacked,
         font=args.font)
