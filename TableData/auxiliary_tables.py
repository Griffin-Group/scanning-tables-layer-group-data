#!/usr/bin/env python3
# Copyright 2024 Bernard Field

import re

from hermann_mauguin import layernumber_to_HM, texify_layer_HM
from load_json import load_json

oblique_scanning = load_json('ObliqueScanning.json')

cvec = r"\mathbf{c}"

def aux_table(num:int, H:int|None=None) -> str:
    """
    Auxiliary table for scanned group num and scanning group H.
    If H not provided, read from data.
    """
    data = aux_table_json(num, H)
    return aux_table_json2tex(data)

def aux_table_json2tex(table:list) -> str:
    """
    Takes output of aux_table_json, makes LaTeX table.
    """
    # Determine if any c's ore not [u,v,0]
    # Only in this case shall we make a dedicated column for c.
    c_column = False
    for row in table:
        if row['c'] != '[u,v,0]':
            c_column = True
            break
    # Grab d_form. It's the same for all rows
    d_form = table[0]['d_form']
    # Create the table header
    if c_column:
        mystr = '\n'.join([r"\begin{tabular}{|c|c|c|c|}",
                          r"\hline",
                          r"\rule{0pt}{1.1em}\unskip",
                          r"Orientation & $"+cvec+r"$ & Scanning & Scanning \\",
                          r"$[uv0]$ & & direction & group \\",
                          r" & & $\mathbf{d} = "+d_form+r"$ & $("+cvec+r",\mathbf{d},\mathbf{z})$ \\",
                          r"\hline",
                          r"\rule{0pt}{1.1em}\unskip"]) + '\n'
    else:
        mystr = '\n'.join([r"\begin{tabular}{|c|c|c|}",
                        r"\hline",
                        r"\rule{0pt}{1.1em}\unskip",
                        r"Orientation & Scanning & Scanning \\",
                        r"$[uv0]="+cvec+r"$ & direction & group \\",
                        r" & $\mathbf{d} = "+d_form+r"$ & $("+cvec+r",\mathbf{d},\mathbf{z})$ \\",
                        r"\hline",
                        r"\rule{0pt}{1.1em}\unskip"]) + '\n'
    # Go over each line
    for idx, row in enumerate(table):
        # First column: orientation
        mystr += _texify_uvpq(row["orientation"]) + " & "
        if c_column:
            # Next column: form of c, wrapped in maths
            mystr += "$" + row["c"] + "$ & "
        # Next column: d
        if (idx > 0) and (table[idx-1]["d"] == row["d"]) and (table[idx-1]["H"] == row["H"]):
            # This is a row with the same d and H as the row above. Don't copy myself.
            mystr += "& "
        else:
            mystr += _texify_uvpq(row["d"]) + " & "
        # Final column: scanning group
        if (idx > 0) and (table[idx-1]["H"] == row["H"]):
            # This is a row with the same H above. Don't copy myself.
            pass # Add empty string.
        else:
            mystr += texify_layer_HM(row["H"]["symbol"]) + r" \hfill L" + str(row["H"]["number"])
        # End this row.
        if (idx == len(table)-1) or (table[idx+1]["H"] != row['H']):
            mystr += r"\\" + '\n' + r"\hline" + '\n'
            if idx < len(table)-1:
                # Tweak line spacing if not at the end.
                mystr += r"\rule{0pt}{1.1em}\unskip" + "\n"
        else:
            # The next row has the same scanning group, so don't draw a line.
            mystr += r"\\" + '\n'
    # End the table.
    mystr += end_table()
    return mystr
    
def _texify_uvpq(name:str) -> str:
    """Adds TeX maths to orientation or d strings"""
    pattern = re.compile(r'(\(p \+/- q\)|\b(?:u,v|u|v|p,q|p|q)\b)')
    # This pattern captures any instances of "(p +/- q)", "u,v", "p,q", "u", "v", "p", "q" that aren't embedded in other words.
    name = re.sub(pattern, r"$\1$", name)
    # Wrap those elements in maths delimiters.
    return name.replace("+/-", r"\pm") # Change +/- to LaTeX macro.

def aux_table_json(num:int, H:int|None=None) -> list:
    """
    JSON version of auxiliary table for scanned group num and scanning group H.
    If H not provided, read from data.
    """
    if H is None:
        H = oblique_scanning[num-1]
    table = []
    # Get some of the shared data
    if num in [10,13,18,22,26,35,36,47,48]: # Centering groups
        d_form = r"[(p+q)/2,(p-q)/2,0]"
    else:
        d_form = "[pq0]"
    c_form = "[u,v,0]" # Default form of c. Some cases are different.
    if H in [5,7]:
        # p11a and p112/a have 3 settings
        Hgroups = [dict(number=H, setting=str(i), origin=[0,0,0],
                      symbol=layernumber_to_HM(H,str(i),tex=False)) for i in range(1,4)]
    else:
        # The other oblique groups have only one setting.
        Hgroup = dict(number=H, setting='1', origin=[0,0,0],
                      symbol=layernumber_to_HM(H,'1',tex=False))
    # Create the tables
    if num in [5,7,31,33,38,41,43,45]:
        # aux_table_a, z-normal a-glide
        table.append(dict(orientation='Odd u, even v', c=c_form, d='Odd q',
                          d_form=d_form, H=Hgroups[0]))
        table.append(dict(orientation='Any u, odd v', c=c_form, d='Even q',
                          d_form=d_form, H=Hgroups[2]))
        table.append(dict(orientation='Any u, odd v', c=c_form, d='Odd q',
                          d_form=d_form, H=Hgroups[1]))
    elif num in [28,30]:
        # aux_table_b, z-normal b-glide
        table.append(dict(orientation='Even u, odd v', c=c_form, d='Odd p',
                          d_form=d_form, H=Hgroups[0]))
        table.append(dict(orientation='Odd u, any v', c=c_form, d='Even p',
                          d_form=d_form, H=Hgroups[2]))
        table.append(dict(orientation='Odd u, any v', c=c_form, d='Odd p',
                          d_form=d_form, H=Hgroups[1]))
    elif num in [32,34,39,42,46,52,62,64]:
        # aux_table_n, z-nomral n-glide
        table.append(dict(orientation='Odd u, odd v', c=c_form, d='Any p,q',
                          d_form=d_form, H=Hgroups[0]))
        table.append(dict(orientation='Even u OR even v', c=c_form, d='Odd p,q',
                          d_form=d_form, H=Hgroups[2]))
        table.append(dict(orientation='Even u, odd v', c=c_form, d='Even q',
                          d_form=d_form, H=Hgroups[1]))
        table.append(dict(orientation='Odd u, even v', c=c_form, d='Even p',
                          d_form=d_form, H=Hgroups[1]))
    elif num in [36,48]:
        # aux_table_c, centring group with z-normal glide
        table.append(dict(orientation='Odd u, even v OR even u, odd v', c=c_form, d='Any p,q',
                          d_form=d_form, H=Hgroups[0]))
        table.append(dict(orientation='Odd u,v', c=c_form+'/2', d='Even (p +/- q)',
                          d_form=d_form, H=Hgroups[2]))
        table.append(dict(orientation='Odd u,v', c=c_form+'/2', d='Odd (p +/- q)',
                          d_form=d_form, H=Hgroups[1]))
    elif num in [10,13,18,22,26,35,47]:
        # aux_table_any_c, centring group without z-normal glide
        table.append(dict(orientation='Odd u, odd v', c=c_form+'/2', d='Any p,q',
                          d_form=d_form, H=Hgroup))
        table.append(dict(orientation='Even u OR even v', c=c_form, d='Any p,q',
                          d_form=d_form, H=Hgroup))
    else:
        # aux_table_any, primitive group without z-normal glide
        table.append(dict(orientation='Any u,v', c=c_form, d='Any p,q',
                          d_form=d_form, H=Hgroup))
    return table


def end_table() -> str:
    return r"\end{tabular}" + "\n\n"
