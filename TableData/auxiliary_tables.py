#!/usr/bin/env python3
# Copyright 2025 Bernard Field

import os
import re
import itertools

from hermann_mauguin import layernumber_to_HM, texify_layer_HM, texify_rod_HM
from load_json import load_json

oblique_scanning = load_json(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ObliqueScanning.json'))

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
                          r"Penetration & $"+cvec+r"$ & Scanning & Scanning \\",
                          r"direction & & direction & group \\",
                          r"$[uv0]$ & & $\mathbf{d} = "+d_form+r"$ & $("+cvec+r",\mathbf{d},\mathbf{z})$ \\",
                          r"\hline",
                          r"\rule{0pt}{1.1em}\unskip"]) + '\n'
    else:
        mystr = '\n'.join([r"\begin{tabular}{|c|c|c|}",
                        r"\hline",
                        r"\rule{0pt}{1.1em}\unskip",
                        r"Penetration & Scanning & Scanning \\",
                        r"direction & direction & group \\",
                        r"$[uv0]="+cvec+r"$ & $\mathbf{d} = "+d_form+r"$ & $("+cvec+r",\mathbf{d},\mathbf{z})$ \\",
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

def format_orientation(v:list) -> str:
    v = [str(x) for x in v]
    if '/2' in v[0] and '/2' in v[1]:
        return '\\ensuremath{{[{}{}{}]/2}}'.format(v[0][:-2],v[1][:-2],v[2]).replace('-',r'\bar')
    else:
        return '\\ensuremath{{[{}{}{}]}}'.format(*v).replace('-',r'\bar')

def aux_table_json2tex_long(table:list) -> str:
    """
    Takes output of expand_aux_table_json, makes LaTeX table.
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
    # Compactify d_form if it's long.
    if d_form == "[(p+q)/2,(p-q)/2,0]":
        d_form = r"[p+q,p-q,0]/2" # There, that saves a few characters
    # Create the table header
    if c_column:
        mystr = '\n'.join([r"\begin{tabular}{|c|c|c|c|c|c|}",
                          r"\hline",
                          r"\rule{0pt}{1.1em}\unskip",
                          r"Penetration & $"+cvec+r"$ & Scanning & Scanning & Location & Penetration \\",
                          r"direction & & direction & group & $s\mathbf{d}$ & rod group \\",
                          r"$[uv0]$ & & $\mathbf{d} = "+d_form+r"$ & $("+cvec+r",\mathbf{d},\mathbf{z})$ & & $(\mathbf{d},\mathbf{z},"+cvec+r")$ \\",
                          r"\hline"]) + '\n'
    else:
        mystr = '\n'.join([r"\begin{tabular}{|c|c|c|c|c|}",
                        r"\hline",
                        r"\rule{0pt}{1.1em}\unskip",
                        r"Penetration & Scanning & Scanning & Location & Penetration \\",
                        r"direction & direction & group & $s\mathbf{d}$ & rod group \\",
                        r"$[uv0]="+cvec+r"$ & $\mathbf{d} = "+d_form+r"$ & $("+cvec+r",\mathbf{d},\mathbf{z})$ & & $(\mathbf{d},\mathbf{z},"+cvec+r")$ \\",
                        r"\hline"]) + '\n'
    # We shall first join adjacent rows with the same scanning group.
    new_table = []
    for row in table:
        if len(new_table) > 0 and new_table[-1]["H"] == row["H"]:
            # Matching rows. Append to prior
            old_row = new_table[-1]
            for key in ["orientation","c","d"]:
                old_row[key].append(row[key])
        else:
            # Different rows create new row.
            new_row = dict(orientation=[row["orientation"]],
                           c=[row["c"]],
                           d=[row["d"]],
                           H=row["H"],
                           general=row["general"],
                           general_linear_orbit=row["general_linear_orbit"],
                           special=row["special"])
            new_table.append(new_row)
    # Go over each set of distinct symmetries
    for group in new_table:
        # We'll build a list for each column we want to print, then we can just iterate.
        # Get locations
        ss = [ special["linear_orbit"] for special in group["special"] ]
        ss.append("$" + group["general_linear_orbit"].replace("1/2",r"\tfrac{1}{2}").replace("+/-",r"\pm") + "$")
        # Get rod groups
        rods = []
        for spec in group["special"]: # Special positions
            rodstr = texify_rod_HM(spec["rod"]["symbol"])
            if spec["rod"]["origin"] != [0,0,0]:
                # Get the origin shift (which is along c'), without quotes.
                rodstr += r" $[" + str(spec["rod"]["origin"][2]).replace("'","") + "]$"
            rodstr += r" \hfill R" + str(spec["rod"]["number"])
            if not spec["rod"]["standard_setting"]:
                rodstr += r"$^\prime$" # Prime if not standard setting.
            rods.append(rodstr)
        # Rod group, general position
        rodstr = texify_rod_HM(group["general"]["symbol"])
        if group["general"]["origin"] != [0,0,0]:
            # Get the origin shift (which is along c'), without quotes.
            rodstr += r" $[" + str(group["general"]["origin"][2]).replace("'","") + "]$"
        rodstr += r" \hfill R" + str(group["general"]["number"])
        if not group["general"]["standard_setting"]:
            rodstr += r"$^\prime$"
        rods.append(rodstr)
        # Orientations
        orientations = [_texify_uvpq(s) for s in group["orientation"]]
        # c form
        if c_column:
            # Wrap in maths
            cforms = ["$"+c+"$" for c in group["c"]]
        # Scanning directions
        ds = [_texify_uvpq(s) for s in group["d"]]
        # There's a few cases where this gets duplicated. We can cut the duplicates.
        for i in range(len(ds)-1,0,-1):
            if ds[i] == ds[i-1]:
                ds[i] = ""
        # Now, there's also a case where there's a very long orientation.
        # Let's line-break it.
        try:
            idx = orientations.index(_texify_uvpq("Odd u, even v OR even u, odd v"))
        except ValueError:
            pass
        else:
            orientations[idx] = _texify_uvpq("Odd u, even v")
            orientations.insert(idx+1, _texify_uvpq("OR even u, odd v"))
            ds.insert(idx+1, "")
            if c_column:
                cforms.insert(idx+1, "")
        # Finally, the scanning group, of which there is just one.
        Hs = [texify_layer_HM(group["H"]["symbol"]) + r" \hfill L" + str(group["H"]["number"]) + ("" if group["H"]["standard_setting"] else r"$^\prime$")]
        # Now we iterate over all of these.
        # Add a bit of padding at the top of the cell
        mystr += r'\rule{0pt}{1.1em}\unskip'+'\n'
        if c_column:
            iterables = [orientations, cforms, ds, Hs, ss, rods]
        else:
            iterables = [orientations, ds, Hs, ss, rods]
        for row in itertools.zip_longest(*iterables, fillvalue=""):
            mystr += " & ".join(row) + r"\\" + "\n"
        # Finished this one, add a horizontal rule
        mystr += r"\hline" + "\n"
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
        d_form = r"[p,q,0]/2"
    else:
        d_form = "[pq0]"
    c_form = "[u,v,0]" # Default form of c. Some cases are different.
    if H in [5,7]:
        # p11a and p112/a have 3 settings
        Hgroups = [dict(number=H, setting=str(i), origin=[0,0,0],
                      symbol=layernumber_to_HM(H,str(i),tex=False),
                      standard_setting=i==1)
                    for i in range(1,4)]
    else:
        # The other oblique groups have only one setting.
        Hgroup = dict(number=H, setting='1', origin=[0,0,0],
                      symbol=layernumber_to_HM(H,'1',tex=False),
                      standard_setting=True)
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
        table.append(dict(orientation='Odd u, even v OR even u, odd v', c=c_form, d='Even p,q OR odd p,q',
                          d_form=d_form, H=Hgroups[0]))
        table.append(dict(orientation='Odd u,v', c=c_form+'/2', d='Even p,q',
                          d_form=d_form, H=Hgroups[2]))
        table.append(dict(orientation='Odd u,v', c=c_form+'/2', d='Odd p,q',
                          d_form=d_form, H=Hgroups[1]))
    elif num in [10,13,18,22,26,35,47]:
        # aux_table_any_c, centring group without z-normal glide
        table.append(dict(orientation='Odd u,v', c=c_form+'/2', d='Even p,q OR odd p,q',
                          d_form=d_form, H=Hgroup))
        table.append(dict(orientation='Even u OR even v', c=c_form, d='Even p,q OR odd p,q',
                          d_form=d_form, H=Hgroup))
    else:
        # aux_table_any, primitive group without z-normal glide
        table.append(dict(orientation='Any u,v', c=c_form, d='Any p,q',
                          d_form=d_form, H=Hgroup))
    return table

def expand_aux_table_json(aux_table:list[dict], main_tables:list[dict]) -> list[dict]:
    """
    Takes a JSON auxiliary table, then expands it using data from main_tables
    to explicitly include location and rod group information.
    Modifies aux_table in-place.
    """
    # Get the auxiliary scanning group number
    H = aux_table[0]['H']['number']
    # Look up the relevant table
    reftable = main_tables[H-1]['table']
    # Add the other information.
    for row in aux_table:
        found = False
        # Find the matching group
        for refrow in reftable:
            # If find it, add the information.
            if row['H'] == refrow['H']:
                found = True
                for key in ['general', 'general_linear_orbit', 'special']:
                    row[key] = refrow[key]
                break
        if not found:
            print("WARNING! Unable to match H ", str(row['H']))
    return aux_table    
        

def end_table() -> str:
    return r"\end{tabular}" + "\n\n"
