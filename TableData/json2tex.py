#!/usr/bin/env python3
# Copyright 2024 Bernard Field

import warnings
import json

from hermann_mauguin import rodnumber_to_HM, layernumber_to_HM, texify_rod_HM, texify_layer_HM, is_standard_setting
from load_json import load_json
from auxiliary_tables import aux_table_json, aux_table_json2tex, expand_aux_table_json, aux_table_json2tex_long

all_orbits = load_json('Orbits.json')

cvec = r"\mathbf{c}"

def process_single_layer_json(original_table:list, layer_num:int) -> dict:
    """
    Takes a single layer from LayerScanningTables.json, returns the corresponding processed layer JSON.
    """
    # Create the main table
    table = []
    for original_scan in original_table:
        new_scan = dict()
        H = dict(number=original_scan['H'],
                 origin=original_scan['H_origin'],
                 setting=original_scan['H_setting'])
        H['symbol'] = layernumber_to_HM(H['number'], H['setting'], tex=False)
        H['standard_setting'] = is_standard_setting(H['setting'])
        new_scan['H'] = H
        new_scan['orientation'] = original_scan['orientation']
        new_scan['d'] = original_scan['d']
        new_scan['general'] = dict(number=original_scan['general'],
                                   origin=original_scan['general_origin'],
                                   setting=original_scan['general_setting'],
                                   symbol=rodnumber_to_HM(original_scan['general'],original_scan['general_setting'],tex=False),
                                   standard_setting=is_standard_setting(original_scan['general_setting']))
        new_scan['general_linear_orbit'] = process_general_s(H['number'], H['setting'], H['origin'], tex=False)
        # Do special positions
        new_scan['special'] = []
        for original_spec in original_scan['special']:
            spec = dict()
            spec['rod'] = dict(number=original_spec['number'],
                               origin=original_spec['origin'],
                               setting=original_spec['setting'],
                               symbol=rodnumber_to_HM(original_spec['number'], original_spec['setting'], tex=False),
                               standard_setting=is_standard_setting(original_spec['setting']))
            spec['s'] = original_spec['s']
            spec['linear_orbit'] = process_spec_s(spec['s'], H['number'], H['setting'], H['origin'])
            new_scan['special'].append(spec)
        table.append(new_scan)
    # Create the auxiliary table
    auxiliary = aux_table_json(layer_num+1) # Because groups are 1-indexed while enumerate is 0-indexed
    # Put it in the big list
    return dict(group=dict(number=layer_num+1, symbol=layernumber_to_HM(layer_num+1,tex=False)),
                          table=table, auxiliary=auxiliary)

def process_layer_json(original_data:list) -> list:
    """
    Takes LayerScanningTables.json, creates LayerScanningTablesProcessed.json
    """
    processed = []
    for layer_num, original_table in enumerate(original_data):
        processed.append(process_single_layer_json(original_table, layer_num))
    # Add the expanded auxiliary table information
    for row in processed:
        expand_aux_table_json(row['auxiliary'], processed)
    return processed


def format_orientation(v:list) -> str:
    v = [str(x) for x in v]
    if '/2' in v[0] and '/2' in v[1]:
        return '\\ensuremath{{[{}{}{}]/2}}'.format(v[0][:-2],v[1][:-2],v[2]).replace('-',r'\bar')
    else:
        return '\\ensuremath{{[{}{}{}]}}'.format(*v).replace('-',r'\bar')


def write_table(f, layer:list, layer_num:int, do_auxiliary:bool=False, long_aux:bool=False):
    """
    Writes the scanning table for layer (with 1-based number layer_num) to file f in LaTeX
    From unprocessed JSON.
    """
    write_table_from_processed_json(f, process_single_layer_json(layer, layer_num),
                                    do_auxiliary=do_auxiliary, long_aux=long_aux)

def write_table_from_processed_json(f, layer:dict, do_auxiliary:bool=True, long_aux:bool=False):
    """
    Writes the scanning table for layer to file f in LaTeX
    Uses the format of process_layer_json.
    """
    # Write the table title
    f.write(r"\section*{" + texify_layer_HM(layer['group']['symbol']) + " No. "+str(layer['group']['number'])+"}\n\n")
    # Write the table header
    f.write('\n'.join([r"\begin{tabular}{|c|c|c|c|c|}",
            r"\hline",
            r"\rule{0pt}{1.1em}\unskip",
            r"Penetration & Scanning & Scanning & Location & Penetration \\",
            r"direction & direction $\mathbf{d}$ & group & $s\mathbf{d}$ & rod group \\",
            r"$[uv0]="+cvec+r"$ & & $("+cvec+r",\mathbf{d},\mathbf{z})$ & & $(\mathbf{d},\mathbf{z},"+cvec+r")$ \\"
            r"\hline"]) + '\n')
    for group in layer['table']:
        # Go over each set of distinct symmetries
        orientations = group["orientation"]
        ds = group["d"]
        specials = group["special"]
        # Set up progress flags (because must go line by line).
        done_general = False
        done_scanning_group = False
        do_scanning_group_origin = False
        # Add a bit of padding at the top of the cell
        f.write(r'\rule{0pt}{1.1em}\unskip'+'\n')
        # Print the entry line by line
        idx_orient = 0
        idx_special = 0
        while idx_orient < len(orientations) or not done_general or not done_scanning_group:
            if idx_orient < len(orientations):
                line = format_orientation(orientations[idx_orient])
                line += " & " + format_orientation(ds[idx_orient]) + " & "
                f.write(line)
                idx_orient += 1
            else:
                f.write(" & & ")
            if not done_scanning_group:
                if not do_scanning_group_origin:
                    line = texify_layer_HM(group["H"]["symbol"]) + r" \hfill L" + str(group["H"]["number"]) + ("" if group["H"]["standard_setting"] else r"$^\prime$") + " & "
                    if group["H"]["origin"] != [0,0,0]:
                        do_scanning_group_origin = True
                    else:
                        done_scanning_group = True
                else:
                    # Strip quotes from the origin
                    line = " $" + str(group["H"]["origin"]).replace("'","") + "$ & "
                    done_scanning_group = True
            else:
                line = " & "
            f.write(line)
            if idx_special < len(specials):
                spec = specials[idx_special]
                line = spec["linear_orbit"] + " & " + texify_rod_HM(spec["rod"]["symbol"])
                if spec["rod"]["origin"] != [0,0,0]:
                    # Get the origin shift (which is along c'), without quotes.
                    line += r" $[" + str(spec["rod"]["origin"][2]).replace("'","") + "]$"
                line += r" \hfill R" + str(spec["rod"]["number"])
                if not spec["rod"]["standard_setting"]:
                    line += r"$^\prime$" # Prime if not standard setting.
                idx_special += 1
            elif not done_general:
                # TeXify the general linear orbit
                gen_lin_orb = "$" + group["general_linear_orbit"].replace("1/2",r"\tfrac{1}{2}").replace("+/-",r"\pm") + "$"
                line = gen_lin_orb + " & " + texify_rod_HM(group["general"]["symbol"])
                if group["general"]["origin"] != [0,0,0]:
                    # Get the origin shift (which is along c'), without quotes.
                    line += r" $[" + str(group["general"]["origin"][2]).replace("'","") + "]$"
                line += r" \hfill R" + str(group["general"]["number"])
                if not group["general"]["standard_setting"]:
                    line += r"$^\prime$" # Prime if not standard setting.
                done_general = True
            else:
                line = r" & "
            f.write(line)
            f.write(r"\\" + "\n")
        f.write(r"\hline" + '\n')
    f.write(r"\end{tabular}"+"\n")
    if do_auxiliary:
        f.write(r"\nopagebreak" + "\n\n" + r"\noindent")
        if long_aux:
            f.write(aux_table_json2tex_long(layer["auxiliary"]))
        else:
            f.write(aux_table_json2tex(layer["auxiliary"]))
    else:
        f.write("\n")

class OrbitError(ValueError):
    pass

def _get_orbit(H:int, H_setting:str, H_origin:list) -> list:
    if H > 48:
        raise ValueError(f"Your scanning group {H} shouldn't be higher than 48.")
    # Identify the setting
    if H_setting == "1" or H_setting == "a":
        setting = "a"
    elif H_setting == "b" or H_setting == "2" or H_setting == "3":
        # While p11n (setting 2) and p11b (setting 3) are not the same group,
        # they have the same linear orbits. The key is whether or not the
        # orientation is parallel to the first vector.
        setting = "b"
    else:
        raise OrbitError(f"Can't process setting {H_setting}.")
    if H_origin == [0,0,0]:
        pass
    elif H_origin == ["1/4","1/4",0]:
        setting += "o"
    else:
        raise OrbitError(f"Can't process origin {H_origin}.")
    # Grab the set of orbits we need
    try:
        return all_orbits[H-1][setting]
    except KeyError:
        raise OrbitError(f"Can't find setting {setting} in orbit for L{H}.")

def process_spec_s(s_list:list, H:int, H_setting:str, H_origin:list) -> str:
    """
    For a list of special locations s, convert to string which presents the linear orbits.
    """
    try:
        # Get the s values which are in the same orbit.
        orbits = _get_orbit(H, H_setting, H_origin)
    except OrbitError as e:
        warnings.warn(str(e))
        # Fall-back to no special presentation.
        # Remove quotes around fractions, strip out the square braces from list
        return str(s_list).replace("'","").strip("[]")
    # We go through the s's and find which orbit they are in.
    indices = []
    for s in s_list:
        for i in range(len(orbits)):
            if s in orbits[i]:
                indices.append(i)
                break
    # What if we did not find an s in orbits? That's a problem.
    if len(s_list) != len(indices):
        warnings.warn(f"Error in matching locations {f} to orbits {orbits} of L{H}: got indices {indices}.")
        return str(s_list).replace("'","").strip("[]")
    # A check: do we have any s's in the same orbit that aren't in s_list?
    # In the event that this happens, I'll need to update my logic.
    for idx in indices:
        if len(orbits[idx]) != indices.count(idx):
            print(f"NOTICE! For H=L{H}, we found that s_list={s_list} did not span the whole orbit={orbits}")
            break
    # Consider our cases:
    # 1) All the s's are in the same orbit.
    if indices.count(indices[0]) == len(indices):
        # They are to be wrapped in a single set of square braces.
        # Because s_list is a list, it comes with square braces.
        return str(s_list).replace("'","")
    # 2) All the s's are in different orbits.
    elif len(set(indices)) == len(indices):
        # s as a list without square braces.
        return str(s_list).replace("'","").strip("[]")
    # 3) Mixed case.
    else:
        # Use a dictionary for sorting.
        # We group the s values into groups with the same orbit (identified by index)
        dd = dict()
        for s,idx in zip(s_list,indices):
            if idx in dd:
                dd[idx].append(s)
            else:
                dd[idx] = [s]
        # At the end of the day, we don't care which orbit they're in, just the grouping.
        # Use [1:-1] to strip the first and last character: extra square braces.
        return str(list(dd.values())).replace("'","")[1:-1]

def process_general_s(H:int, H_setting:str, H_origin:list, tex:bool=True) -> str:
    """
    For a given scanning group H, make a string for the general position's location/orbit.
    """
    try:
        orbit = _get_orbit(H, H_setting, H_origin)[-1]
    except OrbitError as e:
        warnings.warn(str(e))
        if tex:
            return "$s$"
        else:
            return "s"
    # As a personal trick, I generated the general orbit from the orbit of s=1/7
    # and stuck it in the last position in the orbits list.
    # (So this is fragile to implementation.)
    # There are a finite number of cases.
    if orbit == ["1/7"]:
        name = "s"
    elif orbit == ["1/7","6/7"]:
        name = r"[s, -s]"
    elif orbit == ["1/7","9/14"]:
        name = "[s, (s+1/2)]"
    elif orbit == ["1/7","5/14"]:
        name = "[s, (1/2-s)]"
    elif orbit == ["1/7","5/14","9/14","6/7"]:
        name = "[+/- s, (1/2 +/- s)]"
    else:
        warnings.warn(f"Cannot identify general orbit {orbit}.")
        name = "s"
    if tex:
        return "$" + name.replace("1/2",r"\tfrac{1}{2}").replace("+/-",r"\pm") + "$"
    else:
        return name



if __name__ == "__main__":
    tables_json = load_json('LayerScanningTables.json')
    
    processed_json = process_layer_json(tables_json)
    with open('LayerScanningTablesProcessed.json','w') as f:
        json.dump(processed_json, f, indent=0)

    with open('LayerScanningTables.tex','w') as f:
        for layer in processed_json:
            write_table_from_processed_json(f, layer, long_aux=True)

