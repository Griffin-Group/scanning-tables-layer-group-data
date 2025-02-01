#!/usr/bin/env python3
# Copyright 2024 Bernard Field

from load_json import load_json

names = load_json('groupnames.json')

def texify_rod_HM(name:str) -> str:
    return r'\ensuremath{' + name.replace('-',r'\bar').replace('p',r'\mathscr{p}') + r'}'

def rodnumber_to_HM(i:int, setting:str="1", tex:bool=True) -> str:
    """Rod group number (1-based) to HM symbol, LaTeXified"""
    if (3 <= i <= 7 or 17 <= i <= 19 or i == 22) and setting == "bac":
        # Special setting-dependent logic
        # We've swapped the a and b vectors.
        if i == 6:
            name = "p12/m1"
        elif i == 7:
            name = "p12/c1"
        else:
            # Swap the 1st and 2nd element in the HM symbol.
            name_tmp = names["rod"][i-1]
            name = name_tmp[0] + name_tmp[2] + name_tmp[1] + name_tmp[3:]
    else:
        name = names["rod"][i-1]
    if tex:
        return texify_rod_HM(name)
    else:
        return name

def texify_layer_HM(name:str) -> str:
    return r'\ensuremath{' + name.replace('-',r'\bar') + r'}'

def layernumber_to_HM(i:int, setting:str="1", tex:bool=True) -> str:
    """Layer group number (1-based) to HM symbol, LaTeXified"""
    if (i == 5 or i == 7) and setting != "1":
        name = names["layer"][i-1]
        if setting == "2":
            name = name.replace('a','n')
        elif setting == "3":
            name = name.replace('a','b')
    elif (8 <= i <= 48) and setting == 'b':
        # Swap the 1st and 2nd element in the HM symbol, and transform a <-> b.
        name_tmp = names["layer"][i-1]
        # Break apart the HM symbol into its components
        elements = [name_tmp[0],'','','']
        j = 1
        for i in range(1,len(name_tmp)):
            elements[j] += name_tmp[i]
            # Increment to next element if we don't have a connector coming up.
            if name_tmp[i] != '_' and name_tmp[i] != '/' and (i+1 == len(name_tmp) or (name_tmp[i+1] != '_' and name_tmp[i+1] != '/')):
                j += 1
        # Interchange a and b
        for i in range(1,len(elements)):
            if 'a' in elements[i]:
                elements[i] = elements[i].replace('a','b')
            elif 'b' in elements[i]:
                elements[i] = elements[i].replace('b','a')
        # Swap 1st and 2nd element
        name = elements[0] + elements[2] + elements[1] + elements[3]
    else:
        name = names["layer"][i-1]
    if tex:
        return texify_layer_HM(name)
    else:
        return name

def is_standard_setting(setting:str) -> bool:
    """
    Is the setting a standard setting (1, a, abc)
    (I ignore the ones with multiple origin choices because those are captured separately.)
    """
    return (setting == "abc" or setting == "1" or setting == "a")