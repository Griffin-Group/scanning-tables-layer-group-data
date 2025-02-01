#!/usr/bin/env python3
# Copyright 2025 Bernard Field

import json

def load_json(fname):
    """Loads a JSON from file"""
    with open(fname, 'r') as f:
        return json.load(f)

