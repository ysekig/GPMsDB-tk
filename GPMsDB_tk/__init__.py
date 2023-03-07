#!/usr/bin/env python

import os

# Set the module version.
with open(os.path.join(__path__[0], 'VERSION'), 'r') as f:
    __version__ = f.readline().strip()
