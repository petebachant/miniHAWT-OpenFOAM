#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates the dynamic mesh dictionary

@author: pete
"""

import foampy

U = 1.0
R = 0.13
meantsr = 3.5

foampy.gen_dynmeshdict(U, R, meantsr, npoints=0, axis="(1 0 0)", direction=-1)
