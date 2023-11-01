#!/usr/bin/env python

import sys, os

if os.path.exists(sys.argv[1]):

	for line in open(sys.argv[1]):

		cmd=line.strip().strip("\"")
		os.system(cmd)
