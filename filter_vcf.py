#!/usr/bin/python

from lib import testlib
import sys

testlib.filter_vcf(sys.argv[2], sys.argv[1])

