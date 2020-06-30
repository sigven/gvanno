#!/usr/bin/env python

import os

wd = os.getcwd()

cmd_grch37 = "python gvanno.py " + os.path.join(wd,"examples","example.grch37.vcf.gz") + " " + wd + " " + wd  + " grch37 gvanno.toml example --container docker --no_vcf_validate --force_overwrite"
cmd_grch38 = "python gvanno.py " + os.path.join(wd,"examples","example.grch38.vcf.gz") + " " + wd + " " + wd  + " grch38 gvanno.toml example --container docker --no_vcf_validate --force_overwrite"

os.system(cmd_grch37)
os.system(cmd_grch38)


