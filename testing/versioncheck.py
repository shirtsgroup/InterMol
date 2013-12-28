#!/bin/env python

import sys, os, subprocess, re

def check_Gromacs():
    try:
        output = subprocess.Popen('grompp_d -version', shell=True, stdout=subprocess.PIPE)
    except:
        print "Gromacs is not installed"
        return 1
    versionRegEx = re.compile(":-\)  VERSION (?P<version>.*)  \(-:")
    for line in output.stdout:
        version = versionRegEx.match(line.strip())
        if version:
            return version.group('version')  
 

if __name__ == '__main__':
    print "Welcome to check version"
    gmxversion = check_Gromacs()
    
    
