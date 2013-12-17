#!/usr/bin/env python
'''
@author: manju
'''


import sys
import fileinput
import json
import traceback

from tempfile import NamedTemporaryFile
from subprocess import Popen

if __name__=="__main__":

    claspre2_path = "/home/manju/workspace2/claspre2/build/release/bin/claspre"
    
    #===========================================================================
    # if len(sys.argv) > 1:
    #     print("Reading from: "+sys.argv[1])
    #     input_ = open(sys.argv[1],"r")
    # else:
    #===========================================================================
    print("Reading from: STDIN")
    input_ = sys.stdin

    outfile = NamedTemporaryFile(prefix="claspre_out_", suffix=".tmp", dir=".", delete=True)
    cmd = [claspre2_path]
    cmd.extend(sys.argv[1:])
    print(" ".join(cmd))
    popen_= Popen(cmd,stdin=input_,stdout=outfile)
    popen_.communicate()
    
    outfile.flush()
    outfile.seek(0)
    
    try:
        feature_dict = json.load(outfile)
    except:
        print("Features: -512")
        print("CRASHED")
        #traceback.print_exc()
        sys.exit(-1)
    
    preprocessing_feats = feature_dict["Static"]
    preprocessing_time = feature_dict["Static-Time"]
    dynamic_feats = []
    dynamic_times = []
    index = 1
    
    while True:
        restart_feats = feature_dict.get("Dynamic-%d" % (index))
        restart_time = feature_dict.get("Dynamic-Time-%d" % (index))
        if restart_feats:
            dynamic_feats.extend(restart_feats)
            dynamic_times.append(restart_time)
        else:
            break
        index += 1
    
    opt_pre_feats = feature_dict.get("Static_Optimization")
    opt_dynamic_feats = []
    index = 1
    
    while True:
        restart_feats = feature_dict.get("Optimization-%d" % (index))
        if restart_feats:
            opt_dynamic_feats.extend(restart_feats)
        else:
            break
        index += 1
     
    flat_feats = []
    flat_feats.extend([y for (x,y) in preprocessing_feats])
    flat_feats.extend([y for (x,y) in dynamic_feats])
    #===========================================================================
    # if opt_pre_feats:
    #     flat_feats.extend([y for (x,y) in opt_pre_feats])
    #     flat_feats.extend([y for (x,y) in opt_dynamic_feats])
    #===========================================================================
        
    times = [preprocessing_time]
    times.extend(dynamic_times)
    print("Features: %s" %(",".join(map(str,flat_feats))))
    print("Featuretimes: %s" %(",".join(map(str,times))))
    print("%s" %(feature_dict["Status"]))
    
