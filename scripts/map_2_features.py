#!/usr/bin/env python
'''
@author: manju
'''


import sys
import fileinput
import json

from tempfile import NamedTemporaryFile
from subprocess import Popen

if __name__=="__main__":

    claspre2_path = "/home/manju/workspace2/claspre2/build/release/bin/claspre"
    
    if len(sys.argv) > 2:
        print("Reading from: "+sys.argv[1])
        input_ = open(sys.argv[1],"r")
        indicators = sys.argv[2]
    else:
        print("Reading from: STDIN")
        input_ = sys.stdin
        indicators = sys.argv[1]
        
    indicators = map(int,indicators.split(","))

    outfile = NamedTemporaryFile(prefix="claspre_out_", suffix=".tmp", dir=".", delete=True)
    popen_= Popen([claspre2_path],stdin=input_,stdout=outfile)
    popen_.communicate()
    
    outfile.flush()
    outfile.seek(0)
    
    try:
        feature_dict = json.load(outfile)
    except:
        print("Features: -512")
        sys.exit(-1)
    
    preprocessing_feats = feature_dict["After_Preprocessing"]
    dynamic_feats = []
    index = 1
    
    while True:
        restart_feats = feature_dict.get("Dynamic-%d" % (index))
        if restart_feats:
            new_restart_feats = []
            for meta,f in restart_feats:
                new_restart_feats.append([meta+"_"+str(index),f])
            dynamic_feats.extend(new_restart_feats)
        else:
            break
        index += 1
    
    opt_feats = feature_dict.get("Optimization")
     
    flat_feats = []
    flat_names = []
    flat_feats.extend([y for (x,y) in preprocessing_feats])
    flat_names.extend([x for (x,y) in preprocessing_feats])
    flat_feats.extend([y for (x,y) in dynamic_feats])
    flat_names.extend([x for (x,y) in dynamic_feats])
    if opt_feats:
        flat_feats.extend([y for (x,y) in opt_feats])
        flat_names.extend([x for (x,y) in opt_feats])
        
    index = 0
    for i in indicators:
        if i == 1:
            print(flat_names[index])
        index += 1