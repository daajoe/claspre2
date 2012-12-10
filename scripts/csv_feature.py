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

    feature_name = sys.argv[1]
    
    rest_feature_list_dic = {}
    for index in range(2,len(sys.argv)):
        input_ = open(sys.argv[index],"r")
    
        try:
            feature_dict = json.load(input_)
        except:
            print("Features: -512 - %s" % (sys.argv[index]))
            sys.exit(-1)
        input_.close()
            
        restart = 1
        while True:
            restart_feats = feature_dict.get("Dynamic-%d" % (restart))
            if not restart_feats:
                break
            for name, value in restart_feats:
                if name == feature_name:
                    f_list = rest_feature_list_dic.get(restart)
                    if not f_list:
                        rest_feature_list_dic[restart] = []
                        f_list = rest_feature_list_dic[restart]
                    f_list.append(value)
            restart += 1
                    
    for restart,f_list in rest_feature_list_dic.items():
        print("%d, %s" %(restart, ",".join(map(str,f_list))))
