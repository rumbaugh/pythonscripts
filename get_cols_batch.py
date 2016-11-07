import numpy as np
import sys

def get_cols_batch(cr,colnum):
    test = get_col_names(cr)
    if colnum > len(test): sys.exit("More columns than column names in get_cols_batch")
    temparr = []
    for igcb in range(0,colnum):
        tempholder = get_colvals(cr,'col%i'%(igcb+1))
        temparr.append(tempholder)
        temptup = tuple(temparr)
    return temptup
        
