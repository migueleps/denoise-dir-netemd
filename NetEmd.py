import os
import sys
from Orbcnt import count_directed, count_undirected
from ORBemd import compute_dist

def sanitize(user_input, allowed_inputs, index):
    if user_input not in allowed_inputs:
        raise ValueError("input [{}] in argument [{}] must be one of: {}\nPlease consult the README in https://github.com/migueleps/denoise-dir-netemd for the expected format.".format(user_input,index," OR ".join(allowed_inputs)))

subgraph_sizes = ["2","3","4","5"]
direction = ["undir","dir"]
types = ["allorbs", "weighted", "pca", "ica"]

indir = sys.argv[1]
input_size = sys.argv[2]
sanitize(input_size,subgraph_sizes,2)
input_direction = sys.argv[3]
sanitize(input_direction,direction,3)
if input_size == "2" and input_direction == "undir":
    raise ValueError("NetEmd not implemented for size 2 undirected. Please consult the README in https://github.com/migueleps/denoise-dir-netemd")
input_type = sys.argv[4]
sanitize(input_type,types,4)
n_threads = int(sys.argv[5])
if len(sys.argv) > 6:
    opt_arg = sys.argv[6]
else:
    opt_arg = ""

counting_func = {"dir": count_directed, "undir":count_undirected}

if __name__ == "__main__":
    counting_func[input_direction](indir,n_threads,int(input_size))
    compute_dist(indir,"",n_threads,input_direction == "dir",int(input_size),input_type,opt_arg)
