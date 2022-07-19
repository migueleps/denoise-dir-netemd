
# coding: utf-8

# In[1]:

import os
import sys


# In[ ]:

indir = sys.argv[1]
var=sys.argv[2]
n=sys.argv[3]
if var == "pca":
    os.system('python3 Orbcnt5.py '+indir+' '+n)
    os.system('python3 PCA_ORBemd.py '+indir+' '+indir+'-NetEmdOrb5'+' '+n+' '+sys.argv[4])
elif var == "pca3":
    os.system('python3 Orbcnt5.py '+indir+' '+n)
    os.system('python3 PCA_ORBemd3.py '+indir+' '+indir+'-NetEmdOrb5'+' '+n+" "+sys.argv[4])
elif var == "pca4":
    os.system('python3 Orbcnt5.py '+indir+' '+n)
    os.system('python3 PCA_ORBemd4.py '+indir+' '+indir+'-NetEmdOrb5'+' '+n+" "+sys.argv[4])
elif var == "ica":
    os.system('python3 Orbcnt5.py '+indir+' '+n)
    os.system('python3 ICA_ORBemd.py '+indir+' '+indir+'-NetEmdICA5'+' '+n+' '+sys.argv[4])
elif var == "ica3":
    os.system('python3 Orbcnt5.py '+indir+' '+n)
    os.system('python3 ICA_ORBemd3.py '+indir+' '+indir+'-NetEmdICA3'+' '+n+' '+sys.argv[4])
elif var == "ica4":
    os.system('python3 Orbcnt5.py '+indir+' '+n)
    os.system('python3 ICA_ORBemd4.py '+indir+' '+indir+'-NetEmdOrb5'+' '+n+' '+sys.argv[4])
elif var=='orb5':
    os.system('python3 Orbcnt5.py '+indir+' '+n)
    os.system('python3 ORBemd.py '+indir+' '+indir+'-NetEmdOrb5'+' '+n)
elif var=='orb4':
    os.system('python3 Orbcnt4.py '+indir+' '+n)
    os.system('python3 ORBemd4.py '+indir+' '+indir+'-NetEmdOrb4'+' '+n)
elif var=='ego4':
    os.system('python NEgocount.py '+indir+' '+n)
    os.system('python NEemd.py '+indir+' '+indir+'-NetEmdEgo4'+' '+n)
elif var=='spec':
    os.system('python LaplacianSpectra.py '+indir+' '+n)
    os.system('python EMDSpec.py '+indir+' '+n)
else:
    print 'Invalid argument...'


# In[ ]:



