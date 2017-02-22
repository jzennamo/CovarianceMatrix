import sys
from ROOT import gSystem
gSystem.Load("libCovarianceMatrix_CrossSectionMatrix")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing CrossSectionMatrix..."

sys.exit(0)

