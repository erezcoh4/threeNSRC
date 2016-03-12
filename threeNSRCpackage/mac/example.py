import sys
from ROOT import gSystem
gSystem.Load("libthreeNSRC_threeNSRCpackage")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing threeNSRCpackage..."

sys.exit(0)

