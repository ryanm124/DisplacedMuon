import os
import sys
from itertools import islice
submit = 'universe = vanilla\n' ## writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = submit01.out\n'
submit += 'error = submit01.err\n'
submit += 'log = submit01.log\n'
submit += '+JobFlavour = "tomorrow"\n' 
submit += 'queue 1\n'
submitName = 'condorFiles/submit01.sub'
sub1 = open(submitName,'w')
sub1.write(submit+'\n')
sub1.close() ##finish writing .sub file
create = '#!/bin/bash\n' ##writng .sh file
create += 'export PROJECT=/afs/cern.ch/user/r/rmccarth/private/dispVert/DisplacedMuon/\n' 
create += 'cd $PROJECT\n' ## go to the src directly 
create += 'export X509_USER_PROXY=/afs/cern.ch/user/r/rmccarth/x509up_u120948\n' ## exporting the proxy
create += './RunAll'
createName = 'condorFiles/submit01.sh'
sub2 = open(createName,'w')
sub2.write(create+'\n')
sub2.close() ## finish writing .sh file
os.system('chmod 755 '+createName) ## make .sh file executable
os.system('condor_submit '+ submitName+' executable='+createName) ## submit the job using condor_submit command

