#!/usr/bin/env python
import site
site.addsitedir('')
from optparse import OptionParser
import os
import sys
import VC

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = OptionParser()

    parser.add_option("-d", "--data-file", "--df",
                dest="data_file", default=None,
                help="The full path to the HDF5 Virtual California simulation output file.", 
                metavar="FILE"
                )
    (options, args) = parser.parse_args()
    
    data_file       = options.data_file

    sys_name = os.path.basename(data_file).split('.')[0]
    vc_sys = VC.VCSys(sys_name, data_file)
    vc_sys.webgui_checkCaches()

if __name__ == "__main__": 
    sys.exit(main())