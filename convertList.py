import os
import glob
from optparse import OptionParser

def copy(fpath):
    with open(fpath,"r") as f:
        for line in f:
            cmd = "xrdcp %s ."%line.strip()
            print(cmd)
            os.system(cmd)

def makelist():
    flist = glob.glob("./displacedJet*.root")
    fout = "local_list.txt"
    with open(fout,"w") as f:
        for line in flist:
            f.write(line+"\n")

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', '--inputList', dest='inputList', default='./lists/test.txt', metavar='inputList')
    (options, args) = parser.parse_args()

    copy(options.inputList)
    makelist()
