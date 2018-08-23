
from utils.filter_index import remove_redundancy
import sys
import argparse

parser = argparse.ArgumentParser(description='remove index redundancy')
parser.add_argument('file',
                    help='Input absolute path file containing examples')
args = parser.parse_args(sys.argv[1:])

remove_redundancy(args.file)