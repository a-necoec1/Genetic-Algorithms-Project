You need to have these imports inside of Python for it to work:

import time
import os
import argparse
import copy
import sys
import numpy as np

Running the program (ONLY WORKS IN COMMAND LINE):
Enter the name of the file you'd like to read with the number of chromosomes, number of generations, your choice of elitist or tournament, the percentage of what should be formed using selection in this generation , the selection type, and finally the mutation rate.

Command line example: python main.py genAlgData1.txt 10 5 tournament 0.6 uniform 0.3