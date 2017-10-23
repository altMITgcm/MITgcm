"""Parse verification outputs."""

import re
import glob
import os

def verification_parser(filename, threshold):

    # function should be given a threshold value for each sub test
    (directory, _) = os.path.split(filename)
    # how many additional tests are run with tweaks to this configuration
    num_exps = len(glob.glob(directory+'/input.*'))+1

    assert len(threshold)==num_exps

    # open the testreport output to assess pass/fail            
    with open(filename) as f:
        lines = f.readlines()

        first_match = True

        # extract lines from output
        for i, line in enumerate(lines):
            if line[:5] == '2 d e':
                if first_match:
                    # skip the first match, since it doesn't contain output,
                    # but set to false to catch next matches.
                    first_match = False
                else:
                    # save the line number where the output is found
                    output_line = i
                    break

        # loop through each of the subexperiments:
        for j in xrange(len(threshold)):
            test_results = lines[output_line+2+j]

            # split test_results into a list with values for each number. 
            # this uses spaces and the < > characters to separate the numbers.
            test_results = re.split('[ ><]',test_results)
            # ignore the Genmake, depend, make, and run checks, as well as
            # the "pass" or "fail" and test name at the end of the line
            test_results = test_results[4:-3]
            # convert to floats
            dp_similarity = []
            for i, x in enumerate(test_results):
                try:
                    dp_similarity.append(float(x))
                except ValueError:
                    pass

            assert all(elements > threshold[j] for elements in dp_similarity)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Check that verification simulation passed the test.')

    parser.add_argument('-filename', type=str, 
                        help='path to output file from the verification test')

    parser.add_argument('-threshold',nargs='+', type=int, default=15, 
                        help='number of decimal places of similarity required for test to pass. Requires a value for each sub test. Separate values with a space.')

    args = parser.parse_args()

    verification_parser(**vars(args))
