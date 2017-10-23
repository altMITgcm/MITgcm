"""Parse verification outputs."""

import re

def verification_parser(filename, threshold):
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
                    test_results = lines[i+2]

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

        assert all(elements > threshold for elements in dp_similarity)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Check that verification simulation passed the test.')

    parser.add_argument('-filename', type=str, 
                        help='path to output file from the verification test')

    parser.add_argument('-threshold', type=int, default=15, 
                        help='number of decimal places of similarity required for test to pass')

    args = parser.parse_args()

    verification_parser(**vars(args))
