import argparse 
import os
import pandas as pd

ANNOTATIONS = 'annotations'
OUTPUT = 'output_sleuth_annotations'

def parse_input():
    '''
    Parses the commandline input, returns a dict
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, dest=ANNOTATIONS)
    parser.add_argument('-o', required=True, dest=OUTPUT)
    args = parser.parse_args()
    return vars(args)


if __name__ == '__main__':
    
    args = parse_input()
    original_annotation_file = args[ANNOTATIONS]
    annotations = pd.read_csv(original_annotation_file, sep='\t', names=['sample','condition'])
    wd = os.path.abspath(os.getcwd())
    annotations['path'] = annotations['sample'].apply(lambda x: os.path.join(wd, x))
    exists = annotations['path'].apply(lambda x: os.path.exists(x))
    annotations = annotations.loc[exists]
    annotations.to_csv(args[OUTPUT], sep='\t', index=False)