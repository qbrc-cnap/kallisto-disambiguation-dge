import argparse 
import pandas as pd

INPUT = 'input_table'
OUTPUT = 'output_table'

def parse_input():
    '''
    Parses the commandline input, returns a dict
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, dest=INPUT)
    parser.add_argument('-o', required=True, dest=OUTPUT)
    args = parser.parse_args()
    return vars(args)


if __name__ == '__main__':
    
    args = parse_input()
    original_nc_file = args[INPUT]
    df = pd.read_table(original_nc_file, sep='\t')
    dd=pd.pivot_table(df, values='tpm', index='target_id', columns='sample')
    dd.to_csv(args[OUTPUT], sep='\t', na_rep=0.0)