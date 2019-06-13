import argparse 
import os
import shutil
import sys

BASE_DIR = 'base_dir'

def parse_input():
    '''
    Parses the commandline input, returns a dict
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('file_list', metavar='File', nargs='+')
    parser.add_argument('-b', required=True, dest=BASE_DIR)
    args = parser.parse_args()
    return vars(args)


if __name__ == '__main__':
    
    args = parse_input()
    base_dir = args[BASE_DIR]
    file_list = args['file_list']

    # get a set of the unique contrast names
    contrast_list = [os.path.basename(x).split('.')[0] for x in file_list]
    contrast_set = set([x for x in contrast_list])

    # make a directory for each contrast and create a mapping
    dir_mapping = {}
    for x in contrast_set:
        dir_mapping[x] = os.path.join(base_dir, x)
        try:
            os.mkdir(os.path.join(base_dir, x))
        except OSError as ex:
            if ex.errno == 17:
                pass
            else:
                sys.exit(1)

    # move the files into the appropriate dirs:
    for contrast, f in zip(contrast_list, file_list):
        shutil.move(f, dir_mapping[contrast])
