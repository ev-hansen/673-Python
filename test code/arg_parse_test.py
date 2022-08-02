import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--test', type=float, help='test', required=False, 
                    default=1.0)
parser.add_argument('-t2', '--test2', type=str)
args = parser.parse_args()
args_dict = vars(args)
print(args_dict)
print(args_dict['test'])
