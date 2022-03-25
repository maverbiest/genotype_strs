#!/usr/bin/env python3
import argparse

def parse_cla():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--panel", "-p", type=str, required=True, help="Path to STR panel where STRs with overlapping flanking regions will\
            be detected. NOTE: PANEL MUST BE SORTED e.g. 'sort -k1,1V -k2,2n panel > sorted_panel'"
    )
    parser.add_argument(
        "--flanksize", "-s", type=int, default=20, help="Size of flank to consider"
    )

    return parser.parse_args()


def main():
    args = parse_cla()

    with open(args.panel, 'r') as f:
        prev_str = {'line': [], 'flankstart': 0, 'flankend': 0, 'printed': False}
        for line in f:
            line = line.strip().split("\t")
            cur_str = {'line': line, 'flankstart': int(line[1]) - args.flanksize, 'flankend': int(line[2]) + args.flanksize, 'printed': False}
            try:
                if prev_str['line'][0] == cur_str['line'][0] and prev_str['flankend'] >= int(cur_str['line'][1]):
                    if not prev_str['printed']:
                        print("\t".join(prev_str['line']))                    

                    print("\t".join(cur_str['line']))
                    cur_str['printed'] = True                    
            except IndexError:
                # should only happen on the first row when no 'prev_str' with values exists
                pass

            prev_str = cur_str

if __name__ == "__main__":
    main()
