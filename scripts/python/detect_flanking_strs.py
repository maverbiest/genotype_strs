#!/usr/bin/env python3
import argparse

from match_estrs import get_consensus_unit

ALLOWED_CHECKTYPES = {"simple", "cons", "cons_loose"}

def parse_cla():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--panel", "-p", type=str, required=True, help="Path to STR panel where STRs with overlapping flanking regions will\
            be detected. NOTE: PANEL MUST BE SORTED e.g. 'sort -k1,1V -k2,2n panel > sorted_panel'"
    )
    parser.add_argument(
        "--flanksize", "-s", type=int, default=20, help="Size of flank to consider"
    )
    parser.add_argument(
        '--checktype', "-c", type=str, choices=ALLOWED_CHECKTYPES, default="simple", help="Different options for what \
            type of check will performed on flanking STRs. Options: 'simple' -> no check, report all flank STRs; 'cons\
                ' -> report only flanking STRs that have the same consensus units; cons_loose -> report only flanking \
                    STRs where one consensus unit contains the other"
    )

    return parser.parse_args()

def in_flank(prev_str: dict, cur_str: dict, checktype: str) -> bool:
    if not checktype in ALLOWED_CHECKTYPES:
        raise ValueError(f"checktype must be one of {ALLOWED_CHECKTYPES}")
    if not prev_str['line'][0] == cur_str['line'][0]:
        # previous and current are not not on the same chromosome
        return False
    if not prev_str['flankend'] >= int(cur_str['line'][1]):
        # previous and current are not in eachother's flank
        return False

    # previous and current are in eachother's flank
    if checktype == "simple":
        # no check, always return True
        return True
    elif checktype == "cons":
        # check for matching consensus unit
        cons_prev = get_consensus_unit(prev_str['line'][4])
        cons_cur = get_consensus_unit(cur_str['line'][4])
        if cons_prev == cons_cur:
            return True
        return False
    elif checktype == "cons_loose":
        # check if one consensus unit contains the other
        cons_prev = get_consensus_unit(prev_str['line'][4])
        cons_cur = get_consensus_unit(cur_str['line'][4])
        if cons_prev in cons_cur or cons_cur in cons_prev:
            return True
        return False


def main():
    args = parse_cla()

    with open(args.panel, 'r') as f:
        prev_str = {'line': [], 'flankstart': 0, 'flankend': 0, 'printed': False}
        for line in f:
            line = line.strip().split("\t")
            cur_str = {'line': line, 'flankstart': int(line[1]) - args.flanksize, 'flankend': int(line[2]) + args.flanksize, 'printed': False}
            try:
                if in_flank(prev_str, cur_str, checktype=args.checktype):
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
