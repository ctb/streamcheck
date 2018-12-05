#! /usr/bin/env python
import sourmash
import screed
import sys

REPORT_BP = 10000
CUTOFF = 10


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('reference')
    parser.add_argument('reads', nargs='+')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    args = parser.parse_args()

    ref_mh = sourmash.MinHash(0, args.ksize, scaled=1000)
    reads_mh = sourmash.MinHash(0, args.ksize, scaled=1000, track_abundance=True)

    print('reading reference file: {}'.format(args.reference))
    n = 0
    bp = 0
    for record in screed.open(args.reference):
        ref_mh.add_sequence(record.sequence)
        n += 1
        bp += len(record.sequence)

    n = 0
    bp = 0
    watermark = REPORT_BP
    for filename in args.reads:
        for record in screed.open(filename):
            reads_mh.add_sequence(record.sequence)
            n += 1
            bp += len(record.sequence)

            if bp >= watermark:
                highabund_mh = ref_mh.copy_and_clear()
                for hashval, abund in reads_mh.get_mins(with_abundance=True).items():
                    if abund >= CUTOFF:
                        highabund_mh.add_hash(hashval)
                
                print('{} ref_in_r={:.3f} ref_highr_sim={:.3f} ref_in_highr={:.3f} highr_in_ref={:.3f}'.format(n, ref_mh.contained_by(reads_mh), ref_mh.jaccard(highabund_mh), ref_mh.contained_by(highabund_mh), highabund_mh.contained_by(ref_mh)))
                watermark += REPORT_BP

    return 0
    

if __name__ == '__main__':
    sys.exit(main())
