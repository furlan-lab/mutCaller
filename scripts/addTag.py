#!/usr/bin/env python

import argparse
from sys import stdout

def parse_read(ss, separator="|BARCODE=", cblen=16, umilen=10):
  split=ss.qname.split(separator)
  ss.qname=f'{split[0]}:{split[1]}'
  cb=str(split[int(1)][int(0):int(cblen)])+str("-1")
  umi=split[int(1)][int(cblen):int(umilen)+int(cblen)]
  ss.tags.update({'CB':cb, 'CR':cb, 'UB':umi, 'UR':umi})
  return(ss)

def iterate(args):
  from simplesam import Reader, Writer
  with Reader(args.bam) as bam:
    bam.header.get('@HD')['VN:1.4']=['SO:coordinate']
    stdout_sam = Writer(stdout, bam.header)
    for read in bam:
      stdout_sam.write(parse_read(read, args.sep, args.cb_len, args.umi_len))
    args.bam.close()

def main():
  parser = argparse.ArgumentParser(prog='addTags', description="parse BAM sequence name for barcode and add as bam tags")
  parser.add_argument('bam', type=argparse.FileType('r'), help=" BAM file ")
  parser.add_argument('-u', '--umi_len', default=10, help="length of umi")
  parser.add_argument('-c', '--cb_len', default=18, help="length of cb")
  parser.add_argument('-s', '--sep', default="|BARCODE=", help="string separator")
  parser.set_defaults(func=iterate)
  args = parser.parse_args()
  args.func(args)

if __name__ == "__main__":
  main()

"""
from sys import stdout
#debugging stuff
inbam=open("Aligned.out.sampled.bam", 'r')
from simplesam import Reader, Writer

with Reader(inbam) as bam:
  bam.header.get('@HD')['VN:1.4']=['SO:coordinate']
  stdout_sam = Writer(stdout, bam.header)
  read = next(bam)
  stdout_sam.write(parse_read(read))
  # for read in bam:
  #   parse_
  bam.close()
"""
