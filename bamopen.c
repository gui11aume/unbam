/*  sam_view.c -- SAM<->BAM<->CRAM conversion.

    Copyright (C) 2009-2014 Genome Research Ltd.
    Portions copyright (C) 2009, 2011, 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "sam.h"
#include "hfile.h"
#include "bgzf.h"

#ifndef bam_get_seq
#define bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
#endif

int main_samview(int argc, char *argv[])
{

    BGZF *in = bgzf_open(argv[1], "r");
    bam_hdr_t *header = NULL;
    header = bam_hdr_read(in);

    bam1_t *b = bam_init1();
    int r = 0;
    while ((r = bam_read1(in, b)) >= 0) {
       if (r >= 0) { 
          if (b->core.tid  >= header->n_targets || b->core.tid  < -1 ||
             b->core.mtid >= header->n_targets || b->core.mtid < -1)
             return 1;
       }
       // Process 'b' here.
       printf("%s\n", (char *) bam_get_seq(b));
    }

    bam_destroy1(b);
    bgzf_close(in);
    bam_hdr_destroy(header);
    return 0;
}

int main (int argc, char *argv[]) {
   return main_samview(argc, argv);
}
