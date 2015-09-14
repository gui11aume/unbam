/*  hts.c -- format-neutral I/O, indexing, and iterator API functions.

    Copyright (C) 2008, 2009, 2012-2015 Genome Research Ltd.
    Copyright (C) 2012, 2013 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <zlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include "bgzf.h"
#include "hts.h"
#include "hfile.h"

// Needed.
#include "khash.h"
KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

int hts_verbose = 3;

const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

const int seq_nt16_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

// Decompress up to ten or so bytes by peeking at the file, which must be
// positioned at the start of a GZIP block.
static size_t decompress_peek(hFILE *fp, unsigned char *dest, size_t destsize)
{
    // Typically at most a couple of hundred bytes of input are required
    // to get a few bytes of output from inflate(), so hopefully this buffer
    // size suffices in general.
    unsigned char buffer[512];
    z_stream zs;
    ssize_t npeek = hpeek(fp, buffer, sizeof buffer);

    if (npeek < 0) return 0;

    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = buffer;
    zs.avail_in = npeek;
    zs.next_out = dest;
    zs.avail_out = destsize;
    if (inflateInit2(&zs, 31) != Z_OK) return 0;

    while (zs.total_out < destsize)
        if (inflate(&zs, Z_SYNC_FLUSH) != Z_OK) break;

    destsize = zs.total_out;
    inflateEnd(&zs);

    return destsize;
}


htsFile *hts_hopen(struct hFILE *hfile, const char *fn)
{

    htsFile *fp = (htsFile*) calloc(1, sizeof(htsFile));
    if (fp == NULL) return NULL;


    // Make sure that file is in bam format.
    // Also set fp->format to the right values.
    htsFormat *fmt = &fp->format;
    unsigned char s[21];
    ssize_t len = hpeek(hfile, s, 18);
    if (len < 0) return NULL;

    if (len >= 2 && s[0] == 0x1f && s[1] == 0x8b) {
        // The stream is either gzip-compressed or BGZF-compressed.
        // Determine which, and decompress the first few bytes.
        fmt->compression = (len >= 18 && (s[3] & 4) &&
                            memcmp(&s[12], "BC\2\0", 4) == 0)? bgzf : gzip;
        len = decompress_peek(hfile, s, sizeof s);
    }
    // Not a bam file.
    else return NULL;

    fmt->compression_level = -1;
    fmt->specific = NULL;

    if (len >= 4 && s[3] <= '\4' && memcmp(s, "BAM\1", 4) == 0) {
       fmt->category = sequence_data;
       fmt->format = bam;
       // TODO Decompress enough to pick version from @HD-VN header
       fmt->version.major = 1, fmt->version.minor = -1;
    }
    // Not a bam file.
    else return NULL;

    fp->fn = strdup(fn);
    fp->is_be = ed_is_big();
    fp->is_bin = 1;

    fp->fp.bgzf = bgzf_hopen(hfile, "r");

    if (fp->fp.bgzf == NULL) return NULL;

    return fp;

}

int hts_close(htsFile *fp)
{
    int ret = bgzf_close(fp->fp.bgzf);
    free(fp->fn);
    free(fp->fn_aux);
    free(fp->line.s);
    free(fp);
    return ret;
}
