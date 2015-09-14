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
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "sam.h"
#include "khash.h"
#include "samtools.h"
#include "sam_opts.h"
KHASH_SET_INIT_STR(rg)

typedef khash_t(rg) *rghash_t;

// This structure contains the settings for a samview run
typedef struct samview_settings {
    rghash_t rghash;
    int min_mapQ;
    int flag_on;
    int flag_off;
    int min_qlen;
    int remove_B;
    uint32_t subsam_seed;
    double subsam_frac;
    char* library;
    void* bed;
    size_t remove_aux_len;
    char** remove_aux;
} samview_settings_t;


// TODO Add declarations of these to a viable htslib or samtools header
//extern const char *bam_get_library(bam_hdr_t *header, const bam1_t *b);
//extern int bam_remove_B(bam1_t *b);
//extern char *samfaipath(const char *fn_ref);
//void *bed_read(const char *fn);
//void bed_destroy(void *_h);
//int bed_overlap(const void *_h, const char *chr, int beg, int end);

void print_error(const char *subcommand, const char *format, ...)
{
   return;
}

// Returns 0 to indicate read should be output 1 otherwise
static int process_aln(const bam_hdr_t *h, bam1_t *b, samview_settings_t* settings)
{
//    if (settings->remove_B) bam_remove_B(b);
//    if (settings->min_qlen > 0) {
//        int k, qlen = 0;
//        uint32_t *cigar = bam_get_cigar(b);
//        for (k = 0; k < b->core.n_cigar; ++k)
//            if ((bam_cigar_type(bam_cigar_op(cigar[k]))&1) || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
//                qlen += bam_cigar_oplen(cigar[k]);
//        if (qlen < settings->min_qlen) return 1;
//    }
//    if (b->core.qual < settings->min_mapQ || ((b->core.flag & settings->flag_on) != settings->flag_on) || (b->core.flag & settings->flag_off))
//        return 1;
//    if (settings->bed && (b->core.tid < 0 || !bed_overlap(settings->bed, h->target_name[b->core.tid], b->core.pos, bam_endpos(b))))
//        return 1;
//    if (settings->subsam_frac > 0.) {
//        uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ settings->subsam_seed);
//        if ((double)(k&0xffffff) / 0x1000000 >= settings->subsam_frac) return 1;
//    }
//    if (settings->rghash) {
//        uint8_t *s = bam_aux_get(b, "RG");
//        if (s) {
//            khint_t k = kh_get(rg, settings->rghash, (char*)(s + 1));
//            if (k == kh_end(settings->rghash)) return 1;
//        }
//    }
//    if (settings->library) {
//        const char *p = bam_get_library((bam_hdr_t*)h, b);
//        if (p && strcmp(p, settings->library) != 0) return 1;
//    }
//    if (settings->remove_aux_len) {
//        size_t i;
//        for (i = 0; i < settings->remove_aux_len; ++i) {
//            uint8_t *s = bam_aux_get(b, settings->remove_aux[i]);
//            if (s) {
//                bam_aux_del(b, s);
//            }
//        }
//    }
    return 0;
}

static inline int check_sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b, const char *fname, int *retp)
{
    int r = sam_write1(fp, h, b);
    if (r >= 0) return r;

    if (fname) print_error("view", "writing to \"%s\" failed", fname);
    else print_error("view", "writing to standard output failed");

    *retp = EXIT_FAILURE;
    return r;
}

static void check_sam_close(const char *subcmd, samFile *fp, const char *fname, const char *null_fname, int *retp)
{
    int r = sam_close(fp);
    if (r >= 0) return;

    // TODO Need error infrastructure so we can print a message instead of r
    if (fname) print_error(subcmd, "error closing \"%s\": %d", fname, r);
    else print_error(subcmd, "error closing %s: %d", null_fname, r);

    *retp = EXIT_FAILURE;
}

int main_samview(int argc, char *argv[])
{
    int c, is_header = 0, is_header_only = 0, ret = 0, compress_level = -1, is_count = 0;
    int is_long_help = 0, n_threads = 0;
    samFile *in = 0, *out = 0, *un_out=0;
    bam_hdr_t *header = NULL;
    char out_mode[5] = "w\0\0\0\0", out_un_mode[5], *out_format = "";
    char *fn_out = 0, *fn_list = 0, *q, *fn_un_out = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    samview_settings_t settings = {
        .rghash = NULL,
        .min_mapQ = 0,
        .flag_on = 0,
        .flag_off = 0,
        .min_qlen = 0,
        .remove_B = 0,
        .subsam_seed = 0,
        .subsam_frac = -1.,
        .library = NULL,
        .bed = NULL,
    };

    char *fn_in = argv[1];

    if ((in = sam_open_format(fn_in, "r", &ga.in)) == 0) {
        print_error("view", "failed to open \"%s\" for reading", fn_in);
        ret = 1;
        goto view_end;
    }

    if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", fn_in);
        ret = 1;
        goto view_end;
    }
    if (!is_count) {
        if ((out = sam_open_format(fn_out? fn_out : "-", out_mode, &ga.out)) == 0) {
            print_error("view", "failed to open \"%s\" for writing", fn_out? fn_out : "standard output");
            ret = 1;
            goto view_end;
        }
    }
        bam1_t *b = bam_init1();
        int r;
        while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
            if (!process_aln(header, b, &settings)) {
                if (!is_count) { if (check_sam_write1(out, header, b, fn_out, &ret) < 0) break; }
            } else {
                if (un_out) { if (check_sam_write1(un_out, header, b, fn_un_out, &ret) < 0) break; }
            }
        }
        if (r < -1) {
            fprintf(stderr, "[main_samview] truncated file.\n");
            ret = 1;
        }
        bam_destroy1(b);

view_end:
    // close files, free and return
    if (in) check_sam_close("view", in, fn_in, "standard input", &ret);
    if (out) check_sam_close("view", out, fn_out, "standard output", &ret);

    free(fn_list); free(fn_out); free(settings.library);  free(fn_un_out);
    sam_global_args_free(&ga);
    if ( header ) bam_hdr_destroy(header);
    return ret;
}

int main (int argc, char *argv[]) {
   return main_samview(argc, argv);
}
