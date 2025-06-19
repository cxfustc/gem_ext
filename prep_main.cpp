/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2022-03-24 11:20:35
  *Edit History: 
***********************************************************/

#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "mp.h"
#include "image.h"
#include "utils.h"
#include "str_hash.h"

using namespace cv;
using namespace std;

typedef int (*color_match) (Vec3b &);

extern Vec3b xcol_roi;
extern int is_red (Vec3b & col);
extern int is_green (Vec3b & col);
extern int is_blue (Vec3b & col);
extern int is_col (Vec3b & col);

typedef struct {
	int label;
	int row_min;
	int col_min;
	int row_max;
	int col_max;
	int idx;
	long cnt;
} part_t;

typedef struct {
	int row_min;
	int row_max;
	int col_min;
	int col_max;
} gem_size_t;

static int
cmp (const void * a, const void * b)
{
	part_t * pa = (part_t *) a;
	part_t * pb = (part_t *) b;

	if (pa->cnt < pb->cnt)
		return 1;
	else
		return -1;
}

static int
cmp2 (const void * a, const void * b)
{
	part_t * pa = (part_t *) a;
	part_t * pb = (part_t *) b;

	if (pa->row_min < pb->row_min)
		return -1;
	if (pa->row_min > pb->row_min)
		return 1;

	return  (pa->col_min - pb->col_min);
}

static int
cmp3 (const void * a, const void * b)
{
	part_t * pa = (part_t *) a;
	part_t * pb = (part_t *) b;

	return (pa->label - pb->label);
}

static int
usage (void)
{
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:   gem_ext  prep [options] <mask.tif> <n_component>\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Options: -c  STR  The color used to determine ROI (required)\n");
	fprintf (stderr, "         -b       The color is used in the background (optional)\n");
	fprintf (stderr, "         -s  STR  Sample name (optional, default: prefix of <gem>)\n");
	fprintf (stderr, "         -o  STR  output directory (optional, default: ./\n");
	fprintf (stderr, "\n");
	return 1;
}

int
prep_main (int argc, char * argv[])
{
	if (argc < 3)
		return usage ();

	char * buff;
	char * file;
	char * line;
	char * mask_tif;
	char * outdir;
	char * spl_name;
	char * ss_image;
	char * gene;
	char * ch;
	char * ch2;
	int copt;
	int xmin;
	int ymin;
	int num;
	int i, j, val;
	int x, y;
	int exp_n;
	int n_comps;
	int is_bg;
	int col_set;
	int col_match;
	int keep_extra;
	FILE * cfg_out;
	Mat img;
	Mat mask;
	Mat label;
	Mat show;
	str_t s;
	gzFile in;
	gzFile * out;
	part_t tmp;
	part_t * p;
	part_t * parts;
	gem_size_t * size;
	gem_size_t * gem_sizes;
	color_match col_func;

	buff     = ALLOC_LINE;
	file     = ALLOC_LINE;
	line     = ALLOC_LINE;
	gene     = ALLOC_LINE;
	mask_tif = ALLOC_LINE;
	spl_name = ALLOC_LINE;
	ss_image = ALLOC_LINE;

	is_bg = 0;
	col_set = 0;
	col_func = NULL;
	outdir = NULL;
	keep_extra = 0;
	while ((copt=getopt(argc,argv,"c:bs:o:")) != -1) {
		if (copt == 'c') {
			if (strcmp(optarg,"r")==0 || strcmp(optarg,"red")==0)
				col_func = is_red;
			else if (strcmp(optarg,"g")==0 || strcmp(optarg,"green")==0)
				col_func = is_green;
			else if (strcmp(optarg,"b")==0 || strcmp(optarg,"blue")==0)
				col_func = is_blue;
			else {
				ch = optarg;
				if (*ch == '#')
					ch += 1;
				if (strlen(ch) != 6) {
					fprintf (stderr, "Invalid color: '%s'\n\n", optarg);
					return usage ();
				}
				sscanf (ch, "%x", &num);
				int r = (num>>16) & 0xff;
				int g = (num>>8)  & 0xff;
				int b =  num      & 0xff;
				xcol_roi = Vec3b (b, g, r);
				col_func = is_col;
			}
			col_set = 1;
		} else if (copt == 'b') {
			is_bg = 1;
		} else if (copt == 's') {
			strcpy (spl_name, optarg);
		} else if (copt == 'o') {
			outdir = get_abs_path (optarg);
		}
	}

	if (!col_set) {
		fprintf (stderr, "parameter '-c' for color must be set!\n\n");
		return usage ();
	}

	if (optind >= argc) {
		fprintf (stderr, "parameter <mask.tif> must be set!\n\n");
		return usage ();
	} else
		strcpy (mask_tif, argv[optind++]);

	if (optind >= argc) {
		fprintf (stderr, "parameter <n_component> must be set!");
		return usage ();
	} else
		assert ((exp_n = atoi(argv[optind++])) >= 1);

	if (!outdir)
		outdir = get_abs_path (".");

	if (!*spl_name)
		strcpy (spl_name, "sample");

	xcv_imread (img, mask_tif);
	mask = Mat::zeros (Size(img.cols, img.rows), CV_8UC1);
	show = Mat::zeros (Size(img.cols/10+1, img.rows/10+1), CV_8UC3);
	for (i=0; i<img.rows; ++i) {
		for (j=0; j<img.cols; ++j) {
			Vec3b c = img.at<Vec3b>(i,j);

			col_match = col_func (c);
			if (!is_bg && col_match) {
				mask.data[i*img.cols+j] = 255;
				if (i%10==0 && j%10==0)
					show.at<Vec3b>(i/10,j/10) = Vec3b(255,255,255);
			} else if (is_bg && !col_match) {
				mask.data[i*img.cols+j] = 255;
				if (i%10==0 && j%10==0)
					show.at<Vec3b>(i/10,j/10) = Vec3b(255,255,255);
			}
		}
	}

	n_comps = connectedComponents (mask, label, 8, CV_16U);
	if (n_comps < exp_n)
		err_mesg ("Total components are less than expected! Please Check mask tif!");

	parts = (part_t *) ckalloc (n_comps, sizeof(part_t));
	for (i=0; i<n_comps; ++i) {
		p = parts + i;
		p->label = i;
		p->cnt = 0;
		p->row_min = INT_MAX;
		p->col_min = INT_MAX;
		p->row_max = INT_MIN;
		p->col_max = INT_MIN;
	}

	for (i=0; i<label.rows; ++i) {
		for (j=0; j<label.cols; ++j) {
			unsigned short v = label.at<unsigned short>(i,j);
			if (mask.data[i*label.cols+j] == 255) {
				p = parts + v;
				p->cnt += 1;

				if (i < p->row_min)
					p->row_min = i;
				if (j < p->col_min)
					p->col_min = j;

				if (i > p->row_max)
					p->row_max = i;
				if (j > p->col_max)
					p->col_max = j;
			}
		}
	}

	qsort (parts, n_comps, sizeof(part_t), cmp);
	for (i=exp_n; i<n_comps; ++i)
		parts[i].idx = -1;

	qsort (parts, exp_n, sizeof(part_t), cmp2);
	for (i=0; i<exp_n; ++i)
		parts[i].idx = i;

	qsort (parts, n_comps, sizeof(part_t), cmp3);

	sprintf (file, "%s/%s.mask_with_mark.tif", outdir, spl_name);
	for (i=0; i<n_comps; ++i) {
		p = parts + i;
		fprintf (stderr, "%d\t%ld\t%d-%d\t%d\n", p->label, p->cnt, p->row_min, p->col_min, p->idx);

		if (p->idx < 0)
			continue;

		int baseline;
		int row = (p->row_min + p->row_max)/2/10;
		int col = (p->col_min + p->col_max)/2/10;
		double scale;
		sprintf (buff, "%d", p->idx+1);
		string l (buff);
		Size blk = getTextSize (l, FONT_HERSHEY_SIMPLEX, 1, 2, &baseline);
		scale = (double) 50 / (double)blk.height;
		putText (show, l, Point(col,row), FONT_HERSHEY_SIMPLEX, scale, Scalar(0,0,255), 3);
	}
	imwrite (file, show);

	return 0;
}
