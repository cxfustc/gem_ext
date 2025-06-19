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
} comp_t;

typedef struct {
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
	comp_t * pa = (comp_t *) a;
	comp_t * pb = (comp_t *) b;

	if (pa->cnt < pb->cnt)
		return 1;
	else
		return -1;
}

static int
cmp2 (const void * a, const void * b)
{
	comp_t * pa = (comp_t *) a;
	comp_t * pb = (comp_t *) b;

	if (pa->row_min < pb->row_min)
		return -1;
	if (pa->row_min > pb->row_min)
		return 1;

	return  (pa->col_min - pb->col_min);
}

static int
usage (void)
{
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:   gem_ext  part [options] <mask.tif> <gem> <n_component>\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Options: -c  STR  The color used to determine ROI (required)\n");
	fprintf (stderr, "         -b       The color is used in the background (optional)\n");
	fprintf (stderr, "         -s  STR  Sample name (optional, default: prefix of <gem>)\n");
	fprintf (stderr, "         -o  STR  output directory (optional, default: ./\n");
	fprintf (stderr, "         -p  STR  Matched ssDNA image (optional)\n");
	fprintf (stderr, "         -e       keep extra info in gem file\n");
	fprintf (stderr, "         -r  STR  File for arranging parts\n");
	fprintf (stderr, "         -f       Fix coordinate according to expr data\n");
	fprintf (stderr, "\n");
	return 1;
}
int
part_main (int argc, char * argv[])
{
	if (argc < 4)
		return usage ();

	char * buff;
	char * file;
	char * line;
	char * mask_tif;
	char * in_gem;
	char * outdir;
	char * spl_name;
	char * ss_image;
	char * gene;
	char * ch;
	char * ch2;
	char * arr_file;
	int copt;
	int xmin;
	int ymin;
	int num;
	int pidx;
	int i, j, val;
	int x, y;
	int exp_n;
	int n_comps;
	int n_parts;
	int max_cip;
	int is_bg;
	int col_set;
	int col_match;
	int keep_extra;
	int fix_offset;
	int * comp2part;
	int * labl2part;
	FILE * cfg_in;
	FILE * cfg_out;
	Mat img;
	Mat mask;
	Mat label;
	Mat show;
	str_t s;
	gzFile in;
	gzFile * out;
	part_t * p;
	part_t * parts;
	comp_t * comp;
	comp_t * comps;
	gem_size_t * size;
	gem_size_t * gem_sizes;
	color_match col_func;

	buff     = ALLOC_LINE;
	file     = ALLOC_LINE;
	line     = ALLOC_LINE;
	gene     = ALLOC_LINE;
	mask_tif = ALLOC_LINE;
	in_gem   = ALLOC_LINE;
	spl_name = ALLOC_LINE;
	ss_image = ALLOC_LINE;
	arr_file = ALLOC_LINE;

	is_bg = 0;
	col_set = 0;
	col_func = NULL;
	outdir = NULL;
	keep_extra = 0;
	fix_offset = 0;
	while ((copt=getopt(argc,argv,"c:bs:o:p:er:f")) != -1) {
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
		} else if (copt == 'p') {
			strcpy (ss_image, optarg);
		} else if (copt == 'e') {
			keep_extra = 1;
		} else if (copt == 'r') {
			strcpy (arr_file, optarg);
		} else if (copt == 'f') {
			fix_offset = 1;
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
		fprintf (stderr, "parameter <gem> must be set!\n\n");
		return usage ();
	} else
		strcpy (in_gem, argv[optind++]);

	if (optind >= argc) {
		fprintf (stderr, "parameter <n_component> must be set!");
		return usage ();
	} else
		assert ((exp_n = atoi(argv[optind++])) >= 1);

	if (!outdir)
		outdir = get_abs_path (".");

	if (!*spl_name) {
		if ((ch = strrchr(in_gem,'/')) != NULL)
			ch += 1;
		else
			ch = in_gem;
		strcpy (spl_name, ch);
		str_cut_tail (spl_name, ".gz");
		str_cut_tail (spl_name, ".gem");
	}

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

	comps = (comp_t *) ckalloc (n_comps, sizeof(comp_t));
	for (i=0; i<n_comps; ++i) {
		comp = comps + i;
		comp->label = i;
		comp->cnt = 0;
		comp->row_min = INT_MAX;
		comp->col_min = INT_MAX;
		comp->row_max = INT_MIN;
		comp->col_max = INT_MIN;
	}

	for (i=0; i<label.rows; ++i) {
		for (j=0; j<label.cols; ++j) {
			unsigned short v = label.at<unsigned short>(i,j);
			if (mask.data[i*label.cols+j] == 255) {
				comp = comps + v;
				comp->cnt += 1;

				if (i < comp->row_min)
					comp->row_min = i;
				if (j < comp->col_min)
					comp->col_min = j;

				if (i > comp->row_max)
					comp->row_max = i;
				if (j > comp->col_max)
					comp->col_max = j;
			}
		}
	}

	qsort (comps, n_comps, sizeof(comp_t), cmp);
	qsort (comps, exp_n, sizeof(comp_t), cmp2);

	//qsort (parts, n_comps, sizeof(part_t), cmp3);

	if (*arr_file) {
		max_cip = 0;
		n_parts = 0;
		cfg_in = ckopen (arr_file, "r");
		while (fgets(line,4096,cfg_in)) {
			chomp (line);
			if (is_empty_line(line))
				continue;
			assert ((ch = strtok(line," \t")) != NULL);
			num = atoi (ch);
			assert (num > 0);
			if (num > max_cip)
				max_cip = num;

			while ((ch = strtok(NULL," \t")) != NULL) {
				num = atoi (ch);
				assert (num > 0);
				if (num > max_cip)
					max_cip = num;
			}
			++n_parts;
		}
		if (max_cip > exp_n)
			err_mesg ("Components in arrage file is more than expected!");

		comp2part = (int *) ckalloc (exp_n, sizeof(int));

		n_parts = 0;
		rewind (cfg_in);
		while (fgets(line,4096,cfg_in)) {
			chomp (line);
			if (is_empty_line(line))
				continue;
			assert ((ch = strtok(line," \t")) != NULL);
			num = atoi (ch);
			comp2part[num-1] = n_parts;

			while ((ch = strtok(NULL," \t")) != NULL) {
				num = atoi (ch);
				comp2part[num-1] = n_parts;
			}
			++n_parts;
		}
		fclose (cfg_in);
	} else {
		comp2part = (int *) ckalloc (exp_n, sizeof(int));
		n_parts = exp_n;
		for (i=0; i<exp_n; ++i)
			comp2part[i] = i;
	}

	parts = (part_t *) ckalloc (n_parts, sizeof(part_t));
	for (i=0; i<n_parts; ++i) {
		p = parts + i;
		p->row_min = INT_MAX;
		p->row_max = INT_MIN;
		p->col_min = INT_MAX;
		p->col_max = INT_MIN;
	}

	for (i=0; i<exp_n; ++i) {
		comp = comps + i;
		comp->idx = i;

		pidx = comp2part[i];
		p = parts + pidx;
		if (comp->row_min < p->row_min)
			p->row_min = comp->row_min;
		if (comp->row_max > p->row_max)
			p->row_max = comp->row_max;
		if (comp->col_min < p->col_min)
			p->col_min = comp->col_min;
		if (comp->col_max > p->col_max)
			p->col_max = comp->col_max;
	}

	labl2part = (int *) ckalloc (n_comps, sizeof(int));
	for (i=0; i<n_comps; ++i)
		labl2part[i] = -1;
	for (i=0; i<exp_n; ++i) {
		comp = comps + i;
		labl2part[comp->label] = comp2part[i];
	}

	sprintf (file, "%s/%s.mask_with_mark.tif", outdir, spl_name);
	for (i=0; i<exp_n; ++i) {
		comp = comps + i;

		int baseline;
		int row = (comp->row_min + comp->row_max)/2/10;
		int col = (comp->col_min + comp->col_max)/2/10;
		double scale;
		sprintf (buff, "%d", comp2part[i]+1);
		string l (buff);
		Size blk = getTextSize (l, FONT_HERSHEY_SIMPLEX, 1, 2, &baseline);
		scale = (double) 50 / (double)blk.height;
		putText (show, l, Point(col,row), FONT_HERSHEY_SIMPLEX, scale, Scalar(0,0,255), 3);
	}
	imwrite (file, show);

	in = ckgzopen (in_gem, "r");
	while (gzgets(in,line,4096)) {
		if (*line != '#')
			break;
	}
	if (keep_extra) {
		int n_spaces = 0;
		for (ch=line; *ch!='\0'; ++ch)
			if (isspace(*ch) && (++n_spaces)==3)
				break;
		assert (n_spaces == 3);
	}

	out = (gzFile *) ckmalloc (n_parts * sizeof(gzFile));
	gem_sizes = (gem_size_t *) ckmalloc (n_parts * sizeof(gem_size_t));
	for (i=0; i<n_parts; ++i) {
		sprintf (file, "%s/%s.%d.gem.gz", outdir, spl_name, i+1);
		out[i] = ckgzopen (file, "w");
		gzprintf (out[i], "#FileFormat=GEMv0.1\n");
		gzprintf (out[i], "#SortedBy=None\n");
		gzprintf (out[i], "#BinSize=1\n");
		gzprintf (out[i], "#OffsetX=0\n");
		gzprintf (out[i], "#OffsetY=0\n");
		if (!keep_extra)
			gzprintf (out[i], "geneID\tx\ty\tMIDCounts\n");
		else
			gzprintf (out[i], "geneID\tx\ty%s", ch);

		size = gem_sizes + i;
		size->row_min = INT_MAX;
		size->col_min = INT_MAX;
		size->row_max = INT_MIN;
		size->col_max = INT_MIN;
	}

	// extract gem file
	if (fix_offset) {
		xmin = ymin = INT_MAX;
		while (gzgets(in,line,4096)) {
			sscanf (line, "%*s %d %d", &x, &y);
			if (x < xmin) xmin = x;
			if (y < ymin) ymin = y;
		}
		gzclose (in);
	} else {
		xmin = ymin = 0;
	}
	
	in = ckgzopen (in_gem, "r");
	while (gzgets(in,line,4096)) {
		if (*line != '#')
			break;
	}

	fprintf (stderr, "Offset X: %d\n", xmin);
	fprintf (stderr, "Offset Y: %d\n", ymin);

	while (gzgets(in,line,4096)) {
		sscanf (line, "%s %d %d %d", gene, &x, &y, &val);
		x = x - xmin;
		y = y - ymin;

		if (x<0 || x>label.cols)
			continue;
		if (y<0 || y>label.rows)
			continue;

		unsigned short v = label.at<unsigned short>(y,x);
		assert (v>=0 && v<n_comps);
		if ((pidx=labl2part[v]) < 0)
			continue;
		p = parts + pidx;

		if (!keep_extra) {
			gzprintf (out[pidx], "%s\t%d\t%d\t%d\n", gene, x, y, val);
		} else {
			int n_spaces = 0;
			for (ch=line; *ch!='\0'; ++ch)
				if (isspace(*ch) && (++n_spaces)==3)
					break;
			assert (n_spaces == 3);
			gzprintf (out[pidx], "%s\t%d\t%d%s", gene, x, y, ch);
		}

		size = gem_sizes + pidx;
		if (x < size->col_min) size->col_min = x;
		if (x > size->col_max) size->col_max = x;
		if (y < size->row_min) size->row_min = y;
		if (y > size->row_max) size->row_max = y;
	}
	gzclose (in);

	for (i=0; i<n_parts; ++i)
		gzclose (out[i]);

	//qsort (parts, n_comps, sizeof(part_t), cmp);
	//qsort (parts, exp_n, sizeof(part_t), cmp2);

	sprintf (file, "%s/%s.info.txt", outdir, spl_name);
	cfg_out = ckopen (file, "w");

	for (i=0; i<n_parts; ++i) {
		p = parts + i;
		size = gem_sizes + i;

		fprintf (cfg_out, "%s/%s.%d.gem.gz\t0\t0\t0\t", outdir, spl_name, i+1);
		if (*ss_image)
			fprintf (cfg_out, "%s/%s.%d.part.tif\t", outdir, spl_name, i+1);
		else
			fprintf (cfg_out, "*\t");
		fprintf (cfg_out, "%d\t%d\t%d\t%d\n",
				size->col_min - p->col_min, size->row_min - p->row_min,
				size->col_max - p->col_max, size->row_max - p->row_max);
	}
	fclose (cfg_out);

	if (!*ss_image)
		return 0;

	Mat ss;
	vector<Mat> imgs;

	for (i=0; i<n_parts; ++i) {
		p = parts + i;
		Mat m = Mat::zeros (Size(p->col_max-p->col_min+1, p->row_max-p->row_min+1), CV_8UC3);
		imgs.push_back (m);
	}

	//qsort (parts, n_comps, sizeof(part_t), cmp3);

	/*
	for (i=0; i<n_comps; ++i) {
		p = parts + i;
		if (p->idx < 0)
			continue;

		printf ("%d-%d\t%d-%d\n",
				p->col_max-p->col_min, p->row_max-p->row_min,
				imgs[p->idx].cols, imgs[p->idx].rows);
	}
	*/

	xcv_imread (ss, ss_image);
	assert (ss.rows == label.rows);
	assert (ss.cols == label.cols);

	for (i=0; i<ss.rows; ++i) {
		for (j=0; j<ss.cols; ++j) {
			unsigned short v = label.at<unsigned short>(i,j);
			assert (v>=0 && v<n_comps);
			pidx = labl2part[v];
			if (pidx < 0)
				continue;
			p = parts + pidx;
			imgs[pidx].at<Vec3b>(i-p->row_min,j-p->col_min) = ss.at<Vec3b>(i,j);
		}
	}

	for (i=0; i<n_parts; ++i) {
		sprintf (file, "%s/%s.%d.part.tif", outdir, spl_name, i+1);
		imwrite (file, imgs[i]);
	}

	return 0;
}
