/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2022-05-19 17:42:47
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>

#include "image.h"
#include "utils.h"

using namespace cv;

typedef int (*color_match) (Vec3b &);

extern Vec3b xcol_roi;
extern int is_red (Vec3b & col);
extern int is_green (Vec3b & col);
extern int is_blue (Vec3b & col);
extern int is_col (Vec3b & col);

static int
usage (void)
{
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:  ext_gem bin_extr <mask.tif> <gem> <bin.size>\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Options: -c  STR  The color used to determine ROI (required)\n");
	fprintf (stderr, "         -b       The color is used in the background (optional)\n");
	fprintf (stderr, "         -s  STR  Sample name (optional, default: prefix of <gem>)\n");
	fprintf (stderr, "         -o  STR  Output directory (optional, default: ./\n");
	fprintf (stderr, "\n");

	return 1;
}

int
bin_extr_main (int argc, char * argv[])
{
	if (argc < 3)
		return usage ();

	char * file;
	char * line;
	char * gene;
	char * mask;
	char * in_gem;
	char * outdir;
	char * spl_name;
	char * ch;
	int x, y, v;
	int x0, y0;
	int copt;
	int is_bg;
	int col_set;
	int bin_size;
	int xmin, ymin;
	int rad;
	int col_match;
	unsigned int num;
	Mat img;
	gzFile in;
	gzFile out;
	color_match col_func;

	file     = ALLOC_LINE;
	line     = ALLOC_LINE;
	gene     = ALLOC_LINE;
	mask     = ALLOC_LINE;
	in_gem   = ALLOC_LINE;
	outdir   = ALLOC_LINE;
	spl_name = ALLOC_LINE;

	is_bg = 0;
	col_set = 0;
	col_func = NULL;
	while ((copt = getopt(argc,argv,"c:b")) != -1) {
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
			strcpy (outdir, optarg);
		}
	}

	if (!col_set || !col_func) {
		fprintf (stderr, "parameter '-c' for color must be set!\n\n");
		return usage ();
	}

	if (optind >= argc) {
		fprintf (stderr, "parameter <mask.tif> must be set!\n\n");
		return usage ();
	} else
		strcpy (mask, argv[optind++]);

	if (optind >= argc) {
		fprintf (stderr, "parameter <gem> must be set!\n\n");
		return usage ();
	} else
		strcpy (in_gem, argv[optind++]);

	if (optind >= argc) {
		fprintf (stderr, "parameter <bin.size> must be set!\n\n");
		return usage ();
	} else
		bin_size = atoi (argv[optind++]);

	if (!*outdir)
		strcpy (outdir, ".");

	if (!*spl_name) {
		if ((ch = strrchr(in_gem,'/')) != NULL)
			ch += 1;
		else
			ch = in_gem;
		strcpy (spl_name, ch);
		str_cut_tail (spl_name, ".gz");
		str_cut_tail (spl_name, ".gem");
	}

	in = ckgzopen (argv[2], "r");
	while (gzgets(in,line,4096))
		if (*line != '#')
			break;

	xmin = ymin = INT_MAX;
	while (gzgets(in,line,4096)) {
		sscanf (line, "%*s %d %d", &x, &y);
		if (x < xmin) xmin = x;
		if (y < ymin) ymin = y;
	}
	gzclose (in);

	fprintf (stderr, "Input X offset: %d\n", xmin);
	fprintf (stderr, "Input Y offset: %d\n", ymin);

	xcv_imread (img, mask);
	sprintf (file, "%s/%s.register.gem.gz", outdir, spl_name);
	out = ckgzopen (file, "w");
	gzprintf (out, "#FileFormat=GEMv0.1\n");
	gzprintf (out, "#SortedBy=None\n");
	gzprintf (out, "#BinSize=1\n");
	gzprintf (out, "#OffsetX=0\n");
	gzprintf (out, "#OffsetY=0\n");
	gzprintf (out, "geneID\tx\ty\tMIDCounts\n");
	rad = 5;
	while (gzgets(in,line,4096)) {
		sscanf (line, "%s %d %d %d", gene, &x0, &y0, &v);
		x0 = x0 - xmin;
		y0 = y0 - ymin;

		x = x0 / bin_size * 2 * rad;
		y = y0 / bin_size * 2 * rad;

		if (x<0 || x>=img.cols)
			continue;
		if (y<0 || y>=img.rows)
			continue;

		Vec3b col = img.at<Vec3b>(y,x);

		col_match = col_func (col);
		if (!is_bg && col_match) {
			gzprintf (out, "%s\t%d\t%d\t%d\n", gene, x0, y0, v);
		} else if (is_bg && !col_match) {
			gzprintf (out, "%s\t%d\t%d\t%d\n", gene, x0, y0, v);
		}
	}
	gzclose (out);

	return 0;
}
