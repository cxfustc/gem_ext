/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2022-03-14 17:08:13
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"
#include "utils.h"

using namespace cv;
using namespace std;

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
	fprintf (stderr, "Usage:   ext_gem extr <mask.tif> <gem>\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Options: -c  STR  The color used to determine ROI (required)\n");
	fprintf (stderr, "         -b       The color is used in the background (optional)\n");
	fprintf (stderr, "         -s  STR  Sample name (optional, default: prefix of <gem>)\n");
	fprintf (stderr, "         -o  STR  Output directory (optional, default: ./\n");
	fprintf (stderr, "         -p  STR  Matched ssDNA image (optional)\n");
	fprintf (stderr, "         -e       Keep extra info in gem file\n");
	fprintf (stderr, "         -f       Fix coordinate according to expr data\n");
	fprintf (stderr, "\n");
	return 1;
}

int
extract_main (int argc, char * argv[])
{
	if (argc < 3)
		return usage ();

	char * file;
	char * line;
	char * mask;
	char * in_gem;
	char * outdir;
	char * spl_name;
	char * ss_image;
	char * gene;
	char * ch;
	char * ch2;
	int copt;
	int x, y, v;
	int xmin;
	int ymin;
	int rmin, rmax;
	int cmin, cmax;
	int is_bg;
	int col_set;
	int col_match;
	int fix_offset;
	int min[2];
	int max[2];
	int keep_extra;
	unsigned int num;
	Mat img;
	color_match col_func;
	gzFile in;
	gzFile out;

	file     = ALLOC_LINE;
	gene     = ALLOC_LINE;
	line     = ALLOC_LINE;
	mask     = ALLOC_LINE;
	in_gem   = ALLOC_LINE;
	outdir   = ALLOC_LINE;
	spl_name = ALLOC_LINE;
	ss_image = ALLOC_LINE;

	is_bg = 0;
	col_set = 0;
	keep_extra = 0;
	fix_offset = 0;
	col_func = NULL;
	while ((copt=getopt(argc,argv,"c:bs:o:p:ef")) != -1) {
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
		} else if (copt == 'p') {
			strcpy (ss_image, optarg);
		} else if (copt == 'e') {
			keep_extra = 1;
		} else if (copt == 'f') {
			fix_offset = 1;
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

	in = ckgzopen (in_gem, "r");

	xcv_imread (img, mask);

	// extract gem file
	while (gzgets(in,line,4096)) {
		if (*line != '#')
			break;
	}

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
	if (keep_extra) {
		int n_spaces = 0;
		for (ch=line; *ch!='\0'; ++ch)
      if (isspace(*ch) && (++n_spaces)==3)
        break;
    assert (n_spaces == 3);
	}

	fprintf (stderr, "Input Offset X: %d\n", xmin);
	fprintf (stderr, "Input Offset Y: %d\n", ymin);

	sprintf (file, "%s/%s.register.gem.gz", outdir, spl_name);
	out = ckgzopen (file, "w");
	gzprintf (out, "#FileFormat=GEMv0.1\n");
	gzprintf (out, "#SortedBy=None\n");
	gzprintf (out, "#BinSize=1\n");
	gzprintf (out, "#OffsetX=0\n");
	gzprintf (out, "#OffsetY=0\n");
	if (!keep_extra)
		gzprintf (out, "geneID\tx\ty\tMIDCounts\n");
	else
		gzprintf (out, "geneID\tx\ty%s", ch);

	min[0] = min[1] = INT_MAX;
	max[0] = max[1] = INT_MIN;
	while (gzgets(in,line,4096)) {
		sscanf (line, "%s %d %d %d", gene, &x, &y, &v);
		int xi = x - xmin;
		int yi = y - ymin;

		if (xi<0 || xi>=img.cols)
			continue;
		if (yi<0 || yi>=img.rows)
			continue;

		Vec3b col = img.at<Vec3b>(yi,xi);

		col_match = col_func (col);

		/*
		if (!is_bg && col_match) {
			gzprintf (out, "%s\t%d\t%d\t%d\n", gene, x, y, v);
			if (x < min[0]) min[0] = x;
			if (x > max[0]) max[0] = x;
			if (y < min[1]) min[1] = y;
			if (y > max[1]) max[1] = y;
		} else if (is_bg && !col_match) {
			gzprintf (out, "%s\t%d\t%d\t%d\n", gene, x, y, v);
			if (x < min[0]) min[0] = x;
			if (x > max[0]) max[0] = x;
			if (y < min[1]) min[1] = y;
			if (y > max[1]) max[1] = y;
		}
		*/

		if ((!is_bg && col_match)
		 || (is_bg && !col_match)) {
			if (!keep_extra) {
				gzprintf (out, "%s\t%d\t%d\t%d\n", gene, x, y, v);
			} else {
				int n_spaces = 0;
				for (ch=line; *ch!='\0'; ++ch)
        if (isspace(*ch) && (++n_spaces)==3)
          break;
				assert (n_spaces == 3);
				gzprintf (out, "%s\t%d\t%d%s", gene, x, y, ch);
			}

			if (x < min[0]) min[0] = x;
			if (x > max[0]) max[0] = x;
			if (y < min[1]) min[1] = y;
			if (y > max[1]) max[1] = y;
		}
	}
	gzclose (out);

	if (!*ss_image)
		return 0;

	int i, j;
	Mat ss;
	Mat out2;

	max[0] += 1;
	max[1] += 1;
	out2 = Mat::zeros (Size(max[0]-min[0], max[1]-min[1]), CV_8UC3);

	min[0] -= xmin;
	min[1] -= ymin;
	max[0] -= xmin;
	max[1] -= ymin;

	xcv_imread (ss, ss_image);
	assert (img.rows == ss.rows);
	assert (img.cols == ss.cols);

	for (i=0; i<ss.rows; ++i) {
		if (i<min[1] || i>=max[1])
			continue;
		for (j=0; j<ss.cols; ++j) {
			if (j<min[0] || j>=max[0])
				continue;

			Vec3b col = img.at<Vec3b>(i,j);
			col_match = col_func (col);

			if (!is_bg && col_match)
				out2.at<Vec3b>(i-min[1],j-min[0]) = ss.at<Vec3b>(i,j);
			else if (is_bg && !col_match)
				out2.at<Vec3b>(i-min[1],j-min[0]) = ss.at<Vec3b>(i,j);
		}
	}
	sprintf (file, "%s/%s.register.ssDNA.tif", outdir, spl_name);
	imwrite (file, out2);

	return 0;
}
