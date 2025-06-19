/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2022-03-25 09:26:31
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int extract_main (int argc, char * argv[]);
extern int prep_main (int argc, char * argv[]);
extern int part_main (int argc, char * argv[]);
extern int merge_main (int argc, char * argv[]);
extern int bin_extr_main (int argc, char * argv[]);

static int
usage (void)
{
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:   gem_ext  <command>  [options]\n");
	fprintf (stderr, "Version: v3.0\n");
	fprintf (stderr, "Contact: chenxi1@genomics.cn\n");
	fprintf (stderr, "Command: extr    Extract the matched ssDNA and GEM part for ROI\n");
	fprintf (stderr, "         prep    Draw mask image with parts\n");
	fprintf (stderr, "         part    Extract matched ssDNA and GEM parts for ROIs\n");
	fprintf (stderr, "         merge   Merge ssDNAs and GEMs from command <part>\n");
	fprintf (stderr, "         bin_ext Extract GEM data on binned heatmap\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Updates: v2.0 -> v3.0\n");
	fprintf (stderr, "         <extr> and <part>:\n");
	fprintf (stderr, "         1. Add paramter '-p' to process ssDNA\n");
	fprintf (stderr, "         <part>:\n");
	fprintf (stderr, "         2. Output a list for command 'merge'\n");
	fprintf (stderr, "         3. Output matched ssDNA small pictures\n");
	fprintf (stderr, "         <merge>:\n");
	fprintf (stderr, "         4. Merge matched ssDNA pictures\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "         v3.0 -> v3.1\n");
	fprintf (stderr, "         <merge>\n");
	fprintf (stderr, "         1. remove restrictions for input gem files\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "         v3.1 -> v3.2\n");
	fprintf (stderr, "         <merge>\n");
	fprintf (stderr, "         1. fix one bug in row/col_min/max_os init\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "         v3.2 -> v3.3\n");
	fprintf (stderr, "         <extr>\n");
	fprintf (stderr, "         1. fix one bug in col setting\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "         v3.3 -> v3.4\n");
	fprintf (stderr, "         <part>\n");
	fprintf (stderr, "         1. fix one bug while checking x limit in gem file\n");
	fprintf (stderr, "\n");

	return 1;
}

int
main (int argc, char * argv[])
{
	if (argc < 2)
		return usage ();

	int ret;

	if (strcmp(argv[1],"extr") == 0)
		ret = extract_main (argc-1, argv+1);
	else if (strcmp(argv[1],"prep") == 0)
		ret = prep_main (argc-1, argv+1);
	else if (strcmp(argv[1],"part") == 0)
		ret = part_main (argc-1, argv+1);
	else if (strcmp(argv[1],"merge") == 0)
		ret = merge_main (argc-1, argv+1);
	else if (strcmp(argv[1],"bin_ext") == 0)
		ret = bin_extr_main (argc-1, argv+1);
	else {
		fprintf (stderr, "Invalid command: '%s'\n", argv[1]);
		return usage ();
	}

	return ret;
}
