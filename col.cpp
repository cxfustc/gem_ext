/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2022-03-25 11:01:26
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"

using namespace cv;

Vec3b xcol_roi;

int
is_red (Vec3b & col)
{
	if (col[2] < col[0]+20)
		return 0;
	if (col[2] < col[1]+20)
		return 0;

	return 1;
}

int
is_green (Vec3b & col)
{
	if (col[1] < col[0]+20)
		return 0;
	if (col[1] < col[2]+20)
		return 0;

	return 1;
}

int
is_blue (Vec3b & col)
{
	if (col[0] < col[1]+20)
		return 0;
	if (col[0] < col[2]+20)
		return 0;

	return 1;
}

int
is_col (Vec3b & col)
{
	if (col[0] != xcol_roi[0])
		return 0;
	if (col[1] != xcol_roi[1])
		return 0;
	if (col[2] != xcol_roi[2])
		return 0;

	return 1;
}
